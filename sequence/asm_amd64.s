
//func packedKmerAt(data []byte, offset, k int) int32
TEXT ·packedKmerAt(SB),7,$0
  MOVQ data+0(FP), AX //AX points to the slice data
  MOVQ offset+24(FP), BX

  //turn the offset from bases into bytes
  MOVQ BX, CX
  ANDQ $3, CX //keep remainder (in bases) for later shifting
  SHRQ $2, BX //divide by 4 to get index in bytes

  //load 16 bases, potentially up to 3 extra at the beginning
  ADDQ BX, AX
  MOVQ (AX), X1
  //shuffle big-endian to little-endian, just use 128-bit registers :(
  MOVQ $0x7777777700010203, R8
  MOVQ R8, X0
  PSHUFB X0, X1
  MOVL X1, AX

  //shift to the correct base
  SHLQ $1, CX //remaining bases => remaining bits
  SHLL CL, AX //AX left-most bit is now the first of the k-mer 

  //and shift away the excess (32-k*2)
  MOVQ k+32(FP), BX
  SHLQ $1, BX
  MOVQ $32, CX
  SUBQ BX, CX
  SHRL CL, AX

  //return
  MOVQ AX, ret+40(FP)
  RET

//func packBytes(s, data []byte)
TEXT ·packBytes(SB),7,$0
  MOVQ s+0(FP), AX //AX points to the string data
  MOVQ s+8(FP), R8 //input length
  MOVQ data+24(FP), BX //output data
  
  MOVQ $0x7700770177027703, R9 //for shuffling
  MOVQ R9, X0
  MOVQ $0x0004000400040004, R11
  MOVQ $0x0003000300030003, R12

  nextbyte:
  //load 4 zero-padded bytes into 64-bits
  MOVL (AX), X1
  PSHUFB X0, X1
  MOVQ X1, R10
  //copy and shift by 1
  MOVQ X1, R9
  SHRQ $1, R9
  //original "&4" padded, shifted by 2
  ANDQ R11, R10
  SHRQ $2, R10
  //exclusive or, "&3"
  XORQ R10, R9
  ANDQ R12, R9

  //pull out four bytes, shift and OR
  MOVQ $0, DX
  MOVB R9, DL
  SHRQ $14, R9
  ORB R9, DL
  SHRQ $14, R9
  ORB R9, DL
  SHRQ $14, R9
  ORB R9, DL
  //note a total shift of 56 bits, leaving 8 bits in DL

  //write 1 byte of data
  MOVB DL, (BX)

  ADDQ $4, AX
  ADDQ $1, BX
  SUBQ $4, R8
  CMPQ R8, $0
  JGE nextbyte

  RET

//func packedCountKmers(data []byte, upTo, skipFront, skipBack, k int, seeds []bool) int
TEXT ·packedCountKmers(SB),7,$0
  MOVQ data+0(FP), AX //AX points to the slice data
  MOVQ data+8(FP), R8 //data length
  MOVQ upTo+24(FP), SI
  MOVQ skipBack+40(FP), BX
  MOVQ k+48(FP), DX

  SUBQ $1, R8 //ignore first byte
  SHLQ $2, R8 //now length is in bases
  SUBQ BX, R8 //reduce by the missing tail bases
  SUBQ DX, R8 //reduce by k (i.e. bases->kmers)
  ADDQ $1, R8 //k-1, not k
  MOVQ R8, R15 //R15 will hold the remainder
  //reduce to nearest multiple of 4
  ANDQ $0xFFFFFFFFFFFFFFFC, R8
  ANDQ $0x0000000000000003, R15

  //get the right shifts we'll need to isolate k-mers
  SHLQ $1, DX
  MOVQ $32, CX
  SUBQ DX, CX //CX holds the required shift

  //at this point, BX and DX are free to use

  MOVQ $0x7777777700010203, R9 //for big->little endian
  MOVQ R9, X0

  //add k-mers until we hit data length
  MOVQ seeds+56(FP), DX //DX holds the seeds data
  MOVQ $0, R9 //count

  //do the first byte one k-mer at a time
  MOVQ (AX), X1
  PSHUFB X0, X1
  MOVL X1, R10 //R10 is the first byte, starting at pos 32
  MOVQ skipFront+32(FP), BX
  MOVL $4, R11
  SUBL BX, R11 //bases to process at front
  SHLL $1, BX //bits to remove
  MOVL CX, R13 //temp swap
  MOVL BX, CX
  SHLL CL, R10 //R10 now has the correct k-mer at start
  MOVL CX, BX
  MOVL R13, CX //revert swap

  initial:
  MOVL R10, R12
  SHRL CL, R12
  ADDQ DX, R12 //k-mer as index in seeds
  MOVB (R12), R12
  ADDB R12, R9 //count++ (if seed exists)

  SHLL $2, R10 //move to next k-mer
  ADDL $2, BX
  CMPL BX, $6
  JLE initial
  

  internal:
  ADDQ $1, AX 
  MOVQ (AX), X1
  PSHUFB X0, X1

  //four k-mers: trim off 0-3 bases from the front
  MOVL X1, R10
  SHRL CL, R10

  ADDQ DX, R10 //k-mer as index in seeds
  MOVB (R10), R10

  MOVL X1, R11
  SHLL $2, R11
  SHRL CL, R11
  ADDQ DX, R11
  MOVB (R11), R11
  
  MOVL X1, R12
  SHLL $4, R12
  SHRL CL, R12
  ADDQ DX, R12
  MOVB (R12), R12
  
  MOVL X1, R13
  SHLL $6, R13
  SHRL CL, R13
  ADDQ DX, R13
  MOVB (R13), R13

  //increment count
  ADDB R10, R11
  ADDB R12, R13
  ADDB R11, R13
  ANDQ $0x00000000000000FF, R13 //clear the non-byte part
  ADDQ R13, R9 //minimise dependencies on R9

  CMPQ R9, SI
  JGE endtail
  SUBQ $4, R8
  CMPQ R8, $4
  JGE internal

  //the < 4 bases (R15) in the tail
  ADDQ $1, AX 
  MOVQ (AX), X1
  PSHUFB X0, X1
  MOVL X1, R10

  tail:
  CMPQ R15, $0
  JE endtail

  MOVL R10, R12
  SHRL CL, R12
  ADDQ DX, R12
  MOVB (R12), R12
  ANDQ $0x00000000000000FF, R12 
  ADDQ R12, R9

  SHLL $2, R10 //move to next k-mer
  SUBL $1, R15
  JMP tail

  endtail:
  //return
  MOVQ R9, ret+80(FP)

  RET

//func packedWriteSegments(data []byte, skipFront, skipBack, k int, seeds []bool, segments []int)
TEXT ·packedWriteSegments(SB),7,$0
  //largely a duplicate of packedCountKmers, but with more output
  MOVQ data+0(FP), AX //AX points to the slice data
  MOVQ data+8(FP), R8 //data length
  MOVQ skipBack+32(FP), BX
  MOVQ k+40(FP), DX
  MOVQ segments+72(FP), R14 //output

  //setup some constants in XMM registers
  MOVQ DX, R9 
  NEGQ R9
  MOVQ R9, X2 // -k, used to reset offsets

  SUBQ $1, R8 //ignore first byte
  SHLQ $2, R8 //now length is in bases
  SUBQ BX, R8 //reduce by the missing tail bases
  SUBQ DX, R8 //reduce by k (i.e. bases->kmers)
  ADDQ $1, R8 //k-1, not k
  MOVQ R8, R15 //R15 will hold the remainder
  //reduce to nearest multiple of 4
  ANDQ $0xFFFFFFFFFFFFFFFC, R8
  ANDQ $0x0000000000000003, R15

  //get the right shifts we'll need to isolate k-mers
  SHLQ $1, DX
  MOVQ $32, CX
  SUBQ DX, CX //CX holds the required shift

  //at this point, BX and DX are free to use

  MOVQ $0x7777777700010203, R9 //for big->little endian
  MOVQ R9, X0

  //add k-mers until we hit data length
  MOVQ seeds+48(FP), DX //DX holds the seeds data
  MOVQ $0, R9 //the offset

  //do the first byte one k-mer at a time
  MOVQ (AX), X1
  PSHUFB X0, X1
  MOVL X1, R10 //R10 is the first byte, starting at pos 32
  MOVQ skipFront+24(FP), BX
  MOVL $4, R11
  SUBL BX, R11 //bases to process at front
  SHLL $1, BX //bits to remove
  MOVL CX, R13 //temp swap
  MOVL BX, CX
  SHLL CL, R10 //R10 now has the correct k-mer at start
  MOVL CX, BX
  MOVL R13, CX //revert swap

  initial:
  MOVL R10, R12
  SHRL CL, R12
  MOVL R12, R13
  ADDQ DX, R12 //k-mer as index in seeds
  MOVB (R12), R12
  CMPB R12, $0
  JE nohit

  //write offset and k-mer (so save k-mer above)
  MOVQ R9, (R14)
  MOVQ R13, 8(R14) //64-bit values TODO: move to int32?
  ADDQ $16, R14
  //Then offset needs to become -k
  MOVQ X2, R9

  //Then regardless, increment offset and continue
  nohit:
  ADDQ $1, R9 //increment offset
  SHLL $2, R10 //move to next k-mer
  ADDL $2, BX
  CMPL BX, $6
  JLE initial
  

  internal:
  ADDQ $1, AX 
  MOVQ (AX), X1
  PSHUFB X0, X1

  //four k-mers: trim off 0-3 bases from the front
  MOVL X1, R10
  SHRL CL, R10
  MOVL R10, R13

  ADDQ DX, R10 
  MOVB (R10), R10
  CMPB R10, $0
  JE nohita
  //write
  MOVQ R9, (R14) //offset
  MOVQ R13, 8(R14) //seed
  ADDQ $16, R14
  MOVQ X2, R9 //-k
  
  nohita:
  ADDQ $1, R9 //base increment

  MOVL X1, R11
  SHLL $2, R11
  SHRL CL, R11
  MOVQ R11, R13

  ADDQ DX, R11
  MOVB (R11), R11
  CMPB R11, $0
  JE nohitb

  MOVQ R9, (R14) //offset
  MOVQ R13, 8(R14) //seed
  ADDQ $16, R14
  MOVQ X2, R9 //-k

  nohitb:
  ADDQ $1, R9 //base increment
  
  MOVL X1, R12
  SHLL $4, R12
  SHRL CL, R12
  MOVQ R12, R13

  ADDQ DX, R12
  MOVB (R12), R12
  CMPB R12, $0
  JE nohitc
  
  MOVQ R9, (R14) //offset
  MOVQ R13, 8(R14) //seed
  ADDQ $16, R14
  MOVQ X2, R9 //-k

  nohitc:
  ADDQ $1, R9 //base increment

  MOVL X1, R13
  SHLL $6, R13
  SHRL CL, R13
  MOVQ R13, R10
  ADDQ DX, R13
  MOVB (R13), R13
  CMPB R13, $0
  JE nohitd

  MOVQ R9, (R14) //offset
  MOVQ R10, 8(R14) //seed
  ADDQ $16, R14
  MOVQ X2, R9 //-k

  nohitd:
  ADDQ $1, R9 //increment

  SUBQ $4, R8
  CMPQ R8, $4
  JGE internal

  //the < 4 bases (R15) in the tail
  ADDQ $1, AX 
  MOVQ (AX), X1
  PSHUFB X0, X1
  MOVL X1, R10

  tail:
  CMPL R15, $0
  JE endtail

  MOVL R10, R12
  SHRL CL, R12
  MOVQ R12, R13
  ADDQ DX, R12
  MOVB (R12), R12
  CMPB R12, $0
  JE nohitend
  
  MOVQ R9, (R14) //offset
  MOVQ R13, 8(R14) //seed
  ADDQ $16, R14
  MOVQ X2, R9 //-k

  nohitend:
  ADDQ $1, R9 //increment
  SHLL $2, R10 //move to next k-mer
  SUBL $1, R15
  JMP tail

  endtail:

  MOVQ X2, R8
  SUBQ R8, R9  //turn back into bases
  SUBQ $1, R9
  MOVQ R9, (R14) //final offset

  RET
