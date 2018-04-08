
#define POPCNT_R9 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xC9
#define POPCNT_R11 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xDB
#define POPCNT_R13 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xED
#define POPCNT_R15 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xFF
#define POPCNT_R8 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xC0
#define POPCNT_R10 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xD2
#define POPCNT_R12 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xE4
#define POPCNT_R14 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xF6

#define SETZ_DL BYTE $0x0f; BYTE $0x94; BYTE $0xC2

//func countIntersectionToAsm(a, b []uint64, maxCount int) uint
TEXT ·countIntersectionToAsm(SB),7,$0
  MOVQ a+0(FP), AX //AX is a's data
  MOVQ b+24(FP), BX //BX is b's data
  MOVQ a+8(FP), CX //length
  MOVQ maxCount+48(FP), SI
  XORQ DX, DX //count

  docount:
  
  CMPQ CX, $7
  JLE tail
  CMPQ DX, SI
  JGE finished

  //do eight 64-bit entries per loop
  /*MOVQ (AX), R8
  MOVQ (BX), R9
  MOVQ 8(AX), R10
  MOVQ 8(BX), R11
  MOVQ 16(AX), R12
  MOVQ 16(BX), R13
  MOVQ 24(AX), R14
  MOVQ 24(BX), R15
  ANDQ R8, R9
  ANDQ R10, R11
  ANDQ R12, R13
  ANDQ R14, R15*/
  LDDQU (AX), X0
  LDDQU (BX), X1
  LDDQU 16(AX), X2
  LDDQU 16(BX), X3
  LDDQU 32(AX), X4
  LDDQU 32(BX), X5
  LDDQU 48(AX), X6
  LDDQU 48(BX), X7
  PAND X0, X1
  PAND X2, X3
  PAND X4, X5
  PAND X6, X7

  MOVQ X1, R9
  MOVQ X3, R13
  MOVQ X5, R8
  MOVQ X7, R12
  //apparently PEXTRQ is slow, so shuffle and extract instead
  /*
  PEXTRQ $1, X1, R11
  PEXTRQ $1, X3, R15
  PEXTRQ $1, X5, R10
  PEXTRQ $1, X7, R14
  */
  MOVHLPS X1,X1
  MOVHLPS X3,X3
  MOVHLPS X5,X5
  MOVHLPS X7,X7
  MOVQ X1, R11
  MOVQ X3, R15
  MOVQ X5, R10
  MOVQ X7, R14  


  POPCNT_R9
  POPCNT_R11
  POPCNT_R13
  POPCNT_R15
  POPCNT_R8
  POPCNT_R10
  POPCNT_R12
  POPCNT_R14
  //low dependency additions. Does this actuall make a difference?
  ADDQ R9, R11
  ADDQ R13, R15
  ADDQ R8, R10
  ADDQ R12, R14
  ADDQ R11, R15
  ADDQ R10, R14
  ADDQ R15, DX
  ADDQ R14, DX

  ADDQ $64, AX
  ADDQ $64, BX
  SUBQ $8, CX
  JMP docount

  tail:
  //move one at a time for the final entries

  CMPQ CX, $0
  JLE finished
  
  MOVQ (AX), R8
  MOVQ (BX), R9
  ANDQ R8, R9
  POPCNT_R9
  ADDQ R9, DX

  ADDQ $8, AX
  ADDQ $8, BX
  DECQ CX
  JMP tail

  finished:
  MOVQ DX, ret+56(FP)
  RET


//func getSoftUnion4Asm(vs []uint64) (uint64, uint64, uint64, uint64)
TEXT ·getSoftUnion4Asm(SB),7,$0
  MOVQ vs+0(FP), AX //vs
  MOVQ vs+8(FP), BX //n

  XORQ R8, R8
  XORQ R9, R9
  XORQ R10, R10
  XORQ R11, R11

  CMPQ BX, $3 //don't unroll when there are fewer than 4 values
  JLE nextround

  //Use the same unrolling as for 16s

  //1.
  MOVQ (AX), R8 //v1
  //2.
  MOVQ 8(AX), DX
  MOVQ R8, R9
  ANDQ DX, R9 //v2
  ORQ DX, R8  //v1
  //3.
  MOVQ 16(AX), DX
  MOVQ R9, R10
  ANDQ DX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2
  ORQ DX, R8 //v1
  //4.
  MOVQ 24(AX), DX
  MOVQ R10, R11
  ANDQ DX, R11 //v4
  MOVQ R9, CX
  ANDQ DX, CX
  ORQ CX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2
  ORQ DX, R8 //v1

  //cater for the unrolled part
  SUBQ $4, BX
  ADDQ $32, AX 

  nextround:
  CMPQ BX, $0
  JLE finished

  // shuffle up
  MOVQ (AX), DX 
  MOVQ R10, CX 
  ANDQ DX, CX
  ORQ CX, R11 //v4 updated
  MOVQ R9, CX 
  ANDQ DX, CX
  ORQ CX, R10
  MOVQ R8, CX 
  ANDQ DX, CX
  ORQ CX, R9
  ORQ DX, R8

  DECQ BX
  ADDQ $8, AX

  JMP nextround

  finished:
  MOVQ R8, ret+24(FP)
  MOVQ R9, ret+32(FP)
  MOVQ R10, ret+40(FP)
  MOVQ R11, ret+48(FP)
  RET

//func getSoftUnion8Asm(vs []uint64) (uint64, uint64, uint64, uint64)
TEXT ·getSoftUnion8Asm(SB),7,$0
  MOVQ vs+0(FP), AX //vs
  MOVQ vs+8(FP), BX //n

  //clear registers used for output (in case n is low)
  XORQ R12, R12
  XORQ R13, R13
  XORQ R14, R14
  XORQ R15, R15
  
  CMPQ BX, $5 //don't unroll when there are fewer than 6 values (the number we unroll)
  JLE nextround

  //we can unroll a bunch of the first iterations (which require fewer ops)

  //1.
  MOVQ (AX), R8 //v1
  //2.
  MOVQ 8(AX), DX
  MOVQ R8, R9
  ANDQ DX, R9 //v2
  ORQ DX, R8  //v1
  //3.
  MOVQ 16(AX), DX
  MOVQ R9, R10
  ANDQ DX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2
  ORQ DX, R8 //v1
  //4.
  MOVQ 24(AX), DX
  MOVQ R10, R11
  ANDQ DX, R11 //v4
  MOVQ R9, CX
  ANDQ DX, CX
  ORQ CX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2
  ORQ DX, R8 //v1
  //5.
  MOVQ 32(AX), DX
  MOVQ R11, R12
  ANDQ DX, R12 //v5
  MOVQ R10, CX
  ANDQ DX, CX
  ORQ CX, R11 //v4
  MOVQ R9, CX
  ANDQ DX, CX
  ORQ CX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2
  ORQ DX, R8 //v1
  //6.
  MOVQ 40(AX), DX
  MOVQ R12, R13
  ANDQ DX, R13 //v6
  MOVQ R11, CX
  ANDQ DX, CX
  ORQ CX, R12 //v5
  MOVQ R10, CX
  ANDQ DX, CX
  ORQ CX, R11 //v4
  MOVQ R9, CX
  ANDQ DX, CX
  ORQ CX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2
  ORQ DX, R8 //v1

  //just clear the other two
  XORQ R14, R14
  XORQ R15, R15

  SUBQ $6, BX
  ADDQ $48, AX 

  nextround:
  CMPQ BX, $0
  JLE finished

  MOVQ (AX), DX
  MOVQ R14, CX 
  ANDQ DX, CX
  ORQ CX, R15 //v8 updated
  MOVQ R13, CX
  ANDQ DX, CX
  ORQ CX, R14 //v7
  MOVQ R12, CX
  ANDQ DX, CX
  ORQ CX, R13 //v6
  MOVQ R11, CX
  ANDQ DX, CX
  ORQ CX, R12 //v5
  MOVQ R10, CX
  ANDQ DX, CX
  ORQ CX, R11 //v4
  MOVQ R9, CX
  ANDQ DX, CX
  ORQ CX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2
  ORQ DX, R8 //v1

  ADDQ $8, AX
  DECQ BX

  JMP nextround

  finished:
  MOVQ R12, ret+24(FP) //v5
  MOVQ R13, ret+32(FP)
  MOVQ R14, ret+40(FP)
  MOVQ R15, ret+48(FP) //v8
  RET
 
//func getSoftUnion16Asm(vs []uint64) (uint64, uint64, uint64, uint64)
TEXT ·getSoftUnion16Asm(SB),7,$0
  MOVQ vs+0(FP), AX //vs
  MOVQ vs+8(FP), BX //n

  //zero the xmm registers
  PXOR X0, X0
  PXOR X1, X1
  PXOR X2, X2
  PXOR X3, X3

  //we can unroll a bunch of the first iterations (which require fewer ops)
  //unlike the ones above, the unroll should never go past n

  //1.
  MOVQ (AX), R8 //v1
  //2.
  MOVQ 8(AX), DX
  MOVQ R8, R9
  ANDQ DX, R9 //v2
  ORQ DX, R8  //v1
  //3.
  MOVQ 16(AX), DX
  MOVQ R9, R10
  ANDQ DX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2
  ORQ DX, R8 //v1
  //4.
  MOVQ 24(AX), DX
  MOVQ R10, R11
  ANDQ DX, R11 //v4
  MOVQ R9, CX
  ANDQ DX, CX
  ORQ CX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2
  ORQ DX, R8 //v1
  //5.
  MOVQ 32(AX), DX
  MOVQ R11, R12
  ANDQ DX, R12 //v5
  MOVQ R10, CX
  ANDQ DX, CX
  ORQ CX, R11 //v4
  MOVQ R9, CX
  ANDQ DX, CX
  ORQ CX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2
  ORQ DX, R8 //v1
  //6.
  MOVQ 40(AX), DX
  MOVQ R12, R13
  ANDQ DX, R13 //v6
  MOVQ R11, CX
  ANDQ DX, CX
  ORQ CX, R12 //v5
  MOVQ R10, CX
  ANDQ DX, CX
  ORQ CX, R11 //v4
  MOVQ R9, CX
  ANDQ DX, CX
  ORQ CX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2
  ORQ DX, R8 //v1
  //7.
  MOVQ 48(AX), DX
  MOVQ R13, R14
  ANDQ DX, R14 //v7
  MOVQ R12, CX
  ANDQ DX, CX
  ORQ CX, R13 //v6
  MOVQ R11, CX
  ANDQ DX, CX
  ORQ CX, R12 //v5
  MOVQ R10, CX
  ANDQ DX, CX
  ORQ CX, R11 //v4
  MOVQ R9, CX
  ANDQ DX, CX
  ORQ CX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2
  ORQ DX, R8 //v1
  //8.
  MOVQ 56(AX), DX
  MOVQ R14, R15
  ANDQ DX, R15 //v8
  MOVQ R13, CX
  ANDQ DX, CX
  ORQ CX, R14 //v7
  MOVQ R12, CX
  ANDQ DX, CX
  ORQ CX, R13 //v6
  MOVQ R11, CX
  ANDQ DX, CX
  ORQ CX, R12 //v5
  MOVQ R10, CX
  ANDQ DX, CX
  ORQ CX, R11 //v4
  MOVQ R9, CX
  ANDQ DX, CX
  ORQ CX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2

  SUBQ $8, BX
  ADDQ $64, AX 

  //X0 holds 12,16
  //X1 holds 11,15
  //X2 holds 10,14
  //X3 holds 9,13

  nextround:
  CMPQ BX, $0
  JLE finished
  
  //set:
  //X4 = X1 & m
  //X5 = X2 & m
  //X6 = X3 & m
  // using X7 as m, then
  //X7 = 8,12 & m <- this one can't be SIMD, ran out of registers

  MOVOA X1, X4
  MOVOA X2, X5
  MOVOA X3, X6
  PEXTRQ $1, X0, CX // we need this original v12 value for X7 later
  
  MOVQ (AX), DX
  MOVQ DX, X7
  MOVLHPS X7, X7
  PAND X7, X4
  PAND X7, X5
  PAND X7, X6
  //now use the values in X4-X6, then come back to X7 at the end
  POR X4, X0 //finish up (v12 |= v11 & m) and (v16 |= v15 & m)
  POR X5, X1 //etc
  POR X6, X2
  MOVQ CX, X4 // and v12. Then we'll combine into X7...
  PINSRQ $1, R15, X4 //v8 into X4
  PAND X4, X7
  // and finally combine with X3
  POR X7, X3

  //now do the normal shuffle up the other registers
  MOVQ R14, CX 
  ANDQ DX, CX
  ORQ CX, R15 //v8 updated
  MOVQ R13, CX
  ANDQ DX, CX
  ORQ CX, R14 //v7
  MOVQ R12, CX
  ANDQ DX, CX
  ORQ CX, R13 //v6
  MOVQ R11, CX
  ANDQ DX, CX
  ORQ CX, R12 //v5
  MOVQ R10, CX
  ANDQ DX, CX
  ORQ CX, R11 //v4
  MOVQ R9, CX
  ANDQ DX, CX
  ORQ CX, R10 //v3
  MOVQ R8, CX
  ANDQ DX, CX
  ORQ CX, R9 //v2
  ORQ DX, R8 //v1

  ADDQ $8, AX
  DECQ BX

  JMP nextround

  finished:
  //MOVQ X1, ret+32(FP)
  //PEXTRQ $1, X1, ret+24(FP)
  //MOVQ X0, ret+48(FP)
  //PEXTRQ $1, X0, ret+40(FP)
  MOVQ X3, ret+24(FP) //v13
  MOVQ X2, ret+32(FP) //v14
  MOVQ X1, ret+40(FP) //v15
  MOVQ X0, ret+48(FP) //v16

  RET

