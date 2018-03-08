
#define POPCNT_R9 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xC9
#define POPCNT_R11 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xDB
#define POPCNT_R13 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xED
#define POPCNT_R15 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xFF
#define POPCNT_R8 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xC0
#define POPCNT_R10 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xD2
#define POPCNT_R12 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xE4
#define POPCNT_R14 BYTE $0xF3; BYTE $0x4D; BYTE $0x0F; BYTE $0xB8; BYTE $0xF6

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
  PEXTRQ $0, X1, R9
  PEXTRQ $1, X1, R11
  PEXTRQ $0, X3, R13
  PEXTRQ $1, X3, R15
  PEXTRQ $0, X5, R8
  PEXTRQ $1, X5, R10
  PEXTRQ $0, X7, R12
  PEXTRQ $1, X7, R14

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

/*
//func getSoftUnion4Asm(vs []uint64) (uint64, uint64, uint64, uint64, int)
TEXT ·getSoftUnion4Asm(SB),7,$0
  MOVQ vs+0(FP), AX //vs

  //Use the same unrolling as for 16s

  //1.
  MOVQ (AX), R8 //v1
  //2.
  MOVQ 8(AX), R15
  MOVQ R8, R9
  ANDQ R15, R9 //v2
  ORQ R15, R8  //v1
  //3.
  MOVQ 16(AX), R15
  MOVQ R9, R10
  ANDQ R15, R10 //v3
  MOVQ R8, R14
  ANDQ R15, R14
  ORQ R14, R9 //v2
  ORQ R15, R8 //v1
  //4.
  MOVQ 24(AX), R15
  MOVQ R10, R11
  ANDQ R15, R11 //v4
  MOVQ R9, R14
  ANDQ R15, R14
  ORQ R14, R10 //v3
  MOVQ R8, R14
  ANDQ R15, R14
  ORQ R14, R9 //v2
  ORQ R15, R8 //v1

  //count the number of zeroes in the first 4 values
  XORQ CX, CX // zero count
  TEST R11, R11
  SETZ DL
  ADDQ DX, CX
  TEST R10, R10
  SETZ DL
  ADDQ DX, CX
  TEST R9, R9
  SETZ DL
  ADDQ DX, CX
  TEST R8, R8
  SETZ DL
  ADDQ DX, CX

  MOVQ vs+8(FP), BX //n
  MOVQ BX, DX
  SUBQ $4, DX //max zeroes
  //cater for the unrolled part
  SUBQ $4, BX
  ADDQ $32, AX 

  nextround:
  CMPQ BX, $0
  JLE finished

  MOVQ (AX), R8 //m
  MOVQ R8, R9 //prior v & m will go here

  finished:
  MOVQ R8, ret+24(FP)
  MOVQ R9, ret+32(FP)
  MOVQ R10, ret+40(FP)
  MOVQ R11, ret+48(FP)
  MOVQ CX, ret+56(FP)
  RET

/*
//func getSoftUnion16Asm(vs []uint64) (uint64, int)
TEXT ·getSoftUnion16Asm(SB),7,$0
  MOVQ vs+0(FP), AX //vs

  //we can unroll a bunch of the first iterations (which require fewer ops)

  //1.
  MOVQ (AX), R8 //v1
  //2.
  MOVQ 8(AX), R15
  MOVQ R8, R9
  ANDQ R15, R9 //v2
  ORQ R15, R8  //v1
  //3.
  MOVQ 16(AX), R15
  MOVQ R9, R10
  ANDQ R15, R10 //v3
  MOVQ R8, R14
  ANDQ R15, R14
  ORQ R14, R9 //v2
  ORQ R15, R8 //v1
  //4.
  MOVQ 24(AX), R15
  MOVQ R10, R11
  ANDQ R15, R11 //v4
  MOVQ R9, R14
  ANDQ R15, R14
  ORQ R14, R10 //v3
  MOVQ R8, R14
  ANDQ R15, R14
  ORQ R14, R9 //v2
  ORQ R15, R8 //v1
  //5.
  MOVQ 32(AX), R15
  MOVQ R11, R12
  ANDQ R15, R12 //v5
  MOVQ R10, R14
  ANDQ R15, R14
  ORQ R14, R11 //v4
  MOVQ R9, R14
  ANDQ R15, R14
  ORQ R14, R10 //v3
  MOVQ R8, R14
  ANDQ R15, R14
  ORQ R14, R9 //v2
  ORQ R15, R8 //v1
  //6.
  MOVQ 40(AX), R15
  MOVQ R12, R13
  ANDQ R15, R13 //v6
  MOVQ R11, R14
  ANDQ R15, R14
  ORQ R14, R12 //v5
  MOVQ R10, R14
  ANDQ R15, R14
  ORQ R14, R11 //v4
  MOVQ R9, R14
  ANDQ R15, R14
  ORQ R14, R10 //v3
  MOVQ R8, R14
  ANDQ R15, R14
  ORQ R14, R9 //v2
  ORQ R15, R8 //v1

  //count the number of zeroes in the first 6 values
  XORQ CX, CX // zero count
  TEST R13, R13
  SETZ DL
  ADDQ DX, CX
  TEST R12, R12
  SETZ DL
  ADDQ DX, CX
  TEST R11, R11
  SETZ DL
  ADDQ DX, CX
  TEST R10, R10
  SETZ DL
  ADDQ DX, CX
  TEST R9, R9
  SETZ DL
  ADDQ DX, CX
  TEST R8, R8
  SETZ DL
  ADDQ DX, CX

  MOVQ vs+8(FP), BX //n
  MOVQ BX, DX
  SUBQ $16, DX //max zeroes
  SUBQ $6, BX
  ADDQ $48, AX 

  nextround:
  MOVQ (AX), R8 //m
  MOVQ R8, R9 //prior v & m will go here

  //X0 holds 10,14
  //X1 holds 9,13
  //X2 holds 8,12
  //X3 holds 7,11
  //set:
  //X4 = X1 & m
  //X5 = X2 & m
  //X6 = X3 & m
  //X7 = 6,10 & m <- this one can't be SIMD, ran out of registers
  //alg:
  //X4 to X7 get & m
  // X0 = X0 | X4
  // X1 = X1 | X5 etc.

  //shuffle the first 6 as before

*/
