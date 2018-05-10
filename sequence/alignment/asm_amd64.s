
#define PMINUW_BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3A;
#define PMAXUW_BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3E;
#define PMINUW_X0_X1 BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3A; BYTE $0xC8
#define PMINUW_X0_X2 BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3A; BYTE $0xC2
#define PMINUW_X0_X3 BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3A; BYTE $0xD8
#define PMINUW_X0_X4 BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3A; BYTE $0xE0
#define PMINUW_X0_X5 BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3A; BYTE $0xE8
#define PMINUW_X0_X7 BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3A; BYTE $0xF8
#define PMAXUW_X1_X2 BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3E; BYTE $0xD1
#define PMAXUW_X3_X2 BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3E; BYTE $0xD3
#define PMAXUW_X4_X2 BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3E; BYTE $0xD4
#define PMAXUW_X5_X2 BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3E; BYTE $0xD5
#define PMINUW_X7_X2 BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3A; BYTE $0xD7
#define PMINUW_X2_X6 BYTE $0x66; BYTE $0x0F; BYTE $0x38; BYTE $0x3A; BYTE $0xF2
  
//func updateOffsetsAsm(ds, poffs, offsets []uint16, threshold uint16) uint16
TEXT ·updateOffsetsAsm(SB),7,$0
  MOVW $0xFFFF, AX //a mask for later
  //have R9,R10 point to the slice data for distance and previous offset respectively
  MOVQ ds+0(FP), R9
  MOVQ poffs+24(FP), R10
  //load the previous offset data, compare for step vs stay
  LDDQU (R10), X1 
  LDDQU 2(R10), X0
  PMINUW_X0_X1 //missing operator, defined above
  //now X1 holds the minimum of step/stay at each position
  //next are the skips. We need distance information
  LDDQU (R9), X2
  LDDQU (R10), X0        //load the previous offsets again. 
  PADDUSW X2, X0         //this time add the extra distance for the skipped base
  PSLLDQ $2, X0          //shift these to line up the skips
  PINSRW $0, AX, X0     //and maximise the missing (overhanging) offset value
  PMINUW_X0_X1
  //finally the 2-skips. We add another distance and shift again.
  PADDUSW X2, X0         //note that these adds saturate, so we can avoid overflow checks
  PSLLDQ $2, X0          //another shift for 2-skips
  PINSRW $0, AX, X0          
  PMINUW_X0_X1
  //add the original distances
  PADDUSW X2, X1

  //now repeat for the internal 128-bit chunks using their own registers, leaving X1
  LDDQU 16(R10), X3 //X3 in place of X1 above, holding the step values
  LDDQU 18(R10), X0
  PMINUW_X0_X3
  LDDQU 14(R9), X2
  LDDQU 14(R10), X0 //we don't need to shift now - just load the earlier position
  PADDUSW X2, X0
  PMINUW_X0_X3
  LDDQU 12(R10), X0 //but manually add two distances for the larger skip
  PADDUSW X2, X0
  LDDQU 12(R9), X2
  PADDUSW X2, X0
  PMINUW_X0_X3
  LDDQU 16(R9), X2
  PADDUSW X2, X3
  //second 128-bit chunk, now using X4 in place of X1
  LDDQU 32(R10), X4
  LDDQU 34(R10), X0
  PMINUW_X0_X4
  LDDQU 30(R9), X2
  LDDQU 30(R10), X0
  PADDUSW X2, X0
  PMINUW_X0_X4
  LDDQU 28(R10), X0 //but manually add two distances for the larger skip
  PADDUSW X2, X0
  LDDQU 28(R9), X2
  PADDUSW X2, X0
  PMINUW_X0_X4
  LDDQU 32(R9), X2
  PADDUSW X2, X4

  //and the final 128-bits which have limited skips
  LDDQU 48(R10), X5
  MOVOA X5, X0
  PSRLDQ $2, X0
  PINSRW $7, AX, X0          
  PMINUW_X0_X5
  //then the skips checks
  LDDQU 46(R9), X2
  LDDQU 46(R10), X0
  PADDUSW X2, X0
  PMINUW_X0_X5
  LDDQU 44(R10), X0 //but manually add two distances for the larger skip
  PADDUSW X2, X0
  LDDQU 44(R9), X2
  PADDUSW X2, X0
  PMINUW_X0_X5
  LDDQU 48(R9), X2
  PADDUSW X2, X5

  //find the minimum across all values
  PHMINPOSUW X1, X0
  PHMINPOSUW X3, X2
  PHMINPOSUW X4, X6
  PHMINPOSUW X5, X7
  //and take the minimum of the four, redundant packed version
  PMINUW_X0_X7
  PMINUW_X7_X2
  PMINUW_X2_X6
  //write the return value
  PEXTRW $0, X6, AX
  MOVW AX, ret+80(FP)

  //subtract the minimum from all values
  MOVQ $0x0100010001000100, R8 //a mask for permuting (copying) values, 64-bits
  MOVQ R8, X0
  PINSRQ $1, R8, X0 //and the upper 64-bits of X0 are set too
  PSHUFB X0, X6 //now X6 is populated with 16-bit minimum values
  PSUBUSW X6, X1
  PSUBUSW X6, X3
  PSUBUSW X6, X4
  PSUBUSW X6, X5

  //threshold
  MOVW threshold+72(FP), AX
  PINSRW $0, AX, X7
  PSHUFB X0, X7 //populate X7 with the threshold value
  MOVO X7, X2 //now in X2
  PMAXUW_X1_X2 //X2 contains threshold or above in each value
  PCMPEQW X1, X2 //all 1s and 0s in X2 now, with <=threshold values being 1s
  POR X2, X1 //apply the maximum to X1
  //next register.. X7=thresholds, X6=1s
  MOVO X7, X2 
  PMAXUW_X3_X2 
  PCMPEQW X3, X2
  ORPS X2, X3 //apply the maximum to X1
  MOVO X7, X2 
  PMAXUW_X4_X2  //X2 has the max of value/threshold
  PCMPEQW X4, X2 //those equal to the threshold are lower
  ORPS X2, X4 //or them
  MOVO X7, X2
  PMAXUW_X5_X2
  PCMPEQW X5, X2
  ORPS X2, X5

  //write to output
  MOVQ offsets+48(FP), R11
  MOVO X1, 0(R11)
  MOVO X3, 16(R11)
  MOVO X4, 32(R11)
  MOVO X5, 48(R11)

  //TODO: 
  // - maximise any above the lowest cost threshold
  // - find the index of the zero value?

  RET
