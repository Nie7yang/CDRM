        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 22 01:24:41 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE KINEM3__genmod
          INTERFACE 
            SUBROUTINE KINEM3(KB,NODE,DT,R1,R2,ALPH,AM,DX,H1,W1,H2,W2,HB&
     &,WB,GANMA,TOUSUI,ASOU,ALPHA,GANMAC,TOUSUIC,SLOPE,LENGTH,AREA,     &
     &KEISU2,WIDTH,LANDUSE,NCLASS,TIME,ORDER,GNUM,TH2B)
              INTEGER(KIND=4) :: KB
              INTEGER(KIND=4) :: NODE
              REAL(KIND=4) :: DT
              REAL(KIND=4) :: R1
              REAL(KIND=4) :: R2
              REAL(KIND=4) :: ALPH
              REAL(KIND=4) :: AM
              REAL(KIND=4) :: DX
              REAL(KIND=4) :: H1(*)
              REAL(KIND=4) :: W1(*)
              REAL(KIND=4) :: H2(*)
              REAL(KIND=4) :: W2(*)
              REAL(KIND=4) :: HB(*)
              REAL(KIND=4) :: WB(*)
              REAL(KIND=4) :: GANMA
              REAL(KIND=4) :: TOUSUI
              REAL(KIND=4) :: ASOU
              REAL(KIND=4) :: ALPHA
              REAL(KIND=4) :: GANMAC
              REAL(KIND=4) :: TOUSUIC
              REAL(KIND=4) :: SLOPE
              REAL(KIND=4) :: LENGTH
              REAL(KIND=4) :: AREA
              REAL(KIND=8) :: KEISU2
              REAL(KIND=4) :: WIDTH
              INTEGER(KIND=4) :: LANDUSE
              INTEGER(KIND=4) :: NCLASS
              INTEGER(KIND=4) :: TIME
              INTEGER(KIND=4) :: ORDER
              INTEGER(KIND=4) :: GNUM
              REAL(KIND=4) :: TH2B
            END SUBROUTINE KINEM3
          END INTERFACE 
        END MODULE KINEM3__genmod
