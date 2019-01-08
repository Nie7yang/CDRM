        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 22 01:24:41 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PROCESS3__genmod
          INTERFACE 
            SUBROUTINE PROCESS3(LANDUSE,NSUC,ZONE,ORDER,OUTLET,NZONE,   &
     &TOTAL,MAXORDER,MAXT,RAINT,KB,NODE,NCLASS,ALPH,WIDTH,AREA,NOBASIN, &
     &DAM,THETA,COSTHETA,GRID,DT,MAXRAIN,RINTERVAL,MAXCRAIN,MAXCRETEN,  &
     &FILER,FILER2,FILEORAIN,FILEO,NAMEBASIN,INTCP,FILEWSTORE,RSA,F1,   &
     &UPPERCELL,QI1,QI2,QI3,AREAQI,GANMA,TOUSUI,ASOU,ALPHA,GANMAC,      &
     &TOUSUIC,SLOPE,SLOPE2,LENGTH,TCOL,TROW,CONV)
              INTEGER(KIND=4) :: MAXRAIN
              INTEGER(KIND=4) :: NODE
              INTEGER(KIND=4) :: KB
              INTEGER(KIND=4) :: TOTAL
              INTEGER(KIND=4) :: NZONE
              INTEGER(KIND=4) :: LANDUSE(*)
              INTEGER(KIND=4) :: NSUC(*)
              INTEGER(KIND=4) :: ZONE(*)
              INTEGER(KIND=4) :: ORDER(*)
              INTEGER(KIND=4) :: OUTLET
              INTEGER(KIND=4) :: MAXORDER
              INTEGER(KIND=4) :: MAXT
              INTEGER(KIND=4) :: RAINT
              INTEGER(KIND=4) :: NCLASS
              REAL(KIND=4) :: ALPH(*)
              REAL(KIND=4) :: WIDTH(*)
              REAL(KIND=4) :: AREA(*)
              INTEGER(KIND=4) :: NOBASIN
              INTEGER(KIND=4) :: DAM
              REAL(KIND=4) :: THETA(*)
              REAL(KIND=4) :: COSTHETA(*)
              REAL(KIND=4) :: GRID
              REAL(KIND=4) :: DT
              REAL(KIND=4) :: RINTERVAL
              REAL(KIND=4) :: MAXCRAIN
              REAL(KIND=4) :: MAXCRETEN
              INTEGER(KIND=4) :: FILER
              INTEGER(KIND=4) :: FILER2
              INTEGER(KIND=4) :: FILEORAIN
              INTEGER(KIND=4) :: FILEO
              CHARACTER(LEN=2) :: NAMEBASIN
              CHARACTER(LEN=1) :: INTCP
              INTEGER(KIND=4) :: FILEWSTORE
              REAL(KIND=4) :: RSA
              REAL(KIND=4) :: F1
              INTEGER(KIND=4) :: UPPERCELL(*)
              REAL(KIND=4) :: QI1
              REAL(KIND=4) :: QI2
              REAL(KIND=4) :: QI3
              REAL(KIND=4) :: AREAQI(*)
              REAL(KIND=4) :: GANMA
              REAL(KIND=4) :: TOUSUI
              REAL(KIND=4) :: ASOU
              REAL(KIND=4) :: ALPHA(*)
              REAL(KIND=4) :: GANMAC
              REAL(KIND=4) :: TOUSUIC
              REAL(KIND=4) :: SLOPE(*)
              REAL(KIND=4) :: SLOPE2(*)
              REAL(KIND=4) :: LENGTH(*)
              INTEGER(KIND=4) :: TCOL
              INTEGER(KIND=4) :: TROW
              INTEGER(KIND=4) :: CONV(*)
            END SUBROUTINE PROCESS3
          END INTERFACE 
        END MODULE PROCESS3__genmod
