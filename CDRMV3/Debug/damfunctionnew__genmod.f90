        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 22 01:24:42 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DAMFUNCTIONNEW__genmod
          INTERFACE 
            SUBROUTINE DAMFUNCTIONNEW(RAIN,T,KEISU2,TF1,ST1,QIN1,QOUT1, &
     &HDAM1,ST2,QOUT2,H2)
              REAL(KIND=4) :: RAIN
              INTEGER(KIND=4) :: T
              REAL(KIND=8) :: KEISU2
              REAL(KIND=4) :: TF1
              REAL(KIND=4) :: ST1
              REAL(KIND=4) :: QIN1
              REAL(KIND=4) :: QOUT1
              REAL(KIND=4) :: HDAM1
              REAL(KIND=4) :: ST2
              REAL(KIND=4) :: QOUT2
              REAL(KIND=4) :: H2
            END SUBROUTINE DAMFUNCTIONNEW
          END INTERFACE 
        END MODULE DAMFUNCTIONNEW__genmod
