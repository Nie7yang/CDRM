
      SUBROUTINE KINEM3(KB,NODE,&
                        DT,R1,R2,ALPH,AM,DX,&
                        H1,W1,H2,W2,HB,WB,&
                        GANMA,TOUSUI,ASOU,ALPHA,&
                        GANMAC,TOUSUIC,SLOPE,LENGTH,&
						AREA,KEISU2,WIDTH,LANDUSE,NCLASS,&
						TIME,ORDER,GNUM,TH2B)

      INTEGER KB,NODE
      REAL DT,R1,R2,ALPH,AM,DX
      REAL H1(*),W1(*),H2(*),W2(*),HB(*),WB(*)
      REAL GANMA,TOUSUI,ASOU,ALPHA
      REAL GANMAC,TOUSUIC,SLOPE,LENGTH

	  REAL AREA,WIDTH
	  INTEGER NCLASS,LANDUSE,TIME,ORDER,GNUM
	  DOUBLEPRECISION KEISU2
!
      INTEGER I,J
      REAL DRDT,DELR,DDT,RT,TREM,H,C,DWDX,HM,DWDH,DWDT1,DWDT2,DH1,DH2,TH,TW,D,DS
      REAL DC,MU,BETAC,ALPHC

      REAL WSUB
      REAL TC1A, TC1B, TC1C, TC2A, TC2B 
      REAL TW1A, TW1B, TW1C, TW2A, TW2B 
      REAL TH1A, TH1B, TH2A, TH2B 
      REAL E, TS2B, DDX, TS

      REAL SMALL_VALUE 
      PARAMETER(SMALL_VALUE = 0.001)

      D = ASOU * GANMA
      DC = ASOU * GANMAC
      DS = ASOU * (GANMAC+GANMA)
      BETAC = TOUSUI/TOUSUIC

!     WSUB:  [m3/sec/m]
      WSUB = ((DS-DC)*TOUSUI+DC*TOUSUIC)*SLOPE/(LENGTH*1000.0)
      WSUB = WSUB * AREA * KEISU2 / WIDTH
!
!     (1) DERIVATIVE OF THE INPUT

      IF(  (DT<=0.0) &
			.OR. (KB<1) &
			.OR. (AM<=0.) &
			.OR. (ALPH<=0.)) THEN
!			.OR. (R1<0.) &
!			.OR. (R2<0.) ) THEN
			WRITE(*,*) DT,KB,AM,ALPH,R1,R2
      ENDIF
      IF(  (DT<=0.0) &
			.OR. (KB<1) &
			.OR. (AM<=0.) &
			.OR. (ALPH<=0.)) STOP ' !! Input data for KINEM1 is incorrect !! '
!			.OR. (R1<0.) &
!			.OR. (R2<0.) ) STOP ' !! Input data for KINEM1 is incorrect !! '
      DRDT = (R2-R1)/DT
      DELR = (R2-R1)/REAL(KB)
!
      DDT = DT/REAL(KB)
      TH = H1(1)
      TW = W1(1)
!
      DO I=1,KB
        RT = R1 + (REAL(I)*2.0 - 1.0) * DELR / 2.0
!     (2) TIME INCREMENT H ....Courant偺忦審 
!         TREM... REMAINING TIME
        TREM = DDT
    1   H = TREM
!
        DO J=1,NODE
          CALL HTODQDH(GANMA,TOUSUI,ASOU,GANMAC,TOUSUIC,SLOPE,LENGTH, &
                       ALPH,AM,H1(J),C)
          IF((C>0).AND.(H>DX/C)) H = DX/C
        ENDDO
!
!     (3) PROCEED BY H
        H2(1) = TH + (HB(I)-TH) * (DDT-TREM)/DDT
        W2(1) = TW + (WB(I)-TW) * (DDT-TREM)/DDT
!
        IF(NODE>2) THEN
          DWDX = (W1(2)-W1(1))/DX
          HM = (H1(2)+H1(1))/2.0
          CALL HTODQDH(GANMA,TOUSUI,ASOU,GANMAC,TOUSUIC,SLOPE,LENGTH, &
                         ALPH,AM,HM,DWDH)
          DWDT2 = DWDH * (RT-DWDX)
!
          DO J=2,NODE-1
            DH1 = RT - (W1(J+1)-W1(J-1))/DX/2.0
            DWDT1 = DWDT2
            DWDX = (W1(J+1)-W1(J))/DX
            HM = (H1(J+1)+H1(J))/2.0
            CALL HTODQDH(GANMA,TOUSUI,ASOU,GANMAC,TOUSUIC,SLOPE,LENGTH, &
                         ALPH,AM,HM,DWDH)
            DWDT2 = DWDH * (R2-DWDX)
            DH2 = DRDT - (DWDT2-DWDT1)/DX
            H2(J) = H1(J) + H*DH1 + 0.5*(H**2.0)*DH2
            IF(H2(J)<0.0) H2(J) = 0.0
            CALL HTOQ(GANMA,TOUSUI,ASOU,GANMAC,TOUSUIC,SLOPE,LENGTH, &
                      ALPH,AM,H2(J),W2(J))
          ENDDO
        END IF
        DWDX = (W1(NODE) - W1(NODE-1))/DX
        H2(NODE) = H1(NODE) + H*(RT-DWDX)
        IF(H2(NODE)<0.0) H2(NODE)=0.0
        CALL HTOQ(GANMA,TOUSUI,ASOU,GANMAC,TOUSUIC,SLOPE,LENGTH, &
                      ALPH,AM,H2(NODE),W2(NODE))

      DO J=2,NODE

      TH1A = (H1(J-1) - ASOU * (GANMAC+GANMA)) / 1000.0 ! [ [m]
      TH1B = (H1(J) - ASOU * (GANMAC+GANMA)) / 1000.0
      TH2A = (H2(J-1) - ASOU * (GANMAC+GANMA)) / 1000.0
      TH2B = (H2(J) - ASOU * (GANMAC+GANMA)) / 1000.0
      IF(TH1A<SMALL_VALUE) TH1A = 0.
      IF(TH1B<SMALL_VALUE) TH1B = 0.
      IF(TH2A<SMALL_VALUE) TH2A = 0.
      IF(TH2B<SMALL_VALUE) TH2B = 0.

      TW1A = W1(J-1) * AREA * KEISU2 / WIDTH - WSUB ! [m3/sec/m]
      TW1B = W1(J) * AREA * KEISU2 / WIDTH - WSUB
      TW2A = W2(J-1) * AREA * KEISU2 / WIDTH - WSUB
      TW2B = W2(J) * AREA * KEISU2 / WIDTH - WSUB
      IF(J/=NODE)THEN
			TW1C = W1(J+1) * AREA * KEISU2 / WIDTH - WSUB
      ELSE
			TW1C = W1(J) * AREA * KEISU2 / WIDTH - WSUB
      ENDIF
	  ENDDO

!
!    (4) ALTERNATION OF THE ROLE OF THE WORKING AREA 
        DO J=1,NODE
          H1(J) = H2(J)
          W1(J) = W2(J)
        ENDDO
        TREM = TREM - H
!
!     (5) CHECK OF THE REMAINING TIME 
        IF(TREM>(0.001*DDT)) GOTO 1
!
        TH = HB(I)
        TW = WB(I)
        HB(I) = H2(NODE)
        WB(I) = W2(NODE)
      ENDDO
      RETURN
      END SUBROUTINE KINEM3
