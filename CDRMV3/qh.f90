!****************************************************************************
!**                                                                        **
!**                   分布网格径流模型的q-h关系子程序                      **
!**                                                                        **
!****************************************************************************
!
      SUBROUTINE HTOQ(GANMA,TOUSUI,ASOU,GANMAC,TOUSUIC,SLOPE,LENGTH, &
                      ALPH,AM,H,Q)
      REAL GANMA,TOUSUI,ASOU,GANMAC,TOUSUIC,SLOPE,LENGTH,ALPH,AM,H,Q
      REAL DC,DS,BETAC
!
      DC = ASOU*GANMAC
      DS = ASOU*(GANMAC+GANMA)
      IF(TOUSUIC/=0) BETAC = TOUSUI/TOUSUIC
!
      IF(H<DC)THEN
			Q = DC*TOUSUIC*((H/DC)**BETAC)*SLOPE/(LENGTH*1000.0)
      ELSEIF(H<DS)THEN
			Q = ((H-DC)*TOUSUI+DC*TOUSUIC)*SLOPE/(LENGTH*1000.0)
      ELSE
			Q = ALPH*(H-DS)**(AM)+((H-DC)*TOUSUI+DC*TOUSUIC)*SLOPE/(LENGTH*1000.0)
      ENDIF
      RETURN
      END SUBROUTINE HTOQ
!
      SUBROUTINE HTODQDH(GANMA,TOUSUI,ASOU,GANMAC,TOUSUIC,SLOPE,LENGTH, &
                         ALPH,AM,H,DQDH)
      REAL GANMA,TOUSUI,ASOU,GANMAC,TOUSUIC,SLOPE,LENGTH,ALPH,AM,H,DQDH
      REAL DC,DS,BETAC
!
      DC = ASOU*GANMAC
      DS = ASOU*(GANMAC+GANMA)
      IF(TOUSUIC/=0) BETAC = TOUSUI/TOUSUIC
!
      IF(H<DC)THEN
			DQDH = BETAC*TOUSUIC*SLOPE*(H/DC)**(BETAC-1.0)/(LENGTH*1000.0)
      ELSEIF(H<DS)THEN
			DQDH = TOUSUI*SLOPE/(LENGTH*1000.0)
      ELSE
			DQDH = AM*ALPH*(H-DS)**(AM-1.0)+TOUSUI*SLOPE/(LENGTH*1000.0)
      ENDIF
      RETURN
      END SUBROUTINE HTODQDH
!
      SUBROUTINE QTOH(GANMA,TOUSUI,ASOU,GANMAC,TOUSUIC,SLOPE,LENGTH, &
                      ALPH,AM,Q,H)
      REAL GANMA,TOUSUI,ASOU,GANMAC,TOUSUIC,SLOPE,LENGTH,ALPH,AM,Q,H
      REAL DC,DS,BETAC,QC,QA,H1,Q1,ERROR
      INTEGER COUNT,MAXCOUNT
      PARAMETER (ERROR=0.0001)
      PARAMETER (MAXCOUNT=500000)
!
      DC = ASOU*GANMAC
      DS = ASOU*(GANMAC+GANMA)
      IF(TOUSUIC/=0) BETAC = TOUSUI/TOUSUIC
!
      QC = DC*TOUSUIC*SLOPE/(LENGTH*1000.0)
      QA = ((DS-DC)*TOUSUI+DC*TOUSUIC)*SLOPE/(LENGTH*1000.0)
!
      IF(QC<=ERROR.AND.QA<=ERROR) THEN  ! modified ( for surface flow without A layer )
			H = (Q/ALPH)**(1.0/AM)
      ELSEIF(Q<QC)THEN
			H = (Q*(LENGTH*1000.0)/(DC*TOUSUIC*SLOPE))**(1.0/BETAC)*DC   ! modified (from AM to BETAC)
      ELSEIF(Q<=QA)THEN  ! modified by T.Sayama (from LT. to LE.)
			H = (Q-DC*TOUSUIC*SLOPE/(LENGTH*1000.0))/(TOUSUI*SLOPE)*(LENGTH*1000.0)+DC
      ELSE
!			H1 = ((Q-(DC*TOUSUIC+(DS-DC)*TOUSUI)*SLOPE/(LENGTH*1000.0))/ALPH)**(1.0/AM)+DS  ! modified by T.Sayama
			H1 = ((Q*(LENGTH*1000.0)-(DC*TOUSUIC+(DS-DC)*TOUSUI)*SLOPE)/ALPH)**(1.0/AM)+DS
			Q1 = ALPH*(H1-DS)**(AM)+((H1-DC)*TOUSUI+DC*TOUSUIC)*SLOPE/(LENGTH*1000.0)
			COUNT = 0
      DO WHILE((ABS(Q1-Q)/Q)>ERROR.AND.ABS(Q1-Q)>ERROR)
			H1 = H1-(Q1-Q)/(AM*ALPH*(H1-DS)**(AM-1.0)+TOUSUI*SLOPE/(LENGTH*1000.0))
			IF(H1<=DS)WRITE(*,*)'WARNING H1<=DS IN QTOH'
					Q1 = ALPH*(H1-DS)**(AM)+((H1-DC)*TOUSUI+DC*TOUSUIC)*SLOPE/(LENGTH*1000.0)
					COUNT = COUNT+1
			IF(COUNT>MAXCOUNT)THEN
					WRITE(*,*)'WARNING IN QTOH'
					WRITE(*,*) Q
					STOP
			ENDIF
      ENDDO
			H = H1
      ENDIF
      END SUBROUTINE QTOH
