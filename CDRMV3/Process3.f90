
      SUBROUTINE PROCESS3(LANDUSE,NSUC,ZONE,ORDER,&
                          OUTLET,NZONE,TOTAL,MAXORDER,&
                          MAXT,RAINT,KB,NODE,NCLASS,&
                          ALPH,WIDTH,AREA,nobasin,dam,Theta,Costheta,&
                          GRID,DT,MAXRAIN,RINTERVAL,MAXCRAIN,MAXCRETEN,&
                          FILER,FILER2,FILEORAIN,FILEO,namebasin,INTCP,&
		                  FILEWSTORE,RSA,F1,UPPERCELL,QI1,QI2,QI3,areaQI,&
                          GANMA,TOUSUI,ASOU,ALPHA,&
                          GANMAC,TOUSUIC,SLOPE,SLOPE2,LENGTH,&
						  TCOL,TROW,CONV)

      INTEGER	LANDUSE(*),NSUC(*),ZONE(*),ORDER(*),dam
      INTEGER	OUTLET,NZONE,TOTAL,MAXORDER
      INTEGER	MAXT,RAINT,KB,NODE,NCLASS,UPPERCELL(*)
      REAL		ALPH(*),WIDTH(*),AREA(*),areaQI(*),QI1,QI2,QI3
      INTEGER	MAXRAIN,nobasin
      REAL		GRID,DT,RINTERVAL,MAXCRAIN,MAXCRETEN
      INTEGER	FILER,FILER2,FILEORAIN,FILEO,FILEWSTORE
      REAL		GANMA,TOUSUI,ASOU,ALPHA(*)
      REAL		GANMAC,TOUSUIC,SLOPE(*),SLOPE2(*),LENGTH(*),Theta(*),Costheta(*)
!
      INTEGER	I,J,K,T,Z
      INTEGER	CRAIN1,CRAIN2
	  INTEGER	TCOL,TROW,CONV(*),nslide
      REAL		ERAIN,ERAIN2,DX,CRETEN2,TF1,TF2,TFdam1,TFdam2,TFdamin1,DAMIN1,TFdamin2,DAMIN2
      REAL		AM,INPUTW,TMPH,TMPW,RAIN
      REAL		TNODEH1,TNODEW1,TNODEH2,TNODEW2
      REAL		RETENTION,CRAIN,CRETEN
      REAL		RSA,F1,DS,DC,Dtot

	  REAL		storage1,storage2,storage3
      DIMENSION AM(TOTAL),INPUTW(0:TOTAL,KB),TMPH(KB),TMPW(KB)
      DIMENSION RAIN(NZONE)
      REAL, DIMENSION(:,:), ALLOCATABLE :: ENODEH, ENODEW

	  REAL		ERAINet,COST,SINT,TroughI,LeafI
	  REAL		COVURBAN,COVPADDY,COVFIRLD,COVFOREST,COVORCHA,COVWILDS,COVWATER,COVRIVER
	  REAL		IMAXURBAN,IMAXPADDY,IMAXFIRLD,IMAXFOREST,IMAXORCHA,IMAXWILDS,IMAXWATER,IMAXRIVER
	  REAL		PAURBAN,PAPADDY,PAFIRLD,PAFOREST,PAORCHA,PAWILDS,PAWATER,PARIVER,PA,LA,SA
	  REAL		DTRI,DTRIn,SFI,SFIn,LDI,LDIn,GROSSRAIN,GROSSRAINn,IS,ISIn
	  REAL		TotDTRI,TotSFI,TotLDI,TotIS,TOTRAIN
	  REAL		COV,ICmax
!
      DIMENSION TNODEH1(NODE),TNODEW1(NODE),TNODEH2(NODE),TNODEW2(NODE)
      DIMENSION RETENTION(0:MAXRAIN),CRAIN(TOTAL),CRETEN(TOTAL)
      DOUBLE PRECISION TD1,TD2,DH1,DH2,DQ1,DQ2,KEISU,KEISU2
	  INTEGER	YYYY,MM,DD,HH,RGA

!	  CHARACTER namebasin*no    
	  CHARACTER cy*4
	  CHARACTER namebasin*2,simtime*12,INTCP
!
	  REAL		RAINTOT				
	  INTEGER, DIMENSION(:), ALLOCATABLE :: TAS

!	  Added for Interception Process:
	  REAL, DIMENSION(:), ALLOCATABLE :: IC,DirT,TIF,SF,LD,ICstore

!	  For DAM
	  real, dimension(:,:,:), allocatable :: St,Qinf,Qout,Hdam
	  real		initS,initH,initQin,initQout
	  real		initS2,initH2,initQin2,initQout2
	  real		St2,Qout2,Hdam2,Hdam3
	  integer	TOTALANEDAM,TOTALNYUDAM

!	  For 2-D
      REAL, DIMENSION(:,:), ALLOCATABLE    :: TWOD,TWOD2,TWOD3,TWOD4
!	  INTEGER, DIMENSION(:,:), ALLOCATABLE :: TWOD4
      REAL, DIMENSION(:), ALLOCATABLE      :: TH2B,TS2B,EDEPTH,HDEPTH
      INTEGER COUNT,X,Y

	  ALLOCATE (IC(TOTAL),DirT(TOTAL),TIF(TOTAL),SF(TOTAL),LD(TOTAL),ICstore(TOTAL))
      ALLOCATE (TH2B(TOTAL),TS2B(TOTAL),EDEPTH(TOTAL),HDEPTH(TOTAL),TAS(NZONE))
      ALLOCATE (TWOD(TCOL,TROW),TWOD2(TCOL,TROW),TWOD3(TCOL,TROW),TWOD4(TCOL,TROW))
      ALLOCATE (ENODEH(TOTAL,NODE), ENODEW(TOTAL,NODE))
	  ALLOCATE (St(MAXT,TOTAL,KB),Qinf(MAXT,TOTAL,KB),Qout(MAXT,TOTAL,KB),Hdam(MAXT,TOTAL,KB))

!==========================================================================================
	  IF((INTCP=='N').OR.(INTCP=='n')) GOTO 3005
!	  Interception Parameters
!	  Parameters for interception processes (surface coverege -%-): 
      COVURBAN  =0.0    
      COVPADDY  =0.4     
      COVFIRLD  =0.2		!corn 
      COVFOREST =0.9		!
      COVORCHA  =0.5		!coffee    
      COVWILDS  =0.8		!grass    
      COVWATER  =0.0    
      COVRIVER  =0.0   
!	  Min. Interception Storage (mm)
	  IMAXURBAN  =0.0    
	  IMAXPADDY  =0.8     
	  IMAXFIRLD  =0.7		!corn 
	  IMAXFOREST =2.5		!
	  IMAXORCHA  =0.5		!coffee    
	  IMAXWILDS  =2.3		!grass    
	  IMAXWATER  =0.0    
	  IMAXRIVER  =0.0   
!     Max.Interception Storage
!	  IMAXURBAN  =0.0    
!	  IMAXPADDY  =20.8     
!	  IMAXFIRLD  =20.7		!corn 
!	  IMAXFOREST =40.5		!
!	  IMAXORCHA  =30.5		!coffee    
!	  IMAXWILDS  =20.3		!grass    
!	  IMAXWATER  =0.0    
!	  IMAXRIVER  =0.0    
!	  The average of acute angle (degrees)
      PAURBAN  =0.0    
      PAPADDY  =50.0     
      PAFIRLD  =30.0		!corn 
      PAFOREST =10.0		!
      PAORCHA  =15.0		!coffee    
      PAWILDS  =60.0		!grass    
      PAWATER  =0.0    
      PARIVER  =0.0   

!	  Persen Area for stemflow and Direct throughfall 
	  LA=0.8
	  SA=0.2

3005  CONTINUE
 	  RGA=0
	  DO J=1,NZONE
			TAS(J)=0
			DO I=1,TOTAL
				IF (ZONE(I)==J) THEN
				TAS(J)=TAS(J)+1
				ENDIF
			ENDDO
!			WRITE (*,*) TAS(J) !Area/total grid for each rainfall station
			RGA=RGA+TAS(J)
	  ENDDO

!=========================================================================================
3010  CONTINUE
      IF(FILER2/=0) THEN
        RETENTION(0) = 0.0
        DO I=1,MAXRAIN
          READ(FILER2,*) RETENTION(I)
        ENDDO

        TF1 = MAXCRAIN/RINTERVAL

        CRAIN2 = NINT(TF1)
        CRAIN1 = CRAIN2 - 1
        TF1 = RETENTION(CRAIN2)-RETENTION(CRAIN1)
        TF2 = MAXCRAIN-REAL(CRAIN1)*RINTERVAL
        CRETEN2 = TF1*TF2/RINTERVAL+ RETENTION(CRAIN1)
        TF2 = MAXCRETEN/CRETEN2
        DO I=1,MAXRAIN
          RETENTION(I) = RETENTION(I)*TF2
        ENDDO
      END IF
!
      DX = 1.0 / (REAL(NODE) - 1.0)
      DH1 = 0.1/DBLE(GRID)
      DH2 = 1.0/DBLE(GRID)
      DO I=1,TOTAL
        IF(LANDUSE(I)==NCLASS) THEN
          TD1 = (DBLE(WIDTH(I)) * DH1) / (DBLE(WIDTH(I)) + 2.0*DH1)
          TD2 = (DBLE(WIDTH(I)) * DH2) / (DBLE(WIDTH(I)) + 2.0*DH2)
          DQ1 = TD1**0.66666667 * DH1
          DQ2 = TD2**0.66666667 * DH2
          AM(I) = REAL(LOG(DQ1/DQ2) / LOG(DH1/DH2))
        ELSE
          AM(I) = 5.0/3.0
        END IF

        DO K=1,KB
          INPUTW(I,K) = 0.0
        ENDDO
        DO K=1,NODE
          ENODEH(I,K) = 0.0
          ENODEW(I,K) = 0.0
        ENDDO

      ENDDO
      DO K=1,KB
        TMPH(K) = 0.0
        TMPW(K) = 0.0
        INPUTW(0,K) = 0.0
      ENDDO
      DO K=1,NODE
        TNODEH1(K) = 0.0
        TNODEW1(K) = 0.0
        TNODEH2(K) = 0.0
        TNODEW2(K) = 0.0
      ENDDO
      DO I=1,TOTAL
        CRAIN(I) = 0.0
        CRETEN(I) = 0.0
      ENDDO

!****************************************************
      KEISU = DBLE(1.0/DT)
      KEISU2 = DBLE(GRID*GRID /(1000.0 * 3600.0))
!
!****************************************************
!
      WRITE(*,*) '*********************开始处理....检查***************************'

!++++++++++++++++++++++++++++++++++++++++++Setup QI in Each Grid +++++++++++++++++++++++++++++++++++++

	  TOTALANEDAM=0
	  TOTALNYUDAM=0
	  DO J=1,TOTAL
		if(areaQI(J)==1.0) then
			TOTALANEDAM=TOTALANEDAM+1
		elseif(areaQI(J)==2.0) then
			TOTALNYUDAM=TOTALNYUDAM+1
		else
		endif
	  ENDDO

      DO I=1,TOTAL
		DO J=1,NODE
			if(nobasin==1) then
			  if(areaQI(I)==-9999.0) then		!Whole Anegawa Catchment
				ENODEW(I,J) = QI1/KEISU2/AREA(I)*(UPPERCELL(I)-REAL(NODE-J)/REAL(NODE-1))/REAL(TOTAL)	
!				write (*,*) 'ttttt'			
			  elseif(areaQI(I)==1.0) then		!Anegawa DAM Catchment
				ENODEW(I,J) = QI2/KEISU2/AREA(I)*(UPPERCELL(I)-REAL(NODE-J)/REAL(NODE-1))/REAL(TOTALANEDAM)
			  elseif (areaQI(I)==2.0) then	!Nyu DAM Catchment
				ENODEW(I,J) = QI3/KEISU2/AREA(I)*(UPPERCELL(I)-REAL(NODE-J)/REAL(NODE-1))/REAL(TOTALNYUDAM)
			  else	
				ENODEW(I,J) = QI1/KEISU2/AREA(I)*(UPPERCELL(I)-REAL(NODE-J)/REAL(NODE-1))/REAL(TOTAL)							  							
			  endif
			elseif(nobasin==2) then			!Anegawa DAM Catchment
				ENODEW(I,J) = QI2/KEISU2/AREA(I)*(UPPERCELL(I)-REAL(NODE-J)/REAL(NODE-1))/REAL(TOTAL)
			else								!Nyu DAM Catchment
				ENODEW(I,J) = QI3/KEISU2/AREA(I)*(UPPERCELL(I)-REAL(NODE-J)/REAL(NODE-1))/REAL(TOTAL)			
			endif

				CALL QTOH(GANMA,TOUSUI,ASOU,GANMAC,TOUSUIC,SLOPE(I),LENGTH(I), &
                  ALPH(I),AM(I),ENODEW(I,J),ENODEH(I,J))
		ENDDO
      ENDDO
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      WRITE(*,*) '**************开始玩命运算中**************'

      write(FILEORAIN, *) '年',',','月',',','日',',','时',',','Gross',',','Netto',',','Throughfall',',','LeafDrainage',',','Stemflow',',','StoreInterception'
      write(FILEWSTORE,*) '年',',','月',',','日',',','时',',','RAINTOT',',','Total(m3)',',','Saturated(m3)',',','Unsaturated(m3)'
!      write(FILEO,'(4A5,24A15)') 'YYYY','MM','DD','HH','RAINTOT','H','QoutBsn','HDam1','Dam1in','Dam1out','HDam2','Dam2in','Dam2out'
      write(FILEO,*) '年',',','月',',','日',',','时',',','降雨量',',','H',',','流量',',','HDam1',',','Dam1in',',','Dam1out',',','HDam2',',','Dam2in',',','Dam2out'  !nie
!     Soil Depth
	  DC=ASOU*GANMAC			!Unsaturated Depth Max.
	  DS=ASOU*GANMA				!Saturated Depth Max.[mm]
	  Dtot=ASOU*(GANMAC+GANMA)	!Total Depth Max. (Unsaturated+Saturated Depth) [mm]

!+++++++++++++++++++++++++++++++++++++++  水库初始条件  +++++++++++++++++++++++++++++++++++++++++
!	  Anegawa水库初始条件
	  initQin=QI2
!	  initH=427.57   !event1
	  initH=880.00   !event2
	  																						![m]
	  if(initH<880.00) initH=880.00			
	  InitS=(initH-880.000000)/0.00000413													![m3]
	  
	  if (initH<880.0) then
!		initQout=(5.9833333316*initH)-5007.35										![m3/sec]
        initQout=0
	  elseif(initH<=881.29) then
!		initQout=(2.3750000000*(initH**2.))-(5007.35*initH)+430941.8100000000		![m3/sec]
        initQout=0
	  else
!		initQout=(-1.*0.340679492884281*(initH**2.))+(309.615333627541*initH)-70098.14212149340	![m3/sec]
       initQout=0 
	  endif
	  if(initQout<0.0) initQout=0.0

!	  Niew水库初始条件
	  initQin2=QI3
!	  initH2=345.10  !event1
	  initH2=880.00  !event2
	  																						![m]
	  if(initH2<880.00) initH2=880.00			
	  	InitS2=(initH2-880.00)/0.0000002914													![m3]
	  
	  if (initH2<881.29) then
!		initQout2=(8.2444444445*initH2)-2845.1577777861										![m3/sec]
        initQout2=0
	  elseif(initH2<=881.29) then
!		initQout2=(2.2205357143*(initH2**2.))-(1522.2968214286*initH2)+260888.1569999970 	![m3/sec]
        initQout2=0
	  else
!		initQout2=(1.1669696970*(initH2**2.))-(782.2661515151*initH2)+130936.7934848360		![m3/sec]
        initQout2=0
        
	  endif
	  if(initQout2<0.0) initQout2=0.0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO T=1,MAXT
		IF((INTCP=='N').OR.(INTCP=='n')) GOTO 3020
		DTRIn=0.0
		TotDTRI=0.0
		SFIn=0.0
		TotSFI=0.0
		LDIn=0.0
		TotLDI=0.0
		GROSSRAIN=0.0
		TOTRAIN=0.0
		IS=0.0
		TotIS=0.0
		
		nslide=0
		Nunstable=0.
		TotFS=0.0

		3020  CONTINUE
        WRITE(*,'(3x,A5,I8,A2,I7)') ' time',T,'/',MAXT
        IF(T<=RAINT) THEN
          READ(FILER,*) YYYY,MM,DD,HH,(RAIN(Z), Z=1,NZONE)
        ELSE
          DO Z=1,NZONE
            RAIN(Z) = 0.0
          ENDDO
        END IF
        DO J=1,TOTAL
          DO K=1,KB
            INPUTW(J,K) = 0.0
          ENDDO
        ENDDO

		RAINTOT=0.0

        DO I=MAXORDER,1, -1
          DO J=1,TOTAL
            IF(ORDER(J)==I) THEN
              DO K=1,KB
                TMPW(K) = INPUTW(J,K) / AREA(J)
                CALL QTOH(GANMA,TOUSUI,ASOU,GANMAC,TOUSUIC,SLOPE(J),LENGTH(J), &
                          ALPH(J),AM(J),TMPW(K),TMPH(K))
              ENDDO
			  DO K=1,NODE
                TNODEH1(K) = ENODEH(J,K)
                TNODEW1(K) = ENODEW(J,K)
              ENDDO
              CRAIN(J) = CRAIN(J) + RAIN(ZONE(J))
			  IF(FILER2/=0) THEN
                TF1 = CRAIN(J) / RINTERVAL
                CRAIN2 = NINT(TF1)
                CRAIN1 = CRAIN2 - 1
                TF1 = RETENTION(CRAIN2)-RETENTION(CRAIN1)
                TF2 = CRAIN(J)-REAL(CRAIN1)*RINTERVAL
                CRETEN2 = TF1*TF2/RINTERVAL+ RETENTION(CRAIN1)
                ERAIN = RAIN(ZONE(J)) - (CRETEN2 - CRETEN(J))
                CRETEN(J) = CRETEN2
              ELSEIF(FILER2==0) THEN
				
			    ERAIN = RAIN(ZONE(J))
				IF(CRAIN(J) <= RSA) THEN
				    ERAIN = F1 * ERAIN
					ERAIN2=ERAIN

				IF((INTCP=='N').OR.(INTCP=='n')) GOTO 3030
!				Interception Processes are included:
				   IF(RAIN(ZONE(J))>0.0) THEN
						IF(LANDUSE(J)==1) COV = COVWILDS
						IF(LANDUSE(J)==2) COV = COVFOREST
						IF(LANDUSE(J)==3) COV = COVURBAN
						IF(LANDUSE(J)==4) COV = COVURBAN
						IF(LANDUSE(J)==5) COV = COVFIRLD
						IF(LANDUSE(J)==6) COV = COVPADDY
						IF(LANDUSE(J)==7) COV = COVORCHA
						IF(LANDUSE(J)==8) COV = COVWATER
						IF(LANDUSE(J)==9) COV = COVRIVER

						IF(LANDUSE(J)==1) ICmax = IMAXWILDS
						IF(LANDUSE(J)==2) ICmax = IMAXFOREST
						IF(LANDUSE(J)==3) ICmax = IMAXURBAN
						IF(LANDUSE(J)==4) ICmax = IMAXURBAN
						IF(LANDUSE(J)==5) ICmax = IMAXFIRLD
						IF(LANDUSE(J)==6) ICmax = IMAXPADDY
						IF(LANDUSE(J)==7) ICmax = IMAXORCHA
						IF(LANDUSE(J)==8) ICmax = IMAXWATER
						IF(LANDUSE(J)==9) ICmax = IMAXRIVER

						IF(LANDUSE(J)==1) PA = PAWILDS
						IF(LANDUSE(J)==2) PA = PAFOREST
						IF(LANDUSE(J)==3) PA = PAURBAN
						IF(LANDUSE(J)==4) PA = PAURBAN
						IF(LANDUSE(J)==5) PA = PAFIRLD
						IF(LANDUSE(J)==6) PA = PAPADDY
						IF(LANDUSE(J)==7) PA = PAORCHA
						IF(LANDUSE(J)==8) PA = PAWATER
						IF(LANDUSE(J)==9) PA = PARIVER

						IC(J)=RAIN(ZONE(J))*COV			!intercepted rainfall
						DirT(J)=RAIN(ZONE(J))-IC(J)		!Direct throughfall
					
						ICstore(J)=ICmax*(1-EXP(-1.0*RAIN(ZONE(J))/ICmax))    !HYLUC Model

						IF(IC(J)>ICstore(J)) THEN
							TIF(J)=IC(J)-ICstore(J)
						ELSE
							ICstore(J)=IC(J)
							TIF(J)=0.0
						ENDIF
					
						IF(TIF(J)>0.0) THEN
							IF(LANDUSE(J)== 1) THEN    !For Grass
								PA=PA*3.14/180
								COST=COS(PA)
								SINT=SIN(PA)
								SF(J)=0.25*TIF(J)*(COST*SINT*SINT)       !25% as SF 
								LD(J)=TIF(J)-SF(J)
							ELSE
								PA=PA*3.14/180
								COST=COS(PA)
								SF(J)=0.25*TIF(J)*COST
								LD(J)=TIF(J)-SF(J)
							ENDIF
						ELSE
								SF(J)=0.0
								LD(J)=0.0
						ENDIF

!						TOTAL NET RAINFALL:
						ERAIN=DirT(J)+SF(J)+LD(J)

					ELSE
						IC(J)=0.0
						DirT(J)=0.0
						ICstore(J)=0.0
						TIF(J)=0.0					
						SF(J)=0.0
						LD(J)=0.0
						ERAIN=0.0
					ENDIF
3030		    CONTINUE		
				ENDIF
			  ENDIF
			  
              RAINTOT=RAINTOT+((1./(TOTAL-1))*ERAIN)  !areal rainfall
              ERAIN = ERAIN * KEISU / AREA(J)
!
			  IF((INTCP=='N').OR.(INTCP=='n')) GOTO 3040
!			  Total Rainfall for each Component:
!					Gross Rainfall
					GROSSRAIN=GROSSRAIN+(RAIN(ZONE(J))*(1./(TOTAL-1)))						
!					Direct Through Rainfall
					DTRIn=DTRIn+(DirT(J)*(1./(TOTAL-1)))
!					Leafdrainage
					LDIn=LDIn+(LD(J)*(1./(TOTAL-1)))
!					Stemflow
					SFIn=SFIn+(SF(J)*(1./(TOTAL-1)))
!					Interception store
					IS=IS+(ICstore(J)*(1./(TOTAL-1)))

					ERAINet = ERAIN !* KEISU / AREA(J)
					TroughI = DirT(J) * KEISU / AREA(J)
					LeafI = LD(J) * KEISU / AREA(J)
!
3040		  CONTINUE
              CALL KINEM3(KB,NODE,&
                          DT,ERAIN,ERAIN,ALPH(J),AM(J),DX,&
                          TNODEH1,TNODEW1,TNODEH2,TNODEW2,TMPH,TMPW,&
                          GANMA,TOUSUI,ASOU,ALPHA(J),&
                          GANMAC,TOUSUIC,SLOPE(J),LENGTH(J),&
						  AREA(J),KEISU2,WIDTH(J),LANDUSE(J),NCLASS,&
						  T,I,J,TH2B(J))
!
!	         ************************  水库计算  ************************			  		 
			  if(dam==1) then
					do K=1,KB
						if (nobasin==1.and.order(j)==4681.or.nobasin==1.and.order(j)==18185) then	!整个区域, Anegawa DAM and New DAM	
									
							if (nobasin==1.and.order(j)==4681) then				!Ane DAM					
								TF1 = TMPW(K)*AREA(J)
								if(T==1) then
									call DAMfunction(ERAIN2,T,KEISU2,&
										TF1,initS,initQin,initQout,initH,&
										St2,Qout2,Hdam2)
								else
									call DAMfunction(ERAIN2,T,KEISU2,&
										TF1,St(T-1,J,K),Qinf(T-1,J,K),Qout(T-1,J,K),Hdam(T-1,J,K),&
										St2,Qout2,Hdam2)
								endif
								Qinf(T,J,K)=TF1*KEISU2			![m3/sec]
								St(T,J,K) =St2 					![m3]
								Hdam(T,J,K)=Hdam2
								Qout(T,J,K)=Qout2				![m3/sec]
								TF2=Qout2/KEISU2
								INPUTW(NSUC(J),K) = INPUTW(NSUC(J),K) + TF2

!								水库入流与出流
								if(order(j)==4681) then
									DAMIN1=INPUTW(j,K)
									DAMOUT1=TF2	
								else					
								endif
							else														!Nyu DAM
								TF1 = TMPW(K)*AREA(J)
!								write (*,*) nobasin,landuse(j),order(j),nsuc(j),outlet				  !for debug
!								stop
								if(T==1) then
									call DAMfunctionnew(ERAIN2,T,KEISU2,&
										TF1,initS2,initQin2,initQout2,initH2,&  
										St2,Qout2,Hdam3)
								else
									call DAMfunctionnew(ERAIN2,T,KEISU2,&
										TF1,St(T-1,J,K),Qinf(T-1,J,K),Qout(T-1,J,K),Hdam(T-1,J,K),&  
										St2,Qout2,Hdam3)
								endif
								Qinf(T,J,K)=TF1*KEISU2			![m3/sec]
								St(T,J,K) =St2 					![m3]
								Hdam(T,J,K)=Hdam2
								Qout(T,J,K)=Qout2				![m3/sec]
								TF2=Qout2/KEISU2
								INPUTW(NSUC(J),K) = INPUTW(NSUC(J),K) + TF2

!								DAM Inflow and Outflow
								if(order(j)==18185) then
								DAMIN2=INPUTW(j,K)
								DAMOUT2=TF2						
								else
								endif						
							endif

						elseif (nobasin==2.and.landuse(j)==nclass.and.nsuc(j)==outlet) then !Anegawa DAM Catchment			
							TF1 = TMPW(K)*AREA(J)
!							write (*,*) nobasin,landuse(j),order(j),nsuc(j),outlet				  !for debug
!							stop
							if(T==1) then
								call DAMfunction(ERAIN2,T,KEISU2,&
									TF1,initS,initQin,initQout,initH,&  
									St2,Qout2,Hdam2)
							else
								call DAMfunction(ERAIN2,T,KEISU2,&
									TF1,St(T-1,J,K),Qinf(T-1,J,K),Qout(T-1,J,K),Hdam(T-1,J,K),&  
									St2,Qout2,Hdam2)
							endif
							Qinf(T,J,K)=TF1*KEISU2			![m3/sec]
							St(T,J,K) =St2 					![m3]
							Hdam(T,J,K)=Hdam2
							Qout(T,J,K)=Qout2				![m3/sec]
							TF2=Qout2/KEISU2
							INPUTW(NSUC(J),K) = INPUTW(NSUC(J),K) + TF2

!							DAM Inflow and Outflow
							if(landuse(j)==nclass.and.nsuc(j)==outlet) then
								DAMIN1=INPUTW(j,K)
								DAMOUT1=TF2						
							else
								DAMIN2=INPUTW(j,K)
								DAMOUT2=TF2
							endif						
						elseif (nobasin==3.and.landuse(j)==nclass.and.nsuc(j)==outlet) then !Niew DAM Catchment			
							TF1 = TMPW(K)*AREA(J)
!							write (*,*) nobasin,landuse(j),order(j),nsuc(j),outlet				  !for debug
!							stop
							if(T==1) then
								call DAMfunctionnew(ERAIN2,T,KEISU2,&
									TF1,initS2,initQin2,initQout2,initH2,&  
									St2,Qout2,Hdam3)
							else
								call DAMfunctionnew(ERAIN2,T,KEISU2,&
									TF1,St(T-1,J,K),Qinf(T-1,J,K),Qout(T-1,J,K),Hdam(T-1,J,K),&  
									St2,Qout2,Hdam3)
							endif
							Qinf(T,J,K)=TF1*KEISU2			![m3/sec]
							St(T,J,K) =St2 					![m3]
							Hdam(T,J,K)=Hdam2
							Qout(T,J,K)=Qout2				![m3/sec]
							TF2=Qout2/KEISU2
							INPUTW(NSUC(J),K) = INPUTW(NSUC(J),K) + TF2

!							DAM流入和流出
							if(landuse(j)==nclass.and.nsuc(j)==outlet) then
								DAMIN2=INPUTW(j,K)
								DAMOUT2=TF2						
							else
								DAMIN1=INPUTW(j,K)
								DAMOUT1=TF2
							endif						

						else
							TF1 = TMPW(K)*AREA(J)
							INPUTW(NSUC(J),K) = INPUTW(NSUC(J),K) + TF1
						endif
					enddo
			  else
					do K=1,KB
						TF1 = TMPW(K)*AREA(J)
						INPUTW(NSUC(J),K) = INPUTW(NSUC(J),K) + TF1
					enddo
			  endif              
			  
			  DO K=1,NODE
                ENODEH(J,K) = TNODEH2(K)
                ENODEW(J,K) = TNODEW2(K)
              ENDDO
            ENDIF
          ENDDO
        ENDDO

!       每个网格的平均水深
		DO J=1,TOTAL
		HDEPTH(J)=0.0
		DO K=1,NODE
			HDEPTH(J)=HDEPTH(J)+(ENODEH(J,K)/NODE)
		ENDDO
		ENDDO

!		WATER STORAGE
		STORAGE1 = 0.0
		STORAGE2 = 0.0
		STORAGE3 = 0.0
		do J=1,TOTAL
		do K=1,NODE
!			在不饱和层			
			if(ENODEH(J,K)<dc) then
			STORAGE1 = STORAGE1 + (ENODEH(J,K)*KEISU2*3600.0/NODE)
			else
			STORAGE1 = STORAGE1 + (dc*KEISU2*3600.0/NODE)
			endif
!			在饱和层
			if(ENODEH(J,K)<dc) then
			STORAGE2 = STORAGE2 + 0.0
			elseif(ENODEH(J,K)<dtot.and.ENODEH(J,K)>=dc) then
			STORAGE2 = STORAGE2 + (ENODEH(J,K)*KEISU2*3600.0/NODE)
			else
			STORAGE2 = STORAGE2 + (dtot*KEISU2*3600.0/NODE)
			endif
!			total
			STORAGE3 = STORAGE3 + (ENODEH(J,K)*KEISU2*3600.0/NODE)
		end do
		end do
!________________________________________________________________________________________________
!      流域出口处的流量:
	    if(dam==1) then
			TF1 = INPUTW(OUTLET, KB)
!			Dam流入和流出: 	
			TFdamin1=DAMIN1
			TFdamin2=DAMIN2
			TFdam1 = DAMOUT1
			TFdam2 = DAMOUT2
	    else
			TF1 = INPUTW(OUTLET, KB)
			TFdamin1=0.0
			TFdamin2=0.0
			TFdam1 = 0.0
			TFdam2 = 0.0
	    endif

!      河道流量 
!       write(FILEO,'(4I5,F15.3,F15.3,8F16.3)') YYYY,MM,DD,HH,RAINTOT,TF1/(TOTAL-1),TF1*KEISU2,Hdam2,TFdamin1*KEISU2,TFdam1*KEISU2,&
!			  Hdam3,TFdamin2*KEISU2,TFdam2*KEISU2
        write(FILEO,'(I5,A3,I5,A3,I5,A3,I5,A3,F15.3,A3,F15.3,A3,F16.3,A3,F16.3,A3,F16.3,A3,F16.3,A3,F16.3,A3,F16.3,A3,F16.3)') YYYY,',',MM,',',DD,',',HH,',',RAINTOT,',',TF1/(TOTAL-1),',',TF1*KEISU2,',',Hdam2,',',TFdamin1*KEISU2,',',TFdam1*KEISU2,',',Hdam3,',',TFdamin2*KEISU2,',',TFdam2*KEISU2  !nie
!	   流域储水 
       write(FILEWSTORE,'(I5,A3,I5,A3,I5,A3,I5,A3,F17.5,A3,F17.5,A3,F17.5,A3,F17.5)') YYYY,',',MM,',',DD,',',HH,',',RAINTOT,',',storage3,',',storage2,',',storage1

!      截取
	   IF((INTCP=='N').OR.(INTCP=='n')) GOTO 3060
       WRITE(FILEORAIN, '(I5,A3,I5,A3,I5,A3,I5,A3, F16.5,A3,F16.5,A3,F16.5,A3,F16.5,A3,F16.5,A3,F16.5)') YYYY,',',MM,',',DD,',',HH,',',GROSSRAIN,',',RAINTOT,',',DTRIn,',',LDIn,',',SFIn,',',IS

3060   CONTINUE
!________________________________________________________________________________________________
!   河流等级（River Channel Order）平面图
	if(T==1) then
    open(MAXT+201,FILE='./输出数据/河道.asc',STATUS='Unknown')
      COUNT = 0
      J = 0
      do Y = 1, TROW
			do X = 1, TCOL
				COUNT = COUNT + 1
				TWOD4(X,Y) = -9999
				if(CONV(COUNT)>=0) then
					J = J + 1
					if (landuse(j)==nclass) then
					TWOD4(X,Y)=order(J)	
					endif			
				endif
			enddo
      enddo

	  WRITE(MAXT+201,'(2x,A,2X,A,2X,A,2X,I5)') 'X', 'Y', 'Code' 		
      DO Y = 1, TROW
	  DO X = 1, TCOL
				WRITE(MAXT+201,'(2I5,1X,F10.3)') X,TROW-Y+1,TWOD4(X,Y)
	  ENDDO
      ENDDO

!	  WRITE(MAXT+200,'(A)') 'ncols         539            '
!	  WRITE(MAXT+200,'(A)') 'nrows         744            '
!	  WRITE(MAXT+200,'(A)') 'xllcorner     602630.84259341'
!	  WRITE(MAXT+200,'(A)') 'yllcorner     3914823.670891 '
!	  WRITE(MAXT+200,'(A)') 'cellsize      50             '
!	  WRITE(MAXT+200,'(A)') 'NODATA_value  -9999.00       '

!     do Y = 1,TROW
!        write(MAXT+200,'(<TCOL>F8.2)') (TWOD4(X,Y), X = 1,TCOL)
!     enddo
	  endif
!________________________________________________________________________________________________
!   River Channel Network 2-D 
	 IF(T==1) THEN
     open(MAXT+300,FILE='./输出数据/流域编码.asc',STATUS='Unknown')
      COUNT = 0
      J = 0
      do Y = 1, TROW
			do X = 1, TCOL
				COUNT = COUNT + 1
				TWOD4(X,Y) = -9999
				if(CONV(COUNT)>=0) then
					J = J + 1
					if (landuse(j)==nclass) then
					TWOD4(X,Y)=UPPERCELL(J)	
					endif			
				endif
			enddo
      enddo

	  WRITE(MAXT+300,'(2x,A,2X,A,2X,A,2X,I5)') 'X', 'Y', 'Depth' !,T-1 			
      DO Y = 1, TROW
	  DO X = 1, TCOL
				WRITE(MAXT+300,'(2I5,1X,F10.3)') X,TROW-Y+1,TWOD4(X,Y)
	  ENDDO
      ENDDO

!	  WRITE(MAXT+300,'(A)') 'ncols         539            '
!	  WRITE(MAXT+300,'(A)') 'nrows         744            '
!	  WRITE(MAXT+300,'(A)') 'xllcorner     602630.84259341'
!	  WRITE(MAXT+300,'(A)') 'yllcorner     3914823.670891 '
!	  WRITE(MAXT+300,'(A)') 'cellsize      50             '
!	  WRITE(MAXT+300,'(A)') 'NODATA_value  -9999.00       '
!     do Y = 1,TROW
!        write(MAXT+300,'(<TCOL>F8.2)') (TWOD4(X,Y), X = 1,TCOL)
!     enddo
	  endif
!________________________________________________________________________________________________
!   输出水深分布图 [mm] 
	WRITE(cy,'(I4.4)') T
!	IF(T==MAXT) THEN
    OPEN(MAXT+1., FILE = './输出数据/2D水深分布/水深状态_'//cy//'.ASC',STATUS='Unknown')
      COUNT = 0
      J = 0
      DO Y = 1, TROW
			DO X = 1, TCOL
				COUNT = COUNT + 1
				TWOD(X,Y) = -9999.
				IF(CONV(COUNT) >=0) THEN
!                    write(*,*) conv(count),count
					J = J + 1
					IF(LANDUSE(J)/=NCLASS) THEN          !没有河流
!                        if(j<=TOTAL)then !nie2018.2.1HDEPTH经常莫名越界所以加入限制条件
					        TWOD(X,Y) = HDEPTH(J) 
!					        TWOD(X,Y) = TH2B(J)*1000.		     !Overland
!                       endif
					ELSE
					    TWOD(X,Y) = 0.0					 
                        
					ENDIF
				IF(J==OUTLET) TWOD(X,Y) = 0.0
				ENDIF
			ENDDO
      ENDDO

!	  WRITE(MAXT+1.,'(2x,A,2X,A,2X,A,2X,I5)') 'X', 'Y', 'Depth' !,T-1 			
!     DO Y = 1, TROW
!	  DO X = 1, TCOL
!			IF (TWOD(X,Y)>=0.0) THEN
!				WRITE(MAXT+1.,'(2I5,1X,F10.3)') X,TROW-Y+1, TWOD(X,Y)
!			ENDIF
!	  ENDDO
!     ENDDO
!##############################手动输入原来GIS转ASCLL的开头#########################
	  WRITE(MAXT+1,'(A)') 'ncols         539            '                        !##
	  WRITE(MAXT+1,'(A)') 'nrows         744            '                        !##
	  WRITE(MAXT+1,'(A)') 'xllcorner     602630.84259341'                        !##
	  WRITE(MAXT+1,'(A)') 'yllcorner     3914823.670891 '                        !##
	  WRITE(MAXT+1,'(A)') 'cellsize      50             '                        !##
	  WRITE(MAXT+1,'(A)') 'NODATA_value  -9999.00       '                        !##
!###################################################################################
      DO Y = 1,TROW
        WRITE(MAXT+1,'(<TCOL>F10.2)') (TWOD(X,Y), X = 1,TCOL)
      ENDDO
!	  ENDIF
!________________________________________________________________________________________________

      ENDDO

      DEALLOCATE (ENODEH)    !dellocate释放动态数组内存
      DEALLOCATE (ENODEW)
!     DEALLOCATE (HDEPTH)   !nie

      RETURN
      END SUBROUTINE PROCESS3
