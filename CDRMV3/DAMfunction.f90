
!********************************************************************************
! DAM FUNCTION (Apip, 2011)
!********************************************************************************

	SUBROUTINE DAMfunction (RAIN,T,KEISU2,&
							TF1,St1,Qin1,Qout1,Hdam1,&
							St2,Qout2,H2)

!	KEISU2 for converting from mm/hr to m3/sec or from mm to m3 (KEISU2*3600)     
	real	Qin1,Qout1,St1
	real	initS,initH,initQout,TF1,RAIN,Hdam1
	real	Qin2,dS2,St2,H2,Qout2
	integer T
	doubleprecision KEISU2

		Qin2=TF1*KEISU2											![m3/sec]

!		Equation 1 
!		dS2=(Qin2-Qout1)*3600.									![m3]
!		Equation 2
		dS2=(((Qin2+Qin1)/2.)-Qout1)*3600.						![m3]库容变化量
		
		St2= dS2+St1											![m3]

		H2=(0.00000413*St2)+880.00							![m]
	    
		if(H2<881.29) H2=880.00						!Minimum 汛限水位
		if(H2>881.29) H2=880.00						!Maximum 起调水位

		if (H2<881.29) then  !兴利水位
!			Qout2=(5.9833333316*H2)-2557.2766659386								   ![m3/sec]
            Qout2= 0.0
        elseif(H2<=881.29) then
!			Qout2=(2.3750000000*(H2**2.))-(2023.3650000000*H2)+430941.8100000000   ![m3/sec]
            Qout2= 0.0
        else
!			Qout2=(-1.*0.340679492884281*(H2**2.))+(309.615333627541*H2)-70098.14212149340 ![m3/sec]
            Qout2= 0.0
		endif

		if(Qout2<0.0) Qout2=0.0													![m3]

	return
	end subroutine DAMfunction


	SUBROUTINE DAMfunctionnew (RAIN,T,KEISU2,&
							TF1,St1,Qin1,Qout1,Hdam1,&
							St2,Qout2,H2)

!	KEISU2 for converting from mm/hr to m3/sec or from mm to m3 (KEISU2*3600)     
	real	Qin1,Qout1,St1
	real	initS,initH,initQout,TF1,RAIN,Hdam1
	real	Qin2,dS2,St2,H2,Qout2
	integer T
	doubleprecision KEISU2

		Qin2=TF1*KEISU2											![m3/sec]

!		Equation 1 
!		dS2=(Qin2-Qout1)*3600.									![m3]
!		Equation 2
		dS2=(((Qin2+Qin1)/2.)-Qout1)*3600.						![m3]
		
		St2= dS2+St1											![m3]

		H2=(0.0000002914*St2)+880.00						![m]

		if(H2<880.00) H2=880.00						!Minimum		

		if (H2<880.00) then
!			Qout2=(8.2444444445*H2)-2845.1577777861								   ![m3/sec]
            Qout2= 0
        elseif(H2<=881.29) then
!			Qout2=(2.2205357143*(H2**2.))-(1522.2968214286*H2)+260888.1569999970   ![m3/sec]
            Qout2= 0
        else
!			Qout2=(1.1669696970*(H2**2.))-(782.2661515151*H2)+130936.7934848360	   ![m3/sec]
            Qout2= 0
		endif

		if(Qout2<0.0) Qout2=0.0													![m3]

	return
	end subroutine DAMfunctionnew
