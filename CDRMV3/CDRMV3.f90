!****************************************************************************!
!*                                                                          *!
!*    Modified Cell Distributed Rainfall-Runoff Model (CDRM V3)				*!
!*                                                                          *!
!*	  原始模型来自于：児島利治(Kojima Toshiharu)博士  京都大学防灾研究所    *!
!*	  2003.02.02 修改程序 立川康人(Tachikawa Yasuto)            			*!
!*    2003.04.09 加入地下径流,引入初始条件 佐山敬洋(Sayama Takahiro)        *!
!*    2011.01.11 加入截留,平均降雨量,流域信息,存储和坝模块 Apip             *!
!*    2018.06.01 优化结构与输出 聂启阳                                            *!
!*            				                                                *!    
!****************************************************************************!


      PROGRAM	CDRMV3
      INTEGER	COL,ROW,MAXT,NCLASS,KB,NODE
      INTEGER	TOTAL,TOTAL2,OUTLET,OUTLET2,RAINT
      REAL		GRID,DT,TOTAL2_AREA
      REAL		SLP,LengthSLP,SLPav,LengthSLPav,Basinarea,Max,Min
	  INTEGER	SUMNw,SUMNU1,SUMNU2,SUMNU3,SUMNp,SUMNf,SUMNf1,SUMNf2,SUMNf3,SUMNf4,SUMNRv
      REAL		wild,urban1,urban2,urban3,paddy,forest,forest1,forest2,forest3,forest4,river
      REAL		Nw,Nu1,Nu2,Nu3,Np,Nf,Nf1,Nf2,Nf3,Nf4,NRv
      REAL		GANMA,GANMAS,TOUSUIMS,ASOU
      REAL		GANMAC,BETAC
      REAL		QI1,QI2,QI3,QIdist
      INTEGER	MAXRAIN,BL,RADAR,NZONE,DAM
      REAL		RINTERVAL,MAXCRAIN,MAXCRETEN
      REAL		RSA,F1
	  REAL		par1,par2,par3,par4,par5,par6,par7,par8,par9
      CHARACTER NIEQIYANG*2

!	  COL             : 网格数据的列
!	  ROW             : 网格数据的行
!	  TOTAL           : 总网格数 (=COL*ROW)
!	  TOTAL2          : 流域中的网格数
!	  TOTAL2_AREA     : 流域面积 (km^2)
!	  DT              : 计算的时间步长 (小时)
!	  MAXT            : 时间步数 (总计算时间 = DT*MAXT)
!	  RAINT           : 以DT为步长的降雨时长
!	  GRID            : 正方形网格边长(m)
!	  KB              : 运动波方程在一个网格中的时长？
!	  NODE            : N一个网格中的节点
!	  NCLASS          : 班号
!	  IRIVER_THRESHOLD: 根据流量累积图（acc）区分通道的值
!     BL			  : 0表示没有真实的河流宽度数据可用，1表示有
!     RADAR			  : 0为单降雨，1为降雨分布
!	  DAM			  : 1考虑坝，0不考虑坝
      PARAMETER (DT    = 1.0)
      PARAMETER (KB=10,   NODE=10)
      PARAMETER (BL    = 0)
	  PARAMETER (RADAR = 1)
	  PARAMETER (NZONE = 4)			
	  PARAMETER (MAXT  = 124)
	  PARAMETER (RAINT = 124)
      PARAMETER (IRIVER_THRESHOLD = 10)
      PARAMETER (DAM   = 0)
      PARAMETER (NCLASS= 11)

!	  Roughness Manning Coefficients（曼宁粗糙系数）:
!      PARAMETER (Nw	   = 0.2)		!Wilds（荒地）
!      PARAMETER (Nu1   = 0.01)		!Urban1（城市1）
	  PARAMETER (Nu2   = 0.01)		!Urban2（城市2）未分类
      PARAMETER (Nu3   = 0.013935)	!Urban3（城市3）未分类
!	  PARAMETER (Np    = 0.01)		!农地
!     PARAMETER (Nf    = 0.04403)	!Forest（森林）
      PARAMETER (Nf1    = 0.04403)	!Forest1（森林1）未分类
      PARAMETER (Nf2    = 0.04403)	!Forest2（森林2）未分类
      PARAMETER (Nf3    = 0.04403)	!Forest3（森林3）未分类
      PARAMETER (Nf4    = 0.14403)	!Forest4（森林4）未分类
!     PARAMETER (NRv   = 0.05)	    !River（河）

!↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓Soil Property 土壤属性↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
      
!	  这个部分是基于论文"Development of stage-discharge relationship equation incorporationg saturated-unsaturated flow mechanism"
!	  by Yasuto TACHIKAWA et al, Hydrologic Journal of Japan, No48, Feb. 2004.的理论

!	  D   : 土壤总深度 (mm)      
!	  DH  : 地下径流深度 (=DM+DS, mm)
!	  DM  : DH中的非饱和流部分 (mm)
!	  DS  : DH中的饱和流部分 (mm)
!	  KA  : DS的渗透系数  (mm/s, KAMS为 m/s)
!	  KM  : DM的渗透系数 (mm/s)
!	  BETA: 渗透率 (=KA/KM, 一般为 2~10)
!     GANMAS:总土壤孔隙率（%）
!     GANMAC:土壤毛管孔隙率（%）
!                    0 < h < ASOU*GANMAC (DM)  : 非饱和流部分
!     ASOU*GANMAC (DM) < h < ASOU*GANMAS (DA)  : 饱和流（中间流）部分
!     ASOU*GANMAS (DA) < h				       : 表面流部分
!	  
!	  BETA:渗透率(=KA/KM), BETAC 的值(通常为2-10)
!	  TOUSUIMS(KA)/BETAC = KM

!     PARAMETER (TOUSUIMS= 0.002)	!导水率= KA = event 2
!     PARAMETER (ASOU    = 500.)	!影响土壤深度D  da=D*GANMAS dm=D*GANMAC, 1200
      PARAMETER (GANMAS  = 0.65)	!总孔隙度(毛管孔隙度 + 非毛管孔隙度)
      PARAMETER (GANMAC  = 0.05)	!毛管孔隙度
!     PARAMETER (BETAC   = 4)		!event 2

!↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑Soil Property 土壤属性↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑

!	  水流量的初始条件（立方米/秒）
!	  初始条件是假设流域处于均匀状态，所以QI将分布在盆地上，每个单元格的末端都有一个与上层单元格数量成比例的排泄量。
!	  INICON: 决定每个单元格的初始排泄量
!	  0: 按QI均匀状态分配 
!     1: 读取预先模拟的值     

!	  在整个流域的出口处
      PARAMETER (QI1 =200.01)
!	  在Anegawa DAM Catchment的出口处（运行选择不同的流域有不同的QI参数，原模型默认有3的流域信息选择计算）
      PARAMETER (QI2 =0.1)  
!	  在Nyu DAM Catchment的出口处
      PARAMETER (QI3 =0.1) 
!     PARAMETER (F1 =1.0)  
      PARAMETER (MAXRAIN=160)
      PARAMETER (RINTERVAL=3.0)
      PARAMETER (MAXCRAIN=155.75)
      PARAMETER (MAXCRETEN=110.95)
      PARAMETER (RSA=30000.0)

    
	  CHARACTER NAMENS*60,NAMESLP*60,NAMEORD*60,NAMELU*60,NAMEBSN*60,NAMESTORE*60,&
                NAMEBL*60,NAMEW*60,NAMEORAIN*60,NAMEZ*60,NAMER2*60,NAMEFS*60,&
                NAMER*60,NAMEO*60,NAMEAREA*60,NAMESOIL*60,ANS,GEO,HYD,INTCP,STTYPE

      INTEGER   FILER,FILEO,FILEORAIN,FILEWSTORE,FILER2,FILEBSN
      INTEGER   I,K,DUMMY,MAXORDER 
      INTEGER   TMP1,TMP2,TMP3,TMP4,TMP5
      REAL      TGRID,VCOS,FTMP1,FTMP2,FTMP3
      INTEGER   TCOL,TROW,nobasin,notime,yangqinie
      DOUBLE PRECISION DTMP,FN
      Integer, dimension(:), allocatable :: ROWW,COLL,CONV                             !TOTAL
      Integer, dimension(:), allocatable :: ORDER, NSUC, LANDUSE, ZONE, UPPERCELL		 !TOTAL2
      Real,	   dimension(:), allocatable :: ALPH,ALPHA,SLOPE,LENGTH, WIDTH, AREA		 !TOTAL2
      Real,    dimension(:), allocatable :: SLOPE2,THETA,Mcrit,COSTHETA					 !TOTAL2
      Real,    dimension(:), allocatable :: SLOPEMAX,SLOPEMIN,areaQI			         !TOTAL2
      CHARACTER namebasin*2,simtime*12						        

!↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑参数定义完毕↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
      
      WRITE (*,*)
	  WRITE (*,'(10x,A55)') '====================================================='
	  WRITE (*,'(10x,A55)') '=                                                   ='
      WRITE (*,'(10x,A55)') '=                  CDRMV3水文模型                   ='
      WRITE (*,'(10x,A55)') '=                                                   ='
	  WRITE (*,'(10x,A55)') '= = = = = = = = = = = = = = = = = = = = = = = = = = ='
	  WRITE (*,'(10x,A55)') '=                                                   ='
	  WRITE (*,'(10x,A55)') '=              流域名:                              ='  
	  WRITE (*,'(10x,A55)') '=              [01] 灞河流域                        ='
	  WRITE (*,'(10x,A55)') '=              [02] ANEGAWA DAM CATCHMENT           ='
	  WRITE (*,'(10x,A55)') '=              [03] NYU DAM CATCHMENT               ='
	  WRITE (*,'(10x,A55)') '=                                                   ='
	  WRITE (*,'(10x,A55)') '====================================================='


!↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓使用蒙特卡洛校准↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
	
!	  用蒙特卡洛模拟校准最优点...
	 
	  open(89,FILE='./输出数据/MC/蒙特卡罗参数.txt')
      read(89,*)  !nie 敏感性分析
	  read(89,*) par1,par2,par3,par4,par5,par6,par7,par8,par9
		Nu1=par1      !城市
		Nf=par2       !森林
		NRv=par3      !河道
		F1=par4       !降雨除去植被截留与蒸发后的比例系数
		ASOU=par5     !影响土壤深度
		TOUSUIMS=par6 !导水率
		Betac=par7
        np=par8       !农地
        nw=par9       !荒地
	  close (89)
!↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑使用蒙特卡洛校准↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
 nobasin=1
open(618,FILe='土地利用类型年份读取.dat')
READ(618,*) yangqinie
close(618)
WRITE(*,*)'*******landuse代号:',yangqinie,'*******'

      IF(yangqinie==1)  then 
          NIEQIYANG='00'
      elseIF(yangqinie==2)  then 
          NIEQIYANG= '05'
      elseIF(yangqinie==3) then 
          NIEQIYANG= '10' 
      elseIF(yangqinie==4) then 
          NIEQIYANG= '15' 
      elseIF(yangqinie==5) then 
          NIEQIYANG= '80' 
      elseIF(yangqinie==6) then 
          NIEQIYANG= '95' 
      endif
!######################  输入数据  ######################
      
!     流域流向数据文件
	  NAMENS = './地形数据/bsnFlowDirection.dat'
!     流域面积数据文件
      NAMEAREA = './地形数据/流域网格数与面积.dat'
!     次序数据文件
      NAMEORD = './地形数据/order.dat'             !查询位置order号可以查询输出河道号
!     流量累积路数据文件
      NAMEBL = './地形数据/drainage.dat'
!     坡度坡长数据文件
      NAMESLP = './地形数据/Slope.dat'
!     土地利用数据文件
      NAMELU = './地形数据/landuse.dat'

!	  用于设置每个目标点/位置的QI
	  open(45,FILE='./地形数据/QIarea.dat')

!     河道宽度文件 (如果'bl'为0，则不需要输入)
      NAMEW = './地形数据/河宽.dat'

!     保留数据文件（如果'filer2'为0，则不需要输入。）)
      NAMER2= './降雨数据/retention.dat'
!     0表示没有可用的保留数据，1表示有。
      FILER2= 0

!     Radar等于0为单降雨，等于1为空间分布降雨
	  NAMEZ = './降雨数据/降雨分区.dat'
!     降雨量数据输入
      NAMER = './降雨数据/降雨.dat'

!######################  输出数据  ######################
      
!	  截取 
      NAMEORAIN = './输出数据/流域截留.csv'
!	  流域信息
      NAMEBSN   = './输出数据/流域信息'//NIEQIYANG//'.asc'
!	  径流排放
      NAMEO     = './输出数据/流域出流'//NIEQIYANG//'.csv'
!	  流域蓄水
	  NAMESTORE = './输出数据/流域蓄水.csv'
      
!##################  检查数据是否存在  ##################
      
	  WRITE(*,*)
 2002 WRITE(*,'(1x,A)') '是否使用初始文件名称：'
	  ANS='Y'
      IF((ANS=='Y').OR.(ANS=='y')) THEN
			OPEN(10,FILE=NAMENS,ERR=2010,STATUS='OLD')
			WRITE(*,*) '打开文件 <',NAMENS,'> 成功.  '
			OPEN(15,FILE=NAMEAREA,ERR=2010,STATUS='OLD')
			WRITE(*,*) '打开文件 <',NAMEAREA,'> 成功.'
			OPEN(20,FILE=NAMESLP,ERR=2010,STATUS='OLD')
			WRITE(*,*) '打开文件 <',NAMESLP,'> 成功. '
			OPEN(30,FILE=NAMEORD,ERR=2010,STATUS='OLD')
			WRITE(*,*) '打开文件 <',NAMEORD,'>成功. '
			OPEN(40,FILE=NAMELU,ERR=2010,STATUS='OLD')
			WRITE(*,*) '打开文件 <',NAMELU,'>成功.  '
			OPEN(70,FILE=NAMEBL,ERR=2010,STATUS='OLD')
			WRITE(*,*) '打开文件 <',NAMEBL,'> 成功.  '
			
			IF(BL==1) THEN
			OPEN(80,FILE=NAMEW,ERR=2010,STATUS='OLD')
			WRITE(*,*) '打开文件 <',NAMEW,'>成功 '
			END IF

			IF(RADAR==1) THEN
			OPEN(90,FILE=NAMEZ,ERR=2010,STATUS='OLD')
			WRITE(*,*) '打开文件<',NAMEZ,'> 成功. '
			END IF

			IF(FILER2==1) THEN
			OPEN(100,FILE=NAMER2,ERR=2010,STATUS='OLD')
			WRITE(*,*) '打开文件 <',NAMER2,'> 成功.'
			FILER2 = 100
			END IF

			OPEN(50,FILE=NAMER,ERR=2010,STATUS='OLD')
			WRITE(*,*) '打开文件 <',NAMER,'> 成功.  '
			FILER = 50
			
			OPEN(65,FILE=NAMEORAIN,ERR=2010,STATUS='UNKNOWN')
			WRITE(*,*) '打开文件 <',NAMEORAIN,'> 成功.'
			FILEORAIN = 65

			OPEN(63,FILE=NAMEBSN,ERR=2010,STATUS='UNKNOWN')
			WRITE(*,*) '打开文件 <',NAMEBSN,'> 成功.'
			FILEBSN=63

			OPEN(60,FILE=NAMEO,ERR=2010)
			WRITE(*,*) '打开文件 <',NAMEO,'> 成功.  '
			FILEO = 60

			OPEN(71,FILE=NAMESTORE,ERR=2010)
			WRITE(*,*) '打开文件 <',NAMESTORE,'> 成功.  '			
			FILEWSTORE = 71

2012		CONTINUE
			GOTO 2003
      ELSE 
			WRITE(*,*) ' !! 请输入Y !!'
			GOTO 2002
      END IF

! ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓读取检验数据集↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
 2003 CONTINUE
      WRITE(*,*) ' '
      WRITE(*,*) '**************************读取数据中**************************'

!     数据头文件检查
      READ(10,*) GRID,GRID,COL,ROW   !读流向文件抬头的网格尺寸与行列号
	  TOTAL = COL*ROW
      ALLOCATE (CONV(TOTAL))      !配置内存空间

      READ(15,*) TOTAL2, TOTAL2_AREA   !读取流域面积数据文件，流域内网格数与流域面积
   
      ALLOCATE (ORDER(TOTAL2), NSUC(TOTAL2), LANDUSE(TOTAL2), ZONE(TOTAL2), ALPH(TOTAL2), ALPHA(TOTAL2))
      ALLOCATE (SLOPE(TOTAL2), LENGTH(TOTAL2), WIDTH(TOTAL2), AREA(TOTAL2), UPPERCELL(TOTAL2))
      ALLOCATE (SLOPE2(TOTAL2), ROWW(TOTAL2),COLL(TOTAL2),THETA(TOTAL2),Mcrit(TOTAL2))
      ALLOCATE (SLOPEMAX(TOTAL2), COSTHETA(TOTAL2),SLOPEMIN(TOTAL2),areaQI(TOTAL2))
	  
      !读取其他文件抬头的网格尺寸和行列号与流向文件中的抬头是否一致
      READ(20,*) TGRID,TGRID,TCOL,TROW     
      IF(((TGRID<GRID-0.00001).OR.(GRID+0.00001<TGRID)).OR.&
         (TCOL/=COL).OR.&
         (TROW/=ROW)) STOP '!! 山坡数据出错 !!'

      READ(30,*) TGRID,TGRID,TCOL,TROW
      IF(((TGRID<GRID-0.00001).OR.(GRID+0.00001<TGRID)).OR.&
         (TCOL/=COL).OR.&
         (TROW/=ROW)) STOP '!! Order排序数据出错 !!'

      READ(40,*) TGRID,TGRID,TCOL,TROW
      IF(((TGRID<GRID-0.00001).OR.(GRID+0.00001<TGRID)).OR.&
         (TCOL/=COL).OR.&
         (TROW/=ROW)) STOP '!! 土地利用数据出错 !!'

      READ(70,*) TGRID,TGRID,TCOL,TROW
      IF(((TGRID<GRID-0.00001).OR.(GRID+0.00001<TGRID)).OR.&
         (TCOL/=COL).OR.&
         (TROW/=ROW)) STOP '!! 流量累积数据出错 !!'

      IF(BL==1) THEN
        READ(80,*) TGRID,TGRID,TCOL,TROW
        IF(((TGRID<GRID-0.00001).OR.(GRID+0.00001<TGRID)).OR.&
           (TCOL/=COL).OR.&
           (TROW/=ROW)) STOP '!! 河宽数据出错 !!'
      END IF

	  IF(RADAR==1) THEN
		READ(90,*) TGRID,TGRID,TCOL,TROW
		IF(((TGRID<GRID-0.00001).OR.(GRID+0.00001<TGRID)).OR.&
		   (TCOL/=COL).OR.&
		   (TROW/=ROW)) STOP '!! 降雨分区数据出错 !!'
	  ELSE
		TMP5 = 1
	  END IF
! ↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑读取检验数据集↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
2021  CONTINUE
!     数据编号转换为“CONV（单元号）=流域水文单元号”
      MAXORDER = 0   !初始化参数
      K = 0          !初始化参数
      DO I=1,TOTAL
        READ(30,*) TMP1   !读取ORDER文件
        IF(TMP1>=0) THEN  
          K = K+1
          CONV(I) = K     !K累积ORDER非负值个数，最大值即为流域内网格数，也是流域水文编号，CONV（I）有值的为流域水文编号
          IF(TMP1==0) THEN
            OUTLET = I     !ORDER值为0的是流域出口
!			WRITE (*,'(A,I6)') 'Order= ',I
!			GOTO 2017 
          END IF
        else                 !nie
            CONV(I) = -9999  !nie  非流域内点conv均为-9999，避免为定义值的部分被使用到时发生未知错误
        END IF
        IF(MAXORDER<TMP1) MAXORDER = TMP1   !MAXORDER是OEDER数据中的最大值
      ENDDO
      CLOSE(30)

!	  流域面积 (Km2)
	  Basinarea=K*GRID*GRID/1000000.

      WRITE(*,'(1X,A,I8)') '出口位置		      :',OUTLET
      WRITE(*,'(1X,A,I8)') '最大Order值           :',MAXORDER
      WRITE(*,'(1X,A,I8)') '流域总单元数          :',K

      WRITE(FILEBSN,'(1X,A35,I8)')   '流域代号                            :',nobasin
      WRITE(FILEBSN,*) 
      WRITE(FILEBSN,'(1X,A35,I8)')   '出口位置(行列号一维化之后的编号)    :',OUTLET
      WRITE(FILEBSN,'(1X,A35,I8)')   '最大Order值                         :',MAXORDER
      WRITE(FILEBSN,'(1X,A35,I8)')   '流域总单元数                        :',K
      WRITE(FILEBSN,'(1X,A35,F8.3)') '流域面积 [km2]                      :',Basinarea
      WRITE(FILEBSN,*) 

      IF(K/=TOTAL2) STOP '!! 水文单元总数不正确。 !!'    !K与从流域面积文件里读取的流域内网格数TOTAL2不等则报错
      
!     Outlet No Conversion出口不转换   CONV即为converting缩写，意为转换
      OUTLET2 = CONV(OUTLET)  !流域出口的位置转化为流域水文编号
      K = 0
      !打开ORDER文件和QIarea文件读取抬头为正式读取内容做准备，其他文件在检验抬头那一步就已经将读取位跳到正文所以不用这步
      OPEN(30,FILE=NAMEORD)
      READ(30,*) TGRID,TGRID,TCOL,TROW
      READ(45,*) TGRID,TGRID,TCOL,TROW

      DO I=1,TOTAL      !遍历所有网格

        READ(30,*) TMP1                      !ORDER数据     TMP1 ：ORDER
        READ(20,*) FTMP1,FTMP2               !斜坡数据      FTMP1：坡度sin值， FTMP2：坡长
        READ(10,*) DUMMY,DUMMY,TMP2,DUMMY    !流向数据      TMP2 ：流向
        READ(40,*) TMP3                      !土地利用数据  TMP3 ：土地利用类型编号
        READ(70,*) DUMMY,DUMMY,TMP4,DUMMY    !流量累积数据  TMP4 ：流量累积数
		READ(45,*) QIdist                    !QIarea数据   
        IF(BL==1)		READ(80,*) FTMP3     !河宽数据      FTMP3：河网宽度
		IF(RADAR==1)	READ(90,*) TMP5      !降雨分区数据  TMP5 ：降雨分区号
        IF(TMP1>=0)	THEN
			K=K+1
			ORDER(K) = TMP1     !ORDER有数值的重新排序
			SLOPE(K) = FTMP1    !坡度有数值的也重新排序          
			SLOPE2(K)= asin(FTMP1)*180./3.14     !坡度由sin值转换成角度
			UPPERCELL(K) = TMP4      !流量累积数存在的进行排序 
			areaQI(K)=QIdist         !QI也进行同样排序
			DTMP = DBLE(FTMP1)*DBLE(FTMP1)   !DBLE把数据转换为双精度  sin平方
			DTMP = 1.0 - DTMP               !cos平方

			IF (DTMP < 0.0) stop 'DTMP < 0.0'

			VCOS = REAL(SQRT(DTMP))    !VCOS为坡度的cos值
			LENGTH(K) = FTMP2          !坡长排序
			NSUC(K) = CONV(TMP2)

			IF(BL==0) THEN     !无河宽数据的时候
					AREA(K) = 1.0
					WIDTH(K) = GRID*GRID / (VCOS*LENGTH(K))
!				IF(TMP4==0) THEN
				IF(TMP4 <= IRIVER_THRESHOLD) THEN    !当流量累计数小于等于河流阈值，前面设置的阈值常数为10   
					LANDUSE(K) = TMP3                !则该单元的土地利用类型为原来的土地利用类型号
					IF(LANDUSE(K)<0) LANDUSE(K)=6    !如果土地利用类型号小于0则赋为6号
				ELSE                                  
					LANDUSE(K) = NCLASS              !当流量累计数大于等于河流阈值，该单元土地利用类型设置为NCLASS=11
				ENDIF
			ELSE                !有河宽数据
!				IF(TMP4==0) THEN
				IF(FTMP3<=0) THEN    
					LANDUSE(K) = TMP3      !河宽小于等于0，不是被河部分，土地利用类型就是原地利用类型
					IF(LANDUSE(K)<0) LANDUSE(K)=6     !如果土地利用类型号小于0则赋为6号
					WIDTH(K) = GRID*GRID / (VCOS*LENGTH(K))
					AREA(K) = 1.0
				ELSE
					LANDUSE(K) = NCLASS
!					WIDTH(K) = FTMP3
					WIDTH(K) =GRID*GRID / (VCOS*LENGTH(K))
					AREA(K) = VCOS*LENGTH(K)*WIDTH(K)/(GRID*GRID)
				END IF
			ENDIF
					ZONE(K) = TMP5
      ENDIF
      ENDDO

      CLOSE(10)	
      CLOSE(20)
      CLOSE(30)
      CLOSE(40)
	  CLOSE(45)
      CLOSE(70)
      IF(BL==1) CLOSE(80)	  
      IF(RADAR==1) CLOSE(90)

      WRITE(*,*) '****************************初始化****************************'

	  SLP=0.0
	  LengthSLP=0.0
	  SUMNw=0
	  SUMNU1=0
	  SUMNU2=0
	  SUMNU3=0
	  SUMNp=0
	  SUMNf=0
	  SUMNf1=0
	  SUMNf2=0
	  SUMNf3=0
	  SUMNf4=0
	  SUMNRv=0
!****************************判断网格Land Use赋值FN***********************************
      DO I=1,TOTAL2

			IF(LANDUSE(I)==1) FN = Nw         !荒地    
			IF(LANDUSE(I)==2) FN = Nu1        !城市
			IF(LANDUSE(I)==3) FN = Nu2        !城市
			IF(LANDUSE(I)==4) FN = Nu3        !城市
			IF(LANDUSE(I)==5) FN = Np         !农地
			IF(LANDUSE(I)==6) FN = Nf         !森林
			IF(LANDUSE(I)==7) FN = Nf1        !森林
			IF(LANDUSE(I)==8) FN = Nf2        !森林
			IF(LANDUSE(I)==9) FN = Nf3        !森林
			IF(LANDUSE(I)==10) FN = Nf4       !森林
			IF(LANDUSE(I)==11) FN = Nrv       !河

			THETA(I)=asin(SLOPE(I))*180/3.14
			Costheta(I)=Cos(THETA(I)*3.14/180.)
			Theta(I)=Tan(THETA(I)*3.14/180.)

!			(m^(-1/3)/s) change to (mm^(-1/3)/hr)
			FN = FN/3600.0/10.0
!			From (m/sec) to (mm/hr)
			TOUSUI  = TOUSUIMS*1000.0*3600.0
			TOUSUIC = TOUSUI/BETAC
			IF(ASOU==0)THEN
				TOUSUI  = 0
				TOUSUIC = 0
			END IF
			IF(GANMAC==0)THEN
				TOUSUIC = 0
			ENDIF

			IF (SLOPE(I) < 0.0) then 
                stop 'slope < 0.0'!nie
!            SLOPE(I)=0.01 !nie
            endif!nie
			DTMP = SQRT(SLOPE(I)) 
			ALPH(I) = REAL(DTMP/FN/DBLE(LENGTH(I))/1000.0)

!			Cumulative Slope, Slope Length, and Land Use Area for Each Type
			SLP=SLP+SLOPE2(I)
			LengthSLP=LengthSLP+LENGTH(I)

			IF (LANDUSE(I)==1)     THEN
				SUMNw=SUMNw+1
			ELSEIF (LANDUSE(I)==2) THEN
				SUMNu1=SUMNu1+1
			ELSEIF (LANDUSE(I)==3) THEN
				SUMNu2=SUMNu2+1
			ELSEIF (LANDUSE(I)==4) THEN
				SUMNu3=SUMNu3+1
			ELSEIF (LANDUSE(I)==5) THEN
				SUMNp=SUMNp+1
			ELSEIF (LANDUSE(I)==6) THEN
				SUMNf=SUMNf+1
			ELSEIF (LANDUSE(I)==6) THEN
				SUMNf1=SUMNf1+1
			ELSEIF (LANDUSE(I)==6) THEN
				SUMNf2=SUMNf2+1
			ELSEIF (LANDUSE(I)==6) THEN
				SUMNf3=SUMNf3+1
			ELSEIF (LANDUSE(I)==6) THEN
				SUMNf4=SUMNf4+1
			ELSE
				SUMNRv=SUMNRv+1
			ENDIF
			SLOPEMAX(I)=SLOPE2(I)
			SLOPEMIN(I)=SLOPE2(I)
      ENDDO

!			Maximum and Minimum Slope
			DO I=2,TOTAL2
				IF(SLOPEMAX(I)>=SLOPEMAX(I-1)) then
					SLOPEMAX(I)=SLOPEMAX(I)
				ELSE
					SLOPEMAX(I)=SLOPEMAX(I-1)
				ENDIF
					Max=SLOPEMAX(I)
			ENDDO
			DO I=2,TOTAL2
				IF(SLOPEMIN(I)<=SLOPEMIN(I-1)) then
					SLOPEMIN(I)=SLOPEMIN(I)
				ELSE
					SLOPEMIN(I)=SLOPEMIN(I-1)
				ENDIF
					Min=SLOPEMIN(I) 
			ENDDO

!			平均坡度
			SLPav=SLP/TOTAL2	
			LengthSLPav=LengthSLP/TOTAL2
!			Area for Each Land Use Unit:
			wild=SUMNw*GRID*GRID/1000000.
			urban1=SUMNu1*GRID*GRID/1000000.
			urban2=SUMNu2*GRID*GRID/1000000.
			urban3=SUMNu3*GRID*GRID/1000000.
			paddy=SUMNp*GRID*GRID/1000000.
			forest=SUMNf*GRID*GRID/1000000.
			forest1=SUMNf1*GRID*GRID/1000000.
			forest2=SUMNf2*GRID*GRID/1000000.
			forest3=SUMNf3*GRID*GRID/1000000.
			forest4=SUMNf4*GRID*GRID/1000000.
			river=SUMNrv*GRID*GRID/1000000.

			WRITE(FILEBSN,'(1X,A30,F8.4)') '最大坡度      (度):',Max
			WRITE(FILEBSN,'(1X,A30,F8.5)') '最小坡度      (度):',Min
			WRITE(FILEBSN,'(1X,A30,F8.4)') '平均坡度      (度):',SLPav
			WRITE(FILEBSN,'(1X,A30,F8.2)') '平均坡长      (米):',LengthSLPav
			WRITE(FILEBSN,*)
			WRITE(FILEBSN,'(1X,A30)') '各土地利用面积 (km2)    '
			WRITE(FILEBSN,'(1X,A30,F8.2)') '荒地 (1)               :',wild
			WRITE(FILEBSN,'(1X,A30,F8.2)') '城市 (2)               :',urban1
			WRITE(FILEBSN,'(1X,A30,F8.2)') '未分类Urban2 (3)       :',urban2
			WRITE(FILEBSN,'(1X,A30,F8.2)') '未分类Urban3 (4)       :',urban3
			WRITE(FILEBSN,'(1X,A30,F8.2)') '农地 (5)               :',paddy
			WRITE(FILEBSN,'(1X,A30,F8.2)') '林地 (6)               :',forest
			WRITE(FILEBSN,'(1X,A30,F8.2)') '未分类Forest (7)       :',forest1
			WRITE(FILEBSN,'(1X,A30,F8.2)') '未分类Forest (8)       :',forest2
			WRITE(FILEBSN,'(1X,A30,F8.2)') '未分类Forest (9)       :',forest3
			WRITE(FILEBSN,'(1X,A30,F8.2)') '未分类Forest (10)      :',forest4
			WRITE(FILEBSN,'(1X,A30,F8.2)') '河流 (7)               :',river
			WRITE(FILEBSN,*)
!            WRITE(*,*) Betac !NIE

			WRITE(*,*) '***********************降雨径流计算开始***********************'

			GANMA = GANMAS-GANMAC

			CALL PROCESS3(LANDUSE,NSUC,ZONE,ORDER,&
					OUTLET2,NZONE,TOTAL2,MAXORDER,&
                    MAXT,RAINT,KB,NODE,NCLASS,&
                    ALPH,WIDTH,AREA,nobasin,dam,Theta,Costheta,&
                    GRID,DT,MAXRAIN,RINTERVAL,MAXCRAIN,MAXCRETEN,&
                    FILER,FILER2,FILEORAIN,FILEO,namebasin,INTCP,&
                    FILEWSTORE,RSA,F1,UPPERCELL,QI1,QI2,QI3,areaQI,&
                    GANMA,TOUSUI,ASOU,ALPHA,&
                    GANMAC,TOUSUIC,SLOPE,SLOPE2,LENGTH,&
					TCOL,TROW,CONV)
	  WRITE(*,*)
      WRITE(*,*) '************************计算完毕么么哒************************'
 
	  IF(FILER2/=0) CLOSE(100)
      CLOSE(50)
      CLOSE(60)
      STOP
 2010 CONTINUE
      WRITE(*,*) '**********************初始化文件名称错误**********************'
      READ(*,*)
      STOP
      END PROGRAM CDRMV3
