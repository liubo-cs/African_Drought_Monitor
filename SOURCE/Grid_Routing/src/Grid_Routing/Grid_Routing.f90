! Routs the runoff and baseflow from the model to calculate the flow
! at every grid cell
! Josh Roundy (jroundy@princeton.edu)
! originally written 6 July 2012
! Updated with bug fixes and improvements 21 Aug 2013
program Grid_Routing
implicit none
!Library functions
integer::isleap,daysMonth
character(len=2) :: zeroMonth

!dimension variables
integer :: nx,ny,mnx,mny
integer :: i,j,hunk,nth,mxg,it,k,ig,irec,dx,dy
integer :: syr,smn,sdy,m,y,d,nvr,rid,bid,ist,ost,tid,nt,mtd,ivar,oid,mdy,nd
integer,dimension(:,:,:),allocatable :: direc  ! Adjusted direction file to account for flows to ocean and dead end grids used in routing
integer,dimension(:,:,:),allocatable :: mid !id of model grid cell that contributes to the grid 
integer,dimension(:,:),allocatable :: gid !Index of downstream grids
     
real :: xres,yres,xmin,ymin,mxmin,mymin,mxres,myres
real :: trnf,m3s,osum,rsum,tsum,bpr,ttsum,dsum,tmsum,flow,fmsum,dmsum,oflow
real,dimension(:),allocatable :: pg !percent contribution of flow to each grid
real,dimension(:,:),allocatable :: stm,stm0,stmR !Stream State variables
real,dimension(:,:),allocatable :: flw,flwR !Flow variables
real,dimension(:,:),allocatable :: vel,are,l2g !Routing inputs velociy, area and length to next grid
real,dimension(:,:),allocatable :: runo,base !variables use to read in model inputs

character(len=2) :: strd
character(len=120) :: basin
character(len=120) :: outfile,inputfile,indir,modDir,statedir,streamdir

namelist /Grid_Routing_nml/ basin,modDir,streamdir,indir,statedir,nx,ny,xmin,ymin,xres,yres,&
	mxmin,mymin,myres,mxres,mnx,mny,syr,smn,sdy,nt,tid,ist,ost,mtd,nvr,rid,bid,nth,mxg

!--------------Default Parameters----------------------
syr = 1982
smn = 1 !start month of simulation
sdy = 1 !Start Day of simulation
nt = 2 !number of time steps to run (see tid)
tid = 1 !units of nt, 0-days, 1-months
ist = 1 !1-read in previous state, -1 start with zero
ost = 2 !output of state file 0-none,1-everyday,2-last day of month,3-end of time steps
mtd = 1 !Model output time step, 1-monthly, 2-daily
nvr = 21 !number of variables in the model output file
rid = 3 ! variable number for runoff
bid = 4 ! variable number for baseflow
nth = 8 !number of threads to run
mxg = 200 !maxium number of cells that the flow can be routed through, depends on resolution of topography
!constants
m3s = 1000.0*3600.0*24.0 !convert to m3/s
!--------------End Default Parameters----------------------

!----------------Read Namelist------------------
open(33,file='Grid_Routing.nml', status='old')
read(33,nml=Grid_Routing_nml, end=10)
! write(*,nml=Grid_Routing_nml)
10  close(33)
if(xmin<=0)then
   xmin=xmin+360
endif

if(mxmin<=0)then
   mxmin=mxmin+360
endif
!----------------Read Namelist------------------

write(*,*)'##################################################################'
write(*,*)'Start Grid Routing ....'
write(*,*)

!--------------OMP Parameters----------------------
hunk = 1 !chunk size
call OMP_SET_NUM_THREADS(nth) ! SET Number of threads
!--------------End OMP Parameters----------------------

!--------------Read in Routing Parameters ----------------------
allocate(direc(nx,ny,2))
write(inputfile,'(a,a)')trim(indir),'_ADIR.bin'
open(88,file=trim(inputfile),status='old',access='direct',form='unformatted',recl=nx*ny)
read(88,rec=1)direc(:,:,1)
read(88,rec=2)direc(:,:,2)
close(88)

allocate(mid(nx,ny,2))
write(inputfile,'(a,a)')trim(indir),'_Model_Grid.bin'
open(88,file=trim(inputfile),status='old',access='direct',form='unformatted',recl=nx*ny)
read(88,rec=1)mid(:,:,1)
read(88,rec=2)mid(:,:,2)
close(88)

allocate(vel(nx,ny))
write(inputfile,'(a,a)')trim(indir),'_VEL.bin'
open(88,file=trim(inputfile),status='old',access='direct',form='unformatted',recl=nx*ny)
read(88,rec=1)vel(:,:) !Accumulation Area
close(88)

allocate(are(nx,ny))
write(inputfile,'(a,a)')trim(indir),'_ARE.bin'
open(88,file=trim(inputfile),status='old',access='direct',form='unformatted',recl=nx*ny)
read(88,rec=1)are(:,:) !Area of grid cell
close(88)

allocate(l2g(nx,ny))
write(inputfile,'(a,a)')trim(indir),'_LEN.bin'
open(88,file=trim(inputfile),status='old',access='direct',form='unformatted',recl=nx*ny)
read(88,rec=1)l2g(:,:) !length to next grid
close(88)
!--------------End Read in Routing Parameters ----------------------

!--------------Allocate and initialize variables ----------------------
allocate(runo(mnx,mny)) !variable to read in runoff
allocate(base(mnx,mny)) !variable to read in baseflow
allocate(stm(nx,ny)) !Stream State variable
allocate(stm0(nx,ny)) !Initial Stream State variable
allocate(stmR(nx,ny)) !Masked out initialization of stream state
allocate(flw(nx,ny))  !flow variable
allocate(flwR(nx,ny)) !Initialized flow variable

stmR = -999.9
flwR = -999.9
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,pg,gid,ig,trnf,k)
!$OMP DO SCHEDULE(DYNAMIC,HUNK)
do i = 1,ny
  do j = 1, nx
    if(are(j,i) > 0.0)then
      stmR(j,i) = 0.0
    endif
    if(direc(j,i,1) > -1)then
      flwR(j,i) = 0.0
    endif
  enddo
enddo
!$OMP END DO
!$OMP END PARALLEL
!--------------End Allocate and initialize variables ----------------------

!------------Initialize State------------
if(ist > 0)then
  y = syr
  m = smn
  d = sdy
  d = d -1
  if(d < 1)then
    m = m -1
    if(m < 1)then
      m = 1
      y = y -1
    endif
    d = daysMonth(m)
  endif
  if(m == 2)then
    d = d + isleap(y)
  endif    
  if(d < 10)then
    write(strd,'(a1,i1)')'0',d
  else
    write(strd,'(i2)')d
  endif
  write(outfile,'(a,a,i4,a2,a2,a)')trim(statedir),'/State_',y,zeroMonth(m),strd,'.bin'
  open(99,file=trim(outfile), access='direct',status='unknown', form='unformatted',recl=nx*ny)
  read(99,rec=1)stm0
  close(99)
  !---------------This section just checks sums and figures out accuracy, not neccesary for the alorithm
  dsum = 0.0
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j) REDUCTION(+:dsum)
  !$OMP DO SCHEDULE(DYNAMIC,HUNK)
  do i = 1,ny
    do j = 1, nx
      if(stm0(j,i) > 0.0)then
	dsum = dsum + stm0(j,i)
      endif
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  !---------------End This section
else
  stm0 = stmR
  dsum = 0.0
endif
!------------Done Initialize State------------

!------------Initialize time variables------------
if(tid == 0)then !nt is in days
  nd = nt
elseif(tid == 1)then !nt is in months
  m = smn
  y = syr
  nd = 0
  do i=1,nt
    mdy = daysMonth(m)
    if(m.eq.2) then
	mdy= mdy + isleap(y)
    endif
    nd = nd + mdy
    m = m + 1
    if(m > 12)then
      m = 1
      y = y + 1
    endif
  enddo
else
  print*,'Incorrect time id',tid
  STOP
endif

y = syr
m = smn
d = sdy
mdy = daysMonth(m)
if(m.eq.2) then
    mdy= mdy + isleap(y)
endif
irec = -1
oid = -1
ttsum = 0.0
fmsum = 0.0
dmsum = dsum
!------------End Initialize time variables------------
write(*,*)
write(*,*)'Routing ',trim(basin),' for ',nd,' days .....'
write(*,*)
write(*,'(a15,4(a10,1x),a6)')'','Land','Ocean','Total','UR-Total','% Diff'
!----------------Rout the Model Flow------------------
do it=1,nd
  if(d < 10)then
    write(strd,'(a1,i1)')'0',d
  else
    write(strd,'(i2)')d
  endif
  !--------------- Start Reading Model Data
  if(mtd == 1)then
    if(oid < 0)then
!       print*,'READING',y,zeromonth(m),'0100'
      write(inputfile,'(a,i4,a2,a4)')trim(modDir),y,zeromonth(m),'0100'
      open(47,file=trim(inputfile),status='old',access='direct',form='unformatted',recl=mnx*mny)
      if(irec < 0)then
	irec=(sdy-1)*nvr+1
      endif
      oid = 1
    endif
  elseif(mtd == 2)then
!     print*,'READING',y,zeromonth(m),strd,'00'
    write(inputfile,'(a,i4,a2,a,a)')trim(modDir),y,zeromonth(m),strd,'00'
    open(47,file=trim(inputfile),status='old',access='direct',form='unformatted',recl=mnx*mny)
    irec=1
  endif
  !Read Data
  do ivar = 1,nvr
    if(ivar == rid)then
      read(47,rec=irec)runo   !unit mm/day
    endif
    if(ivar == bid)then
      read(47,rec=irec)base   !unit mm/day
    endif
    irec=irec+1
  enddo
  if(mtd == 2)then
    close(47)
  endif
  !--------------- End Reading Model Data

  !--------------- Start of Routing Flow
  stm = stmR
  flw = flwR
  tsum = 0.0
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,pg,gid,ig,trnf,bpr,k,flow,oflow) REDUCTION(+:tsum,flw,stm)
  allocate(gid(mxg,2)) !index of downstream grids
  allocate(pg(mxg)) !percent contribution to each grid
  !$OMP DO SCHEDULE(DYNAMIC,HUNK)
  do i = 1,ny
    do j = 1,nx
      if(are(j,i) > 0.0)then
! 	if(runo(mid(j,i,1),mid(j,i,2)) > -0.0001 .and. base(mid(j,i,1),mid(j,i,2)) > -0.0001)then !Dont need this if you trust your model to not give negative runoff or baseflow
	  call Rout_Grid(j,i,nx,ny,direc,vel,l2g,mxg,dx,dy,pg,gid,ig)
	  bpr = (runo(mid(j,i,1),mid(j,i,2))+base(mid(j,i,1),mid(j,i,2)))*are(j,i)/m3s !convert m3/s
	  trnf = bpr + stm0(j,i)
	  tsum = tsum + bpr
	  flow = 1.0
	  oflow = 1.0
	  do k = 1,ig
	    flow = flow - pg(k)
	    if(direc(gid(k,1),gid(k,2),1) > 0)then !Check to make sure it is not an end cell
	      stm(gid(k,1),gid(k,2)) = stm(gid(k,1),gid(k,2)) + trnf*pg(k)
	      flw(gid(k,1),gid(k,2)) = flw(gid(k,1),gid(k,2)) + trnf*flow
	    else
	      flw(gid(k,1),gid(k,2)) = flw(gid(k,1),gid(k,2)) + trnf*oflow !The ocean and dead cells are different because everything that enters it is considered flow
	    endif
	    oflow = oflow - pg(k)
	  enddo
! 	endif
      endif
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  !--------------- End of Routing Flow

  !---------------This section just checks sums and figures out accuracy, not neccesary for the alorithm
  rsum = 0.0
  osum = 0.0
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j) REDUCTION(+:rsum,osum)
  !$OMP DO SCHEDULE(DYNAMIC,HUNK)
  do i = 1,ny
    do j = 1, nx
      if(stm(j,i) > 0.0)then
	rsum = rsum + stm(j,i)
      endif
      if(direc(j,i,1) == 0)then
	osum = osum + flw(j,i)
      endif
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  tmsum = rsum+osum
  fmsum = fmsum + osum
  ttsum = ttsum + tsum
  !--------------- End of Check section

  !--------------- Increment variables and write output
  stm0 = stm
  write(outfile,'(a,a,i4,a2,a2,a)')trim(streamdir),'/Streamflow_',y,zeromonth(m),strd,'.bin'
  open(99,file=trim(outfile), access='direct',status='unknown', form='unformatted',recl=nx*ny)
  write(99,rec=1)flw 
  close(99)

  if(ost == 1)then !output state everyday
    write(outfile,'(a,a,i4,a2,a2,a)')trim(statedir),'/State_',y,zeromonth(m),strd,'.bin'
    open(99,file=trim(outfile), access='direct',status='unknown', form='unformatted',recl=nx*ny)
    write(99,rec=1)stm 
    close(99)
  endif

  if(ost == 3 .and. it==nd)then !output state only at the end of the time run
    write(outfile,'(a,a,i4,a2,a2,a)')trim(statedir),'/State_',y,zeromonth(m),strd,'.bin'
    open(99,file=trim(outfile), access='direct',status='unknown', form='unformatted',recl=nx*ny)
    write(99,rec=1)stm 
    close(99)
  endif
  if(d == mdy)then !End of the month 
    write(*,'(i4,a1,a2,a,4(f10.0,1x),f6.3)')y,'-',zeromonth(m),' Totals:',rsum-dmsum,fmsum,rsum-dmsum+fmsum,ttsum,(fmsum+rsum-dmsum)/ttsum
    if(mtd == 1)then
      close(47) !close monthly filename
      oid = -1
      irec = 1
    endif
    if(ost == 2)then !output state
      write(outfile,'(a,a,i4,a2,a2,a)')trim(statedir),'/State_',y,zeromonth(m),strd,'.bin'
      open(99,file=trim(outfile), access='direct',status='unknown', form='unformatted',recl=nx*ny)
      write(99,rec=1)stm 
      close(99)
    endif
    ttsum = 0.0
    fmsum = 0.0
    dmsum = rsum
    !increment to next month
    d = 1
    m = m + 1
    if(m > 12)then
      m = 1
      y = y + 1
    endif
    mdy = daysMonth(m)
    if(m.eq.2) then
	mdy= mdy + isleap(y)
    endif
  else
    d = d + 1
  endif
  !--------------- END Increment variables and write output
enddo

write(*,*)
write(*,*)
write(*,*)'End Grid_Routing'
write(*,*)'##################################################################'

    
endprogram    
