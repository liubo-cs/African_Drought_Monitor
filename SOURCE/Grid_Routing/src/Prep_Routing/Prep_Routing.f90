! This code creates the need input to run the grid routing
! including the mapping to large scale, and creation of parameter files
! Josh Roundy (jroundy@princeton.edu)
! originally writen 6 July 2012
! updated 21 Aug 2013 to improve and fix some bugs
program Prep_Routing
implicit none

!Library functions
real::area_m,dist

!Local variables
integer :: nx,ny,mnx,mny
integer :: status,i,j,hunk,nthr,inx,iny,soid,sid,cnt
integer :: mgid(3)  !unique grid id & x&y index
integer,dimension(:,:,:),allocatable :: direc   !flow direction file, (nx,ny,2) gives x and y of flow direction    
integer,dimension(:,:,:),allocatable :: adir  ! Adjusted direction file to account for flows to ocean and dead end grids used in routing
integer,dimension(:,:,:),allocatable :: msid  !Id of model grid cell that contributes to the grid

real :: mxmin,mymin,mxres,myres,xmin,ymin,xres,yres 
real :: lat,lon,d2s,vmin,vstr,vnstr,mind
real :: mll(2) 
real,dimension(4):: lldeg
real,dimension(:,:),allocatable :: so !Stream order of grid cell
real,dimension(:,:),allocatable :: are !Area of grid cell 
real,dimension(:,:),allocatable :: v !velocity of water in grid cell
real,dimension(:,:),allocatable :: l !distance to next grid cell
real,dimension(:,:),allocatable :: slp !slope of grid cell
real,dimension(:,:),allocatable :: mask !mask of model output 

character(len=120) :: basin         !filename, this file provides the basin list
character(len=120) :: outdir,outfile,inputfile,maskfile

namelist /Prep_Routing_nml/ basin,outdir,maskfile,nx,ny,xmin,ymin,xres,yres,&
	mxmin,mymin,myres,mxres,mnx,mny,vnstr,vstr,soid,sid

write(*,*)'##################################################################'
write(*,*)'Start Prep Routing .....'
write(*,*)
write(*,*)

!default Parameters
vnstr = 0.5 !non stream velocity
vstr = 2 !stream velocity
soid = 1 !cut off for what is a stream (stream order)
sid = 1 !1-use slope in calculating velocity for non-streams
mind = 0.6; !minium travel distance for grid cell

!----------------Read Namelist------------------
open(33,file='Prep_Routing.nml', status='old')
read(33,nml=Prep_Routing_nml, end=10)
10  close(33)

if(xmin<=0)then
   xmin=xmin+360
endif
if(mxmin<=0)then
   mxmin=mxmin+360
endif

d2s = 24*3600 !convert to seconds

! Parallel parameters
hunk = 1 !chunk size
nthr = 8 !number of threads
call OMP_SET_NUM_THREADS(nthr) ! SET Number of threads
!----------------Read Namelist------------------

!------Read in Data ------
allocate(mask(mnx,mny),stat=status)
open(11,file=trim(maskfile), form='unformatted',status='old',access='direct',recl=mnx*mny)
read(11,rec=1)mask
close(11)

allocate(so(nx,ny))
write(inputfile,'(a,a)')trim(outdir),'_SO.bin'
open(88,file=trim(inputfile),status='old',access='direct',form='unformatted',recl=nx*ny)
read(88,rec=1)so(:,:)
close(88)

allocate(slp(nx,ny))
write(inputfile,'(a,a)')trim(outdir),'_SLP.bin'
open(88,file=trim(inputfile),status='old',access='direct',form='unformatted',recl=nx*ny)
read(88,rec=1)slp(:,:)
close(88)

allocate(direc(nx,ny,2))
write(inputfile,'(a,a)')trim(outdir),'_DIR.bin'
open(88,file=trim(inputfile),status='old',access='direct',form='unformatted',recl=nx*ny)
read(88,rec=1)direc(:,:,1)
read(88,rec=2)direc(:,:,2)
close(88)
!-----End of Read in Data ------

!-----Allocate and initialize variables ------
allocate(adir(nx,ny,2))
allocate(msid(nx,ny,2),stat=status)
allocate(v(nx,ny))
allocate(l(nx,ny))
allocate(are(nx,ny))
adir = direc
msid = -999
v = -999.9
l = -999.9
are = -999.9
!-----End Allocate and initialize variables ------

!----------------Map to Course scale and Create Routing Parameters ------------------
cnt = 0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,mgid,mll,lat,lon,inx,iny,lldeg,vmin) REDUCTION(+:cnt)
!$OMP DO SCHEDULE(DYNAMIC,HUNK)
do i = 1,ny
  do j = 1,nx
    if(slp(j,i) > -1)then
      lon = xmin + (float(j)-1)*xres
      if(lon < 0)then
	lon = lon + 360
      endif
      lat = ymin + (float(i)-1)*yres
      call point2grid(lon,lat,mxmin,mymin,mxres,myres,mnx,mny,mgid,mll)
      if(mgid(2) > 0 .and. mgid(2) <= mnx .and. mgid(3) > 0 .and. mgid(3) <= mny)then
	if(mask(mgid(2),mgid(3)) > 0)then
	  are(j,i) = area_m(lldeg(1),yres)
	  msid(j,i,1) = mgid(2)
	  msid(j,i,2) = mgid(3)
	endif
      endif
      if(direc(j,i,1) > 0 .and. direc(j,i,2) > 0)then
	inx = direc(j,i,1)
	iny = direc(j,i,2)
	if(inx > 0 .and. inx <= nx .and. iny > 0 .and. iny <= ny)then
	  if(direc(inx,iny,1) < 1)then
	    adir(inx,iny,1) = 0
	    adir(inx,iny,2) = 0
	  endif
	  lldeg(2) = lon
	  lldeg(1) = lat
	  lldeg(4) = xmin + (float(inx)-1)*xres
	  if(lldeg(4) < 0)then
	    lldeg(4) = lldeg(4) + 360
	  endif
	  lldeg(3) = ymin + (float(iny)-1)*yres
	  l(j,i) = dist(lldeg)
	  if(so(j,i) > soid)then
	    v(j,i) = vstr
	  else
	    if(sid > 0)then
	      v(j,i) = vnstr*sqrt(slp(j,i))
	    else
	      v(j,i) = vnstr
	    endif
	  endif
	  vmin = mind*l(j,i)/d2s
	  if(v(j,i) < vmin)then
	    v(j,i) = vmin
	    cnt = cnt + 1
	  endif
	else
	  !-Grids that flow outside the domain
	  adir(j,i,1) = 0
	  adir(j,i,2) = 0
	endif
      else
	adir(j,i,1) = 0
	adir(j,i,2) = 0
      endif
    endif
  enddo
enddo
!$OMP END DO
!$OMP END PARALLEL
write(*,'(i8,a,f4.1,a)')cnt,' velocities reset to minimum travel distance of ',mind,' of length'
!----------------End Map to Course scale and Create Routing Parameters ------------------

!----------------Write  Parameters ------------------
write(outfile,'(a,a)')trim(outdir),'_Model_Grid.bin'
open(99,file=trim(outfile), access='direct',status='unknown', form='unformatted',recl=nx*ny)
write(99,rec=1)msid(:,:,1)
write(99,rec=2)msid(:,:,2)
close(99)
deallocate(msid)

write(outfile,'(a,a)')trim(outdir),'_ADIR.bin'
open(99,file=trim(outfile), access='direct',status='unknown', form='unformatted',recl=nx*ny)
write(99,rec=1)adir(:,:,1)
write(99,rec=2)adir(:,:,2)
close(99)
deallocate(adir)

write(outfile,'(a,a)')trim(outdir),'_VEL.bin'
open(99,file=trim(outfile), access='direct',status='unknown', form='unformatted',recl=nx*ny)
write(99,rec=1)v
close(99)
deallocate(v)

write(outfile,'(a,a)')trim(outdir),'_ARE.bin'
open(99,file=trim(outfile), access='direct',status='unknown', form='unformatted',recl=nx*ny)
write(99,rec=1)are
close(99)
deallocate(are)

write(outfile,'(a,a)')trim(outdir),'_LEN.bin'
open(99,file=trim(outfile), access='direct',status='unknown', form='unformatted',recl=nx*ny)
write(99,rec=1)l
close(99)
deallocate(l)
!----------------End Write  Parameters ------------------

write(*,*)
write(*,*)
write(*,*)'End of Prep_Routing'
write(*,*)'##################################################################'

    
endprogram    