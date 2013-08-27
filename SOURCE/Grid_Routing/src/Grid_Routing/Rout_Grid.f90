!Calculate the percentage of water that 
!moves to which grids during the times step
!jroundy@princeton.edu 6 July 2012
recursive subroutine Rout_Grid(ii,jj,nx,ny,direc,vel,l2g,mxig,dx,dy,p,gid,ig)
implicit none
!inputs
integer,intent(IN) :: ii,jj,nx,ny,mxig,dx,dy
integer,intent(IN) :: direc(nx,ny,2) !direction matrix of flow
real,intent(IN) :: vel(nx,ny) !velocity matrix
real,intent(IN) :: l2g(nx,ny) !length to next grid matrix

!outputs
real,intent(OUT) :: p(mxig)
integer,intent(OUT) :: ig
integer,intent(OUT) :: gid(mxig,2)

!Library functions
real::dist

!local
integer :: ix,iy,iny,inx,i,eflg,oig
real,dimension(4):: lldeg
real :: s2d,dt,tt,dt0,tt0,xx,yy,xm,d,td,psum,val,zsig,tsig

s2d = 3600*24
ix = ii
iy = jj
ig = 1
gid(1,1) = ix
gid(1,2) = iy
tt = 0
dt = 0
tt0 = 0
dt0 = 0
p = 0
eflg = 0

zsig = 0.0001 !significane level of zero
tsig = 0.01 !significance level of sum check at the end
psum = 0.0
do while (tt < 1)
    ig = ig +1
    if(ig > mxig)then
      print*,'Maxium grid travel is too small',ig,mxig,tt,td
      stop
    endif
    inx = direc(ix,iy,1)
    iny = direc(ix,iy,2)
    d = l2g(ix,iy)
    dt = dt + d/1000
    td = (d/vel(ix,iy))/s2d
    tt = tt + td
!     print*,'loop',ig-1,inx,iny,td
    if(tt > 1)then
        xx = (1-tt)/((tt-tt0)/(dt-dt0))+dt !distance the water travels in one full day
        xm = (dt+dt0)/2 !Distance to middle of the two grids
        yy = 1-tt0
        if(xm < xx)then
	    val = xm/xx*yy
            p(ig-1) = p(ig-1)+val
	    psum = psum + val
	    val = yy-xm/xx*yy
            p(ig) = p(ig) + val
	    psum = psum + val
	    if(abs(psum-1.0) > zsig)then
		p(ig) = p(ig) + 1.0-psum
	    endif
	    eflg = 3
        else
	   ig = ig -1
	   val = yy
           p(ig) = p(ig) + val
	   psum = psum + val
	   if(abs(psum-1.0) > zsig)then
	      p(ig) = p(ig) + 1.0-psum
	   endif
	   eflg = 1
	   goto 85
        endif
    else
	val = td/2.0
        p(ig-1) = p(ig-1) + val
	psum = psum + val
        p(ig) = p(ig) + val
	psum = psum + val
	eflg = 4
    endif
    ix = inx
    iy = iny
    if(ix == 0 .and. iy == 0)then !Reached an ocean or dead cell, so water accumulates here.
      ig = ig -1
      val = 1-psum
      p(ig) = p(ig) + val
      psum = psum + val
      if(abs(psum-1.0) > zsig)then
	p(ig-1) = p(ig-1) + 1.0-psum
      endif
      eflg = 2
      goto 85
    endif
    dt0 = dt
    tt0 = tt
    gid(ig,1) = ix
    gid(ig,2) = iy
enddo
85 continue

oig = ig
psum = 0.0
do i = 1,ig
  psum = psum + p(i)
  if(p(i) < zsig)then
    ig = i-1
    goto 75
  endif
enddo
75 continue

if(abs(psum-1.0) > tsig)then
  print*,'Something is messed up',psum,ii,jj,ig,oig,eflg
!   print*, tt,td,tt0
!   print*, dt, d,dt0
!   print*, xx,xm,yy
  print*, p(1:oig)
  STOP
endif

return
end subroutine