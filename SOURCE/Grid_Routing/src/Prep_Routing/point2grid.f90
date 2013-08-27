!Maps a point (lon,lat) to a grid given by xmin,ymin
! xres, yres,nx,ny
! returns the grid_id,x,y,lat,lon of the grid
! jroundy@princton.edu Feb 2011
recursive subroutine point2grid(lon,lat,xmin,ymin,xres,yres,nx,ny,gid,ll)
IMPLICIT NONE
integer,intent(IN) :: nx,ny
real,intent(IN) :: lon,lat,xmin,ymin,xres,yres
integer,intent(OUT) :: gid(3)!unique grid id & x&y index
real,intent(OUT) :: ll(2) !lat lon of containing grid
!local
integer :: x,y

x = nint((lon-xmin)/xres)+1
y = nint((lat-ymin)/yres)+1
! if(x < 1 .or. x > nx)then
!     print*,'X point out of grid range',x,nx,lon,xmin
!     gid(1:3) = -1
!     ll(1:2) = -999.9
!     goto 12
! endif
! if(y < 1 .or. y > ny)then
!     print*,'Y point out of grid range',y,ny,lat,ymin
!     gid(1:3) = -1
!     ll(1:2) = -999.9
!     goto 12
! endif
gid(1) = nx*(y-1)+x
gid(2) = x
gid(3) = y
ll(1) = (y-1)*yres+ymin
ll(2) = (x-1)*xres+xmin

12 continue
end subroutine