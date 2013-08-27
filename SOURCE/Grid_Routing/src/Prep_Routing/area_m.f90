!Calculates the area of a cell assuming xres = yres
!returns the area in meters squred
recursive real function area_m(lat,dphi)
IMPLICIT NONE
! Calcualtes area and returns it in [m^2]
!input
real ::   lat,dphi
!local
real ::   pi
REAL, PARAMETER :: R = 6378137 ![m]

pi = atan(1.0) * 4.0
area_m = (2.0*pi*R*dphi/360.0)**2.0 * cos((lat)*pi/180.0)

end function