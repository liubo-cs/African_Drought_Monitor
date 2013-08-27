!Calculate distance between grid cells
!using the great circle distance using the Haversine formula
!jroundy@princeton.edu 27 Jan 2011
recursive real function dist(lldeg)
implicit none
!inputs
real,dimension(4),intent(IN):: lldeg
!local
real:: dlat,dlon,dt,d2r,pi
real,dimension(4):: ll
REAL, PARAMETER :: R = 6378137 ![m]

pi = atan(1.0) * 4.0
d2r = pi/180
ll = lldeg*d2r
dlat = ll(3)-ll(1)
dlon = ll(4)-ll(2)
dt = 2*asin(sqrt((sin(dlat/2))**2+cos(ll(1))*cos(ll(3))*(sin(dlon/2))**2))
dist = R * dt;    
end function