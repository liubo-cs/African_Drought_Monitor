!This functions determines if itis a leap year
! returns 1, if yes and 0 if no    
integer function isleap( iyr )
implicit none
integer :: iyr
if( (mod(iyr,4) .eq. 0 .and. mod(iyr,100) .ne.0) &
    .or. mod(iyr,400) .eq. 0) then
    isleap = 1
else
    isleap = 0
endif
end function