!This functions returns the 2 width month string
character(len=2) function zeroMonth(mn)
implicit none
integer :: mn
character(len=2),dimension(12) :: cmonth = (/'01','02','03','04','05','06','07','08','09','10','11','12'/)
zeroMonth = cmonth(mn)
end function