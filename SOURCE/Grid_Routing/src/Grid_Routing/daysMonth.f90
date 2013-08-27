!This functions returns the 3 letter month string  
integer function daysMonth(mn)
implicit none
integer :: mn
integer,dimension(12) :: x = (/31,28,31,30,31,30,31,31,30,31,30,31/)
daysMonth = x(mn)
end function