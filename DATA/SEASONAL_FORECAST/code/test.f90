program postnarrmonth
use netcdf
implicit none
integer,parameter::clon=296,clat=292, dims=3
integer :: status, ncid1
integer :: i,ii,j,jj,k,kk,kkk,varid1
integer :: start(dims),count(dims),dimids(dims)
real mask(clon,clat)
character*255 ext
integer basin_dimid, fcst_dimid, year_dimid, bmon_dimid

ext='/home/freeze/water_monitor/African_Drought_Monitor/DATA/3B42RT_BC/MONTHLY/'
call check( nf90_open(trim(ext)//'3B42RT_200606_monthly_0.250deg.nc', NF90_NOWRITE, ncid1) )
call check( nf90_inq_varid(ncid1, 'prec', varid1) )
start=(/1,1,1/);  count=(/clon,clat,1/)
call check( nf90_get_var(ncid1, varid1, mask, start, count) )
call check( nf90_close(ncid1) )
print*,'mask',mask(100,:)

contains
  subroutine check(status)
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then
     print *, trim(nf90_strerror(status))
     stop "Stopped"
  end if
  end subroutine check
end
