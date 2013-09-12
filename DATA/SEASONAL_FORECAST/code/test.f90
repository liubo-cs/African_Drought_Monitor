program postnarrmonth
use netcdf
implicit none
integer,parameter::clon=296,clat=292, dims=3
integer :: status, ncid1,nday
integer :: i,ii,j,jj,k,kk,kkk,varid1
integer :: start(dims),count(dims),dimids(dims)
real mask(clon,clat)
character*255 ext
integer basin_dimid, fcst_dimid, year_dimid, bmon_dimid
integer day(12)
data day/31,28,31,30,31,30,31,31,30,31,30,31/

ext='/home/stream1/monitor/African_Drought_Monitor/DATA/3B42RT_BC/MONTHLY/'
print*,trim(ext)//'3B42RT_BC_201308_monthly_0.250deg.nc'
call check( nf90_open(trim(ext)//'3B42RT_BC_201308_monthly_0.250deg.nc', NF90_NOWRITE, ncid1) )
call check( nf90_inq_varid(ncid1, 'prec', varid1) )
start=(/1,1,1/);  count=(/clon,clat,1/)
call check( nf90_get_var(ncid1, varid1, mask, start, count) )
call check( nf90_close(ncid1) )
nday=day(8)
if(mod(2013,4)==0 .and. 8==2) nday=nday+1
mask=mask/nday
print*,mask(100,100)

contains
  subroutine check(status)
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then
     print *, trim(nf90_strerror(status))
     stop "Stopped"
  end if
  end subroutine check
end
