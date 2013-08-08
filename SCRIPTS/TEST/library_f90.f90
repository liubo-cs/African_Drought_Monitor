subroutine percentileofscore(scores,n,score,rnd,pct)

implicit none
!Declare variables
real*8 :: scores(n)
real*8 :: score,pct,rnd
integer :: n,i,nstrict,nweak,nrnd
!f2py intent(in) :: scores,score,rnd
!f2py intent(hide) :: n
!f2py intent(out) :: pct
!Find the scores below the given score
nweak = 0
nstrict = 0
do i = 1,n
 if (score .ge. scores(i)) nweak = nweak + 1
 if (score .gt. scores(i)) nstrict = nstrict + 1
enddo
!Compute the random percentile between the limits
nrnd = nstrict + rnd*(nweak - nstrict)

!Compute the percentile
!pct = 100.0*(float(nweak+nstrict)/float(2*n))
pct = 100.0*(float(nrnd)/float(n))

end subroutine percentileofscore

subroutine scoreatpercentile(scores,n,pct,score)

implicit none
!Declare variables
real*8 :: scores(n),quantiles(n),percentiles(n),scores_copy(n)
real*8 :: score,pct
integer :: n,i,pos(1)
!f2py intent(in) :: scores,pct
!f2py intent(hide) :: n
!f2py intent(out) :: score
scores_copy = scores
!Sort the scores
do i = 1,n
 quantiles(i) = minval(scores_copy)
 pos = minloc(scores_copy)
 scores_copy(pos) = 999999.0
 percentiles(i) = 100*(real(i)-1)/(real(n)-1)
 if (i .gt. 1 .and. pct .ge. percentiles(i-1) .and. pct .le. percentiles(i)) then
  score = quantiles(i-1) + (quantiles(i) - quantiles(i-1))*(pct - percentiles(i-1))/(percentiles(i)-percentiles(i-1))
  exit
 endif
enddo

end subroutine scoreatpercentile
