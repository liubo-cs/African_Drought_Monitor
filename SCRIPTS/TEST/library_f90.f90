subroutine percentileofscore(scores,n,score,rnd,pct)

implicit none
!Declare variables
real*8 :: scores(n)
real*8 :: score,pct,rnd,pct_weak,pct_strict
integer :: n,i,nstrict,nweak,nrnd,countn
!f2py intent(in) :: scores,score,rnd
!f2py intent(hide) :: n
!f2py intent(out) :: pct
!Find the scores below the given score
nweak = 0
nstrict = 0
countn = 0
do i = 1,n
 if (score .eq. scores(i)) then
  nweak = nweak + 1
  if (countn .eq. 0) nstrict = nstrict + 1
  countn = countn + 1
 endif
 if (score .gt. scores(i)) then
  nweak = nweak + 1 
  nstrict = nstrict + 1
 endif
enddo
!Compute both percentiles
pct_weak = 100.0*float(nweak-1)/float(n-1)
pct_strict = 100.0*float(nstrict-1)/float(n-1)


!Compute the random percentile between the limits
pct = pct_strict + rnd*(pct_weak - pct_strict)
!nrnd = nstrict + rnd*(nweak - nstrict)

!Compute the percentile
!pct = 100.0*(float(nweak+nstrict)/float(2*n))
!pct = 100.0*(float(nrnd)/float(n))
!print*,nweak,nstrict,rnd,n,nrnd

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
