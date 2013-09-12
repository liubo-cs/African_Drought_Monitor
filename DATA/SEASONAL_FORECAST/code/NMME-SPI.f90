!=======================================================================================================
!
! This package is used to downscale NMME forecast and calculate SPI drought index
! The main program is called "NMMEspi"
!
! Author: Dr. Xing Yuan @ Princeton Univerity, 2013-08-06
! Revised for t2 anomaly and grads input: Xing Yuan, 2013-08-13
! Revised for African forecast at 1/4 degree resolution: Xing Yuan, 2013-08-25
!
! References: Yuan, X., and E. F. Wood, 2013: Multi-model seasonal forecasting of global drought onset.
!
!=======================================================================================================

!===================================================================================
MODULE nrtype
      INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
      INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
      INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
      INTEGER, PARAMETER :: SP = KIND(1.0)
      INTEGER, PARAMETER :: DP = KIND(1.0D0)
      INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
      INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
      INTEGER, PARAMETER :: LGT = KIND(.true.)
      REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
      REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
      REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
      REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
      REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
      REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
      REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
      REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
      TYPE sprs2_sp
            INTEGER(I4B) :: n,len
            REAL(SP), DIMENSION(:), POINTER :: val
            INTEGER(I4B), DIMENSION(:), POINTER :: irow
            INTEGER(I4B), DIMENSION(:), POINTER :: jcol
      END TYPE sprs2_sp
      TYPE sprs2_dp
            INTEGER(I4B) :: n,len
            REAL(DP), DIMENSION(:), POINTER :: val
            INTEGER(I4B), DIMENSION(:), POINTER :: irow
            INTEGER(I4B), DIMENSION(:), POINTER :: jcol
      END TYPE sprs2_dp
END MODULE nrtype
!===================================================================================

!===================================================================================
MODULE nrutil
      USE nrtype
      IMPLICIT NONE
      INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
      INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
      INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
      INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
      INTEGER(I4B), PARAMETER :: NPAR_POLY=8
      INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8
      INTERFACE array_copy
            MODULE PROCEDURE array_copy_r, array_copy_d, array_copy_i
      END INTERFACE
      INTERFACE swap
            MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
                  swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
                  masked_swap_rs,masked_swap_rv,masked_swap_rm
      END INTERFACE
      INTERFACE reallocate
            MODULE PROCEDURE reallocate_rv,reallocate_rm,&
                  reallocate_iv,reallocate_im,reallocate_hv
      END INTERFACE
      INTERFACE imaxloc
            MODULE PROCEDURE imaxloc_r,imaxloc_i
      END INTERFACE
      INTERFACE assert
            MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
      END INTERFACE
      INTERFACE assert_eq
            MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
      END INTERFACE
      INTERFACE arth
            MODULE PROCEDURE arth_r, arth_d, arth_i
      END INTERFACE
      INTERFACE geop
            MODULE PROCEDURE geop_r, geop_d, geop_i, geop_c, geop_dv
      END INTERFACE
      INTERFACE cumsum
            MODULE PROCEDURE cumsum_r,cumsum_i
      END INTERFACE
      INTERFACE poly
            MODULE PROCEDURE poly_rr,poly_rrv,poly_dd,poly_ddv,&
                  poly_rc,poly_cc,poly_msk_rrv,poly_msk_ddv
      END INTERFACE
      INTERFACE poly_term
            MODULE PROCEDURE poly_term_rr,poly_term_cc
      END INTERFACE
      INTERFACE outerprod
            MODULE PROCEDURE outerprod_r,outerprod_d
      END INTERFACE
      INTERFACE outerdiff
            MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
      END INTERFACE
      INTERFACE scatter_add
            MODULE PROCEDURE scatter_add_r,scatter_add_d
      END INTERFACE
      INTERFACE scatter_max
            MODULE PROCEDURE scatter_max_r,scatter_max_d
      END INTERFACE
      INTERFACE diagadd
            MODULE PROCEDURE diagadd_rv,diagadd_r
      END INTERFACE
      INTERFACE diagmult
            MODULE PROCEDURE diagmult_rv,diagmult_r
      END INTERFACE
      INTERFACE get_diag
            MODULE PROCEDURE get_diag_rv, get_diag_dv
      END INTERFACE
      INTERFACE put_diag
            MODULE PROCEDURE put_diag_rv, put_diag_r
      END INTERFACE
CONTAINS
!BL
      SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
      REAL(SP), DIMENSION(:), INTENT(IN) :: src
      REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
      INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
      n_copied=min(size(src),size(dest))
      n_not_copied=size(src)-n_copied
      dest(1:n_copied)=src(1:n_copied)
      END SUBROUTINE array_copy_r
!BL
      SUBROUTINE array_copy_d(src,dest,n_copied,n_not_copied)
      REAL(DP), DIMENSION(:), INTENT(IN) :: src
      REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
      INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
      n_copied=min(size(src),size(dest))
      n_not_copied=size(src)-n_copied
      dest(1:n_copied)=src(1:n_copied)
      END SUBROUTINE array_copy_d
!BL
      SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: src
      INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: dest
      INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
      n_copied=min(size(src),size(dest))
      n_not_copied=size(src)-n_copied
      dest(1:n_copied)=src(1:n_copied)
      END SUBROUTINE array_copy_i
!BL
!BL
      SUBROUTINE swap_i(a,b)
      INTEGER(I4B), INTENT(INOUT) :: a,b
      INTEGER(I4B) :: dum
      dum=a
      a=b
      b=dum
      END SUBROUTINE swap_i
!BL
      SUBROUTINE swap_r(a,b)
      REAL(SP), INTENT(INOUT) :: a,b
      REAL(SP) :: dum
      dum=a
      a=b
      b=dum
      END SUBROUTINE swap_r
!BL
      SUBROUTINE swap_rv(a,b)
      REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
      REAL(SP), DIMENSION(SIZE(a)) :: dum
      dum=a
      a=b
      b=dum
      END SUBROUTINE swap_rv
!BL
      SUBROUTINE swap_c(a,b)
      COMPLEX(SPC), INTENT(INOUT) :: a,b
      COMPLEX(SPC) :: dum
      dum=a
      a=b
      b=dum
      END SUBROUTINE swap_c
!BL
      SUBROUTINE swap_cv(a,b)
      COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
      COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
      dum=a
      a=b
      b=dum
      END SUBROUTINE swap_cv
!BL
      SUBROUTINE swap_cm(a,b)
      COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
      COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
      dum=a
      a=b
      b=dum
      END SUBROUTINE swap_cm
!BL
      SUBROUTINE swap_z(a,b)
      COMPLEX(DPC), INTENT(INOUT) :: a,b
      COMPLEX(DPC) :: dum
      dum=a
      a=b
      b=dum
      END SUBROUTINE swap_z
!BL
      SUBROUTINE swap_zv(a,b)
      COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
      COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
      dum=a
      a=b
      b=dum
      END SUBROUTINE swap_zv
!BL
      SUBROUTINE swap_zm(a,b)
      COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
      COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
      dum=a
      a=b
      b=dum
      END SUBROUTINE swap_zm
!BL
      SUBROUTINE masked_swap_rs(a,b,mask)
      REAL(SP), INTENT(INOUT) :: a,b
      LOGICAL(LGT), INTENT(IN) :: mask
      REAL(SP) :: swp
      if (mask) then
            swp=a
            a=b
            b=swp
      end if
      END SUBROUTINE masked_swap_rs
!BL
      SUBROUTINE masked_swap_rv(a,b,mask)
      REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
      LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
      REAL(SP), DIMENSION(size(a)) :: swp
      where (mask)
            swp=a
            a=b
            b=swp
      end where
      END SUBROUTINE masked_swap_rv
!BL
      SUBROUTINE masked_swap_rm(a,b,mask)
      REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
      LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
      REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
      where (mask)
            swp=a
            a=b
            b=swp
      end where
      END SUBROUTINE masked_swap_rm
!BL
!BL
      FUNCTION reallocate_rv(p,n)
      REAL(SP), DIMENSION(:), POINTER :: p, reallocate_rv
      INTEGER(I4B), INTENT(IN) :: n
      INTEGER(I4B) :: nold,ierr
      allocate(reallocate_rv(n),stat=ierr)
      if (ierr /= 0) call &
            nrerror('reallocate_rv: problem in attempt to allocate memory')
      if (.not. associated(p)) RETURN
      nold=size(p)
      reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
      deallocate(p)
      END FUNCTION reallocate_rv
!BL
      FUNCTION reallocate_iv(p,n)
      INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
      INTEGER(I4B), INTENT(IN) :: n
      INTEGER(I4B) :: nold,ierr
      allocate(reallocate_iv(n),stat=ierr)
      if (ierr /= 0) call &
            nrerror('reallocate_iv: problem in attempt to allocate memory')
      if (.not. associated(p)) RETURN
      nold=size(p)
      reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
      deallocate(p)
      END FUNCTION reallocate_iv
!BL
      FUNCTION reallocate_hv(p,n)
      CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
      INTEGER(I4B), INTENT(IN) :: n
      INTEGER(I4B) :: nold,ierr
      allocate(reallocate_hv(n),stat=ierr)
      if (ierr /= 0) call &
            nrerror('reallocate_hv: problem in attempt to allocate memory')
      if (.not. associated(p)) RETURN
      nold=size(p)
      reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
      deallocate(p)
      END FUNCTION reallocate_hv
!BL
      FUNCTION reallocate_rm(p,n,m)
      REAL(SP), DIMENSION(:,:), POINTER :: p, reallocate_rm
      INTEGER(I4B), INTENT(IN) :: n,m
      INTEGER(I4B) :: nold,mold,ierr
      allocate(reallocate_rm(n,m),stat=ierr)
      if (ierr /= 0) call &
            nrerror('reallocate_rm: problem in attempt to allocate memory')
      if (.not. associated(p)) RETURN
      nold=size(p,1)
      mold=size(p,2)
      reallocate_rm(1:min(nold,n),1:min(mold,m))=&
            p(1:min(nold,n),1:min(mold,m))
      deallocate(p)
      END FUNCTION reallocate_rm
!BL
      FUNCTION reallocate_im(p,n,m)
      INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
      INTEGER(I4B), INTENT(IN) :: n,m
      INTEGER(I4B) :: nold,mold,ierr
      allocate(reallocate_im(n,m),stat=ierr)
      if (ierr /= 0) call &
            nrerror('reallocate_im: problem in attempt to allocate memory')
      if (.not. associated(p)) RETURN
      nold=size(p,1)
      mold=size(p,2)
      reallocate_im(1:min(nold,n),1:min(mold,m))=&
            p(1:min(nold,n),1:min(mold,m))
      deallocate(p)
      END FUNCTION reallocate_im
!BL
      FUNCTION ifirstloc(mask)
      LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
      INTEGER(I4B) :: ifirstloc
      INTEGER(I4B), DIMENSION(1) :: loc
      loc=maxloc(merge(1,0,mask))
      ifirstloc=loc(1)
      if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
      END FUNCTION ifirstloc
!BL
      FUNCTION imaxloc_r(arr)
      REAL(SP), DIMENSION(:), INTENT(IN) :: arr
      INTEGER(I4B) :: imaxloc_r
      INTEGER(I4B), DIMENSION(1) :: imax
      imax=maxloc(arr(:))
      imaxloc_r=imax(1)
      END FUNCTION imaxloc_r
!BL
      FUNCTION imaxloc_i(iarr)
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
      INTEGER(I4B), DIMENSION(1) :: imax
      INTEGER(I4B) :: imaxloc_i
      imax=maxloc(iarr(:))
      imaxloc_i=imax(1)
      END FUNCTION imaxloc_i
!BL
      FUNCTION iminloc(arr)
      REAL(SP), DIMENSION(:), INTENT(IN) :: arr
      INTEGER(I4B), DIMENSION(1) :: imin
      INTEGER(I4B) :: iminloc
      imin=minloc(arr(:))
      iminloc=imin(1)
      END FUNCTION iminloc
!BL
      SUBROUTINE assert1(n1,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, INTENT(IN) :: n1
      if (.not. n1) then
            write (*,*) 'nrerror: an assertion failed with this tag:', &
                  string
            STOP 'program terminated by assert1'
      end if
      END SUBROUTINE assert1
!BL
      SUBROUTINE assert2(n1,n2,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, INTENT(IN) :: n1,n2
      if (.not. (n1 .and. n2)) then
            write (*,*) 'nrerror: an assertion failed with this tag:', &
                  string
            STOP 'program terminated by assert2'
      end if
      END SUBROUTINE assert2
!BL
      SUBROUTINE assert3(n1,n2,n3,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, INTENT(IN) :: n1,n2,n3
      if (.not. (n1 .and. n2 .and. n3)) then
            write (*,*) 'nrerror: an assertion failed with this tag:', &
                  string
            STOP 'program terminated by assert3'
      end if
      END SUBROUTINE assert3
!BL
      SUBROUTINE assert4(n1,n2,n3,n4,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, INTENT(IN) :: n1,n2,n3,n4
      if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
            write (*,*) 'nrerror: an assertion failed with this tag:', &
                  string
            STOP 'program terminated by assert4'
      end if
      END SUBROUTINE assert4
!BL
      SUBROUTINE assert_v(n,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, DIMENSION(:), INTENT(IN) :: n
      if (.not. all(n)) then
            write (*,*) 'nrerror: an assertion failed with this tag:', &
                  string
            STOP 'program terminated by assert_v'
      end if
      END SUBROUTINE assert_v
!BL
      FUNCTION assert_eq2(n1,n2,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, INTENT(IN) :: n1,n2
      INTEGER :: assert_eq2
      if (n1 == n2) then
            assert_eq2=n1
      else
            write (*,*) 'nrerror: an assert_eq failed with this tag:', &
                  string
            STOP 'program terminated by assert_eq2'
      end if
      END FUNCTION assert_eq2
!BL
      FUNCTION assert_eq3(n1,n2,n3,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, INTENT(IN) :: n1,n2,n3
      INTEGER :: assert_eq3
      if (n1 == n2 .and. n2 == n3) then
            assert_eq3=n1
      else
            write (*,*) 'nrerror: an assert_eq failed with this tag:', &
                  string
            STOP 'program terminated by assert_eq3'
      end if
      END FUNCTION assert_eq3
!BL
      FUNCTION assert_eq4(n1,n2,n3,n4,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, INTENT(IN) :: n1,n2,n3,n4
      INTEGER :: assert_eq4
      if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
            assert_eq4=n1
      else
            write (*,*) 'nrerror: an assert_eq failed with this tag:', &
                  string
            STOP 'program terminated by assert_eq4'
      end if
      END FUNCTION assert_eq4
!BL
      FUNCTION assert_eqn(nn,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, DIMENSION(:), INTENT(IN) :: nn
      INTEGER :: assert_eqn
      if (all(nn(2:) == nn(1))) then
            assert_eqn=nn(1)
      else
            write (*,*) 'nrerror: an assert_eq failed with this tag:', &
                  string
            STOP 'program terminated by assert_eqn'
      end if
      END FUNCTION assert_eqn
!BL
      SUBROUTINE nrerror(string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      write (*,*) 'nrerror: ',string
      STOP 'program terminated by nrerror'
      END SUBROUTINE nrerror
!BL
      FUNCTION arth_r(first,increment,n)
      REAL(SP), INTENT(IN) :: first,increment
      INTEGER(I4B), INTENT(IN) :: n
      REAL(SP), DIMENSION(n) :: arth_r
      INTEGER(I4B) :: k,k2
      REAL(SP) :: temp
      if (n > 0) arth_r(1)=first
      if (n <= NPAR_ARTH) then
            do k=2,n
                  arth_r(k)=arth_r(k-1)+increment
            end do
      else
            do k=2,NPAR2_ARTH
                  arth_r(k)=arth_r(k-1)+increment
            end do
            temp=increment*NPAR2_ARTH
            k=NPAR2_ARTH
            do
                  if (k >= n) exit
                  k2=k+k
                  arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
                  temp=temp+temp
                  k=k2
            end do
      end if
      END FUNCTION arth_r
!BL
      FUNCTION arth_d(first,increment,n)
      REAL(DP), INTENT(IN) :: first,increment
      INTEGER(I4B), INTENT(IN) :: n
      REAL(DP), DIMENSION(n) :: arth_d
      INTEGER(I4B) :: k,k2
      REAL(DP) :: temp
      if (n > 0) arth_d(1)=first
      if (n <= NPAR_ARTH) then
            do k=2,n
                  arth_d(k)=arth_d(k-1)+increment
            end do
      else
            do k=2,NPAR2_ARTH
                  arth_d(k)=arth_d(k-1)+increment
            end do
            temp=increment*NPAR2_ARTH
            k=NPAR2_ARTH
            do
                  if (k >= n) exit
                  k2=k+k
                  arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
                  temp=temp+temp
                  k=k2
            end do
      end if
      END FUNCTION arth_d
!BL
      FUNCTION arth_i(first,increment,n)
      INTEGER(I4B), INTENT(IN) :: first,increment,n
      INTEGER(I4B), DIMENSION(n) :: arth_i
      INTEGER(I4B) :: k,k2,temp
      if (n > 0) arth_i(1)=first
      if (n <= NPAR_ARTH) then
            do k=2,n
                  arth_i(k)=arth_i(k-1)+increment
            end do
      else
            do k=2,NPAR2_ARTH
                  arth_i(k)=arth_i(k-1)+increment
            end do
            temp=increment*NPAR2_ARTH
            k=NPAR2_ARTH
            do
                  if (k >= n) exit
                  k2=k+k
                  arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
                  temp=temp+temp
                  k=k2
            end do
      end if
      END FUNCTION arth_i
!BL
!BL
      FUNCTION geop_r(first,factor,n)
      REAL(SP), INTENT(IN) :: first,factor
      INTEGER(I4B), INTENT(IN) :: n
      REAL(SP), DIMENSION(n) :: geop_r
      INTEGER(I4B) :: k,k2
      REAL(SP) :: temp
      if (n > 0) geop_r(1)=first
      if (n <= NPAR_GEOP) then
            do k=2,n
                  geop_r(k)=geop_r(k-1)*factor
            end do
      else
            do k=2,NPAR2_GEOP
                  geop_r(k)=geop_r(k-1)*factor
            end do
            temp=factor**NPAR2_GEOP
            k=NPAR2_GEOP
            do
                  if (k >= n) exit
                  k2=k+k
                  geop_r(k+1:min(k2,n))=temp*geop_r(1:min(k,n-k))
                  temp=temp*temp
                  k=k2
            end do
      end if
      END FUNCTION geop_r
!BL
      FUNCTION geop_d(first,factor,n)
      REAL(DP), INTENT(IN) :: first,factor
      INTEGER(I4B), INTENT(IN) :: n
      REAL(DP), DIMENSION(n) :: geop_d
      INTEGER(I4B) :: k,k2
      REAL(DP) :: temp
      if (n > 0) geop_d(1)=first
      if (n <= NPAR_GEOP) then
            do k=2,n
                  geop_d(k)=geop_d(k-1)*factor
            end do
      else
            do k=2,NPAR2_GEOP
                  geop_d(k)=geop_d(k-1)*factor
            end do
            temp=factor**NPAR2_GEOP
            k=NPAR2_GEOP
            do
                  if (k >= n) exit
                  k2=k+k
                  geop_d(k+1:min(k2,n))=temp*geop_d(1:min(k,n-k))
                  temp=temp*temp
                  k=k2
            end do
      end if
      END FUNCTION geop_d
!BL
      FUNCTION geop_i(first,factor,n)
      INTEGER(I4B), INTENT(IN) :: first,factor,n
      INTEGER(I4B), DIMENSION(n) :: geop_i
      INTEGER(I4B) :: k,k2,temp
      if (n > 0) geop_i(1)=first
      if (n <= NPAR_GEOP) then
            do k=2,n
                  geop_i(k)=geop_i(k-1)*factor
            end do
      else
            do k=2,NPAR2_GEOP
                  geop_i(k)=geop_i(k-1)*factor
            end do
            temp=factor**NPAR2_GEOP
            k=NPAR2_GEOP
            do
                  if (k >= n) exit
                  k2=k+k
                  geop_i(k+1:min(k2,n))=temp*geop_i(1:min(k,n-k))
                  temp=temp*temp
                  k=k2
            end do
      end if
      END FUNCTION geop_i
!BL
      FUNCTION geop_c(first,factor,n)
      COMPLEX(SP), INTENT(IN) :: first,factor
      INTEGER(I4B), INTENT(IN) :: n
      COMPLEX(SP), DIMENSION(n) :: geop_c
      INTEGER(I4B) :: k,k2
      COMPLEX(SP) :: temp
      if (n > 0) geop_c(1)=first
      if (n <= NPAR_GEOP) then
            do k=2,n
                  geop_c(k)=geop_c(k-1)*factor
            end do
      else
            do k=2,NPAR2_GEOP
                  geop_c(k)=geop_c(k-1)*factor
            end do
            temp=factor**NPAR2_GEOP
            k=NPAR2_GEOP
            do
                  if (k >= n) exit
                  k2=k+k
                  geop_c(k+1:min(k2,n))=temp*geop_c(1:min(k,n-k))
                  temp=temp*temp
                  k=k2
            end do
      end if
      END FUNCTION geop_c
!BL
      FUNCTION geop_dv(first,factor,n)
      REAL(DP), DIMENSION(:), INTENT(IN) :: first,factor
      INTEGER(I4B), INTENT(IN) :: n
      REAL(DP), DIMENSION(size(first),n) :: geop_dv
      INTEGER(I4B) :: k,k2
      REAL(DP), DIMENSION(size(first)) :: temp
      if (n > 0) geop_dv(:,1)=first(:)
      if (n <= NPAR_GEOP) then
            do k=2,n
                  geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
            end do
      else
            do k=2,NPAR2_GEOP
                  geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
            end do
            temp=factor**NPAR2_GEOP
            k=NPAR2_GEOP
            do
                  if (k >= n) exit
                  k2=k+k
                  geop_dv(:,k+1:min(k2,n))=geop_dv(:,1:min(k,n-k))*&
                        spread(temp,2,size(geop_dv(:,1:min(k,n-k)),2))
                  temp=temp*temp
                  k=k2
            end do
      end if
      END FUNCTION geop_dv
!BL
!BL
      RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)
      REAL(SP), DIMENSION(:), INTENT(IN) :: arr
      REAL(SP), OPTIONAL, INTENT(IN) :: seed
      REAL(SP), DIMENSION(size(arr)) :: ans
      INTEGER(I4B) :: n,j
      REAL(SP) :: sd
      n=size(arr)
      if (n == 0_i4b) RETURN
      sd=0.0_sp
      if (present(seed)) sd=seed
      ans(1)=arr(1)+sd
      if (n < NPAR_CUMSUM) then
            do j=2,n
                  ans(j)=ans(j-1)+arr(j)
            end do
      else
            ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
            ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
      end if
      END FUNCTION cumsum_r
!BL
      RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
      INTEGER(I4B), OPTIONAL, INTENT(IN) :: seed
      INTEGER(I4B), DIMENSION(size(arr)) :: ans
      INTEGER(I4B) :: n,j,sd
      n=size(arr)
      if (n == 0_i4b) RETURN
      sd=0_i4b
      if (present(seed)) sd=seed
      ans(1)=arr(1)+sd
      if (n < NPAR_CUMSUM) then
            do j=2,n
                  ans(j)=ans(j-1)+arr(j)
            end do
      else
            ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
            ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
      end if
      END FUNCTION cumsum_i
!BL
!BL
      RECURSIVE FUNCTION cumprod(arr,seed) RESULT(ans)
      REAL(SP), DIMENSION(:), INTENT(IN) :: arr
      REAL(SP), OPTIONAL, INTENT(IN) :: seed
      REAL(SP), DIMENSION(size(arr)) :: ans
      INTEGER(I4B) :: n,j
      REAL(SP) :: sd
      n=size(arr)
      if (n == 0_i4b) RETURN
      sd=1.0_sp
      if (present(seed)) sd=seed
      ans(1)=arr(1)*sd
      if (n < NPAR_CUMPROD) then
            do j=2,n
                  ans(j)=ans(j-1)*arr(j)
            end do
      else
            ans(2:n:2)=cumprod(arr(2:n:2)*arr(1:n-1:2),sd)
            ans(3:n:2)=ans(2:n-1:2)*arr(3:n:2)
      end if
      END FUNCTION cumprod
!BL
!BL
      FUNCTION poly_rr(x,coeffs)
      REAL(SP), INTENT(IN) :: x
      REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
      REAL(SP) :: poly_rr
      REAL(SP) :: pow
      REAL(SP), DIMENSION(:), ALLOCATABLE :: vec
      INTEGER(I4B) :: i,n,nn
      n=size(coeffs)
      if (n <= 0) then
            poly_rr=0.0_sp
      else if (n < NPAR_POLY) then
            poly_rr=coeffs(n)
            do i=n-1,1,-1
                  poly_rr=x*poly_rr+coeffs(i)
            end do
      else
            allocate(vec(n+1))
            pow=x
            vec(1:n)=coeffs
            do
                  vec(n+1)=0.0_sp
                  nn=ishft(n+1,-1)
                  vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                  if (nn == 1) exit
                  pow=pow*pow
                  n=nn
            end do
            poly_rr=vec(1)
            deallocate(vec)
      end if
      END FUNCTION poly_rr
!BL
      FUNCTION poly_dd(x,coeffs)
      REAL(DP), INTENT(IN) :: x
      REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs
      REAL(DP) :: poly_dd
      REAL(DP) :: pow
      REAL(DP), DIMENSION(:), ALLOCATABLE :: vec
      INTEGER(I4B) :: i,n,nn
      n=size(coeffs)
      if (n <= 0) then
            poly_dd=0.0_dp
      else if (n < NPAR_POLY) then
            poly_dd=coeffs(n)
            do i=n-1,1,-1
                  poly_dd=x*poly_dd+coeffs(i)
            end do
      else
            allocate(vec(n+1))
            pow=x
            vec(1:n)=coeffs
            do
                  vec(n+1)=0.0_dp
                  nn=ishft(n+1,-1)
                  vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                  if (nn == 1) exit
                  pow=pow*pow
                  n=nn
            end do
            poly_dd=vec(1)
            deallocate(vec)
      end if
      END FUNCTION poly_dd
!BL
      FUNCTION poly_rc(x,coeffs)
      COMPLEX(SPC), INTENT(IN) :: x
      REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
      COMPLEX(SPC) :: poly_rc
      COMPLEX(SPC) :: pow
      COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
      INTEGER(I4B) :: i,n,nn
      n=size(coeffs)
      if (n <= 0) then
            poly_rc=0.0_sp
      else if (n < NPAR_POLY) then
            poly_rc=coeffs(n)
            do i=n-1,1,-1
                  poly_rc=x*poly_rc+coeffs(i)
            end do
      else
            allocate(vec(n+1))
            pow=x
            vec(1:n)=coeffs
            do
                  vec(n+1)=0.0_sp
                  nn=ishft(n+1,-1)
                  vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                  if (nn == 1) exit
                  pow=pow*pow
                  n=nn
            end do
            poly_rc=vec(1)
            deallocate(vec)
      end if
      END FUNCTION poly_rc
!BL
      FUNCTION poly_cc(x,coeffs)
      COMPLEX(SPC), INTENT(IN) :: x
      COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: coeffs
      COMPLEX(SPC) :: poly_cc
      COMPLEX(SPC) :: pow
      COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
      INTEGER(I4B) :: i,n,nn
      n=size(coeffs)
      if (n <= 0) then
            poly_cc=0.0_sp
      else if (n < NPAR_POLY) then
            poly_cc=coeffs(n)
            do i=n-1,1,-1
                  poly_cc=x*poly_cc+coeffs(i)
            end do
      else
            allocate(vec(n+1))
            pow=x
            vec(1:n)=coeffs
            do
                  vec(n+1)=0.0_sp
                  nn=ishft(n+1,-1)
                  vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                  if (nn == 1) exit
                  pow=pow*pow
                  n=nn
            end do
            poly_cc=vec(1)
            deallocate(vec)
      end if
      END FUNCTION poly_cc
!BL
      FUNCTION poly_rrv(x,coeffs)
      REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
      REAL(SP), DIMENSION(size(x)) :: poly_rrv
      INTEGER(I4B) :: i,n,m
      m=size(coeffs)
      n=size(x)
      if (m <= 0) then
            poly_rrv=0.0_sp
      else if (m < n .or. m < NPAR_POLY) then
            poly_rrv=coeffs(m)
            do i=m-1,1,-1
                  poly_rrv=x*poly_rrv+coeffs(i)
            end do
      else
            do i=1,n
                  poly_rrv(i)=poly_rr(x(i),coeffs)
            end do
      end if
      END FUNCTION poly_rrv
!BL
      FUNCTION poly_ddv(x,coeffs)
      REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
      REAL(DP), DIMENSION(size(x)) :: poly_ddv
      INTEGER(I4B) :: i,n,m
      m=size(coeffs)
      n=size(x)
      if (m <= 0) then
            poly_ddv=0.0_dp
      else if (m < n .or. m < NPAR_POLY) then
            poly_ddv=coeffs(m)
            do i=m-1,1,-1
                  poly_ddv=x*poly_ddv+coeffs(i)
            end do
      else
            do i=1,n
                  poly_ddv(i)=poly_dd(x(i),coeffs)
            end do
      end if
      END FUNCTION poly_ddv
!BL
      FUNCTION poly_msk_rrv(x,coeffs,mask)
      REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
      LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
      REAL(SP), DIMENSION(size(x)) :: poly_msk_rrv
      poly_msk_rrv=unpack(poly_rrv(pack(x,mask),coeffs),mask,0.0_sp)
      END FUNCTION poly_msk_rrv
!BL
      FUNCTION poly_msk_ddv(x,coeffs,mask)
      REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
      LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
      REAL(DP), DIMENSION(size(x)) :: poly_msk_ddv
      poly_msk_ddv=unpack(poly_ddv(pack(x,mask),coeffs),mask,0.0_dp)
      END FUNCTION poly_msk_ddv
!BL
!BL
      RECURSIVE FUNCTION poly_term_rr(a,b) RESULT(u)
      REAL(SP), DIMENSION(:), INTENT(IN) :: a
      REAL(SP), INTENT(IN) :: b
      REAL(SP), DIMENSION(size(a)) :: u
      INTEGER(I4B) :: n,j
      n=size(a)
      if (n <= 0) RETURN
      u(1)=a(1)
      if (n < NPAR_POLYTERM) then
            do j=2,n
                  u(j)=a(j)+b*u(j-1)
            end do
      else
            u(2:n:2)=poly_term_rr(a(2:n:2)+a(1:n-1:2)*b,b*b)
            u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
      end if
      END FUNCTION poly_term_rr
!BL
      RECURSIVE FUNCTION poly_term_cc(a,b) RESULT(u)
      COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
      COMPLEX(SPC), INTENT(IN) :: b
      COMPLEX(SPC), DIMENSION(size(a)) :: u
      INTEGER(I4B) :: n,j
      n=size(a)
      if (n <= 0) RETURN
      u(1)=a(1)
      if (n < NPAR_POLYTERM) then
            do j=2,n
                  u(j)=a(j)+b*u(j-1)
            end do
      else
            u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b)
            u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
      end if
      END FUNCTION poly_term_cc
!BL
!BL
      FUNCTION zroots_unity(n,nn)
      INTEGER(I4B), INTENT(IN) :: n,nn
      COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
      INTEGER(I4B) :: k
      REAL(SP) :: theta
      zroots_unity(1)=1.0
      theta=TWOPI/n
      k=1
      do
            if (k >= nn) exit
            zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),SPC)
            zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
                  zroots_unity(2:min(k,nn-k))
            k=2*k
      end do
      END FUNCTION zroots_unity
!BL
      FUNCTION outerprod_r(a,b)
      REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
      REAL(SP), DIMENSION(size(a),size(b)) :: outerprod_r
      outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
            spread(b,dim=1,ncopies=size(a))
      END FUNCTION outerprod_r
!BL
      FUNCTION outerprod_d(a,b)
      REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
      REAL(DP), DIMENSION(size(a),size(b)) :: outerprod_d
      outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
            spread(b,dim=1,ncopies=size(a))
      END FUNCTION outerprod_d
!BL
      FUNCTION outerdiv(a,b)
      REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
      REAL(SP), DIMENSION(size(a),size(b)) :: outerdiv
      outerdiv = spread(a,dim=2,ncopies=size(b)) / &
            spread(b,dim=1,ncopies=size(a))
      END FUNCTION outerdiv
!BL
      FUNCTION outersum(a,b)
      REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
      REAL(SP), DIMENSION(size(a),size(b)) :: outersum
      outersum = spread(a,dim=2,ncopies=size(b)) + &
            spread(b,dim=1,ncopies=size(a))
      END FUNCTION outersum
!BL
      FUNCTION outerdiff_r(a,b)
      REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
      REAL(SP), DIMENSION(size(a),size(b)) :: outerdiff_r
      outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
            spread(b,dim=1,ncopies=size(a))
      END FUNCTION outerdiff_r
!BL
      FUNCTION outerdiff_d(a,b)
      REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
      REAL(DP), DIMENSION(size(a),size(b)) :: outerdiff_d
      outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
            spread(b,dim=1,ncopies=size(a))
      END FUNCTION outerdiff_d
!BL
      FUNCTION outerdiff_i(a,b)
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
      INTEGER(I4B), DIMENSION(size(a),size(b)) :: outerdiff_i
      outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
            spread(b,dim=1,ncopies=size(a))
      END FUNCTION outerdiff_i
!BL
      FUNCTION outerand(a,b)
      LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: a,b
      LOGICAL(LGT), DIMENSION(size(a),size(b)) :: outerand
      outerand = spread(a,dim=2,ncopies=size(b)) .and. &
            spread(b,dim=1,ncopies=size(a))
      END FUNCTION outerand
!BL
      SUBROUTINE scatter_add_r(dest,source,dest_index)
      REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
      REAL(SP), DIMENSION(:), INTENT(IN) :: source
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
      INTEGER(I4B) :: m,n,j,i
      n=assert_eq2(size(source),size(dest_index),'scatter_add_r')
      m=size(dest)
      do j=1,n
            i=dest_index(j)
            if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
      end do
      END SUBROUTINE scatter_add_r
      SUBROUTINE scatter_add_d(dest,source,dest_index)
      REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
      REAL(DP), DIMENSION(:), INTENT(IN) :: source
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
      INTEGER(I4B) :: m,n,j,i
      n=assert_eq2(size(source),size(dest_index),'scatter_add_d')
      m=size(dest)
      do j=1,n
            i=dest_index(j)
            if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
      end do
      END SUBROUTINE scatter_add_d
      SUBROUTINE scatter_max_r(dest,source,dest_index)
      REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
      REAL(SP), DIMENSION(:), INTENT(IN) :: source
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
      INTEGER(I4B) :: m,n,j,i
      n=assert_eq2(size(source),size(dest_index),'scatter_max_r')
      m=size(dest)
      do j=1,n
            i=dest_index(j)
            if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
      end do
      END SUBROUTINE scatter_max_r
      SUBROUTINE scatter_max_d(dest,source,dest_index)
      REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
      REAL(DP), DIMENSION(:), INTENT(IN) :: source
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
      INTEGER(I4B) :: m,n,j,i
      n=assert_eq2(size(source),size(dest_index),'scatter_max_d')
      m=size(dest)
      do j=1,n
            i=dest_index(j)
            if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
      end do
      END SUBROUTINE scatter_max_d
!BL
      SUBROUTINE diagadd_rv(mat,diag)
      REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
      REAL(SP), DIMENSION(:), INTENT(IN) :: diag
      INTEGER(I4B) :: j,n
      n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
      do j=1,n
            mat(j,j)=mat(j,j)+diag(j)
      end do
      END SUBROUTINE diagadd_rv
!BL
      SUBROUTINE diagadd_r(mat,diag)
      REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
      REAL(SP), INTENT(IN) :: diag
      INTEGER(I4B) :: j,n
      n = min(size(mat,1),size(mat,2))
      do j=1,n
            mat(j,j)=mat(j,j)+diag
      end do
      END SUBROUTINE diagadd_r
!BL
      SUBROUTINE diagmult_rv(mat,diag)
      REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
      REAL(SP), DIMENSION(:), INTENT(IN) :: diag
      INTEGER(I4B) :: j,n
      n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
      do j=1,n
            mat(j,j)=mat(j,j)*diag(j)
      end do
      END SUBROUTINE diagmult_rv
!BL
      SUBROUTINE diagmult_r(mat,diag)
      REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
      REAL(SP), INTENT(IN) :: diag
      INTEGER(I4B) :: j,n
      n = min(size(mat,1),size(mat,2))
      do j=1,n
            mat(j,j)=mat(j,j)*diag
      end do
      END SUBROUTINE diagmult_r
!BL
      FUNCTION get_diag_rv(mat)
      REAL(SP), DIMENSION(:,:), INTENT(IN) :: mat
      REAL(SP), DIMENSION(size(mat,1)) :: get_diag_rv
      INTEGER(I4B) :: j
      j=assert_eq2(size(mat,1),size(mat,2),'get_diag_rv')
      do j=1,size(mat,1)
            get_diag_rv(j)=mat(j,j)
      end do
      END FUNCTION get_diag_rv
!BL
      FUNCTION get_diag_dv(mat)
      REAL(DP), DIMENSION(:,:), INTENT(IN) :: mat
      REAL(DP), DIMENSION(size(mat,1)) :: get_diag_dv
      INTEGER(I4B) :: j
      j=assert_eq2(size(mat,1),size(mat,2),'get_diag_dv')
      do j=1,size(mat,1)
            get_diag_dv(j)=mat(j,j)
      end do
      END FUNCTION get_diag_dv
!BL
      SUBROUTINE put_diag_rv(diagv,mat)
      REAL(SP), DIMENSION(:), INTENT(IN) :: diagv
      REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
      INTEGER(I4B) :: j,n
      n=assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
      do j=1,n
            mat(j,j)=diagv(j)
      end do
      END SUBROUTINE put_diag_rv
!BL
      SUBROUTINE put_diag_r(scal,mat)
      REAL(SP), INTENT(IN) :: scal
      REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
      INTEGER(I4B) :: j,n
      n = min(size(mat,1),size(mat,2))
      do j=1,n
            mat(j,j)=scal
      end do
      END SUBROUTINE put_diag_r
!BL
      SUBROUTINE unit_matrix(mat)
      REAL(SP), DIMENSION(:,:), INTENT(OUT) :: mat
      INTEGER(I4B) :: i,n
      n=min(size(mat,1),size(mat,2))
      mat(:,:)=0.0_sp
      do i=1,n
            mat(i,i)=1.0_sp
      end do
      END SUBROUTINE unit_matrix
!BL
      FUNCTION upper_triangle(j,k,extra)
      INTEGER(I4B), INTENT(IN) :: j,k
      INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
      LOGICAL(LGT), DIMENSION(j,k) :: upper_triangle
      INTEGER(I4B) :: n
      n=0
      if (present(extra)) n=extra
      upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
      END FUNCTION upper_triangle
!BL
      FUNCTION lower_triangle(j,k,extra)
      INTEGER(I4B), INTENT(IN) :: j,k
      INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
      LOGICAL(LGT), DIMENSION(j,k) :: lower_triangle
      INTEGER(I4B) :: n
      n=0
      if (present(extra)) n=extra
      lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
      END FUNCTION lower_triangle
!BL
      FUNCTION vabs(v)
      REAL(SP), DIMENSION(:), INTENT(IN) :: v
      REAL(SP) :: vabs
      vabs=sqrt(dot_product(v,v))
      END FUNCTION vabs
!BL
END MODULE nrutil
!===================================================================================

!===================================================================================
! CDF matching module
! This fortran module transfer one distribution to another using cdf mapping method.
!-----------------------------------------------------------------------------------
module cdftransfer
use NRTYPE
implicit none
private

public :: transfer,&
        get_F_from_data_normal,&
        get_F_from_data_weibul,&
        get_F_from_data_EVI,&
        get_data_from_F_normal,&
        get_data_from_F_weibul,&
        get_data_from_F_EVI
        
contains

subroutine transfer(src, src_cdf_ref, n1, quan, desti_cdf_ref, n2, desti, extype,tag)
!This is a general transfer function to transfer value (src) in one reference
!set (src_cdf_ref) to a value (desti) in another referecne set(desti_cdf_ref). 
!The query is src, quan is the quantile value used in the transfer
!desti is the output
!extype is the type of distribution used if the value falls outside of the reference set
use NRTYPE
implicit none
integer n1,n2
real(SP) :: src
real(SP),dimension(n1) :: src_cdf_ref
real(SP),dimension(n2) :: desti_cdf_ref
real(SP),dimension(n1) :: quan1
real(SP),dimension(n2) :: quan2
real(SP) :: quan
real(SP) :: desti
real(SP) :: X
real(SP) :: src_mean, src_sd, src_skew
real(SP) :: desti_mean, desti_sd, desti_skew
real(SP) :: dummy_adev, dummy_var, dummy_curt
character(len=4) :: extype
character(len=*),optional:: tag
character(len=12)::mytag
integer :: i
integer :: ihigh, ilow

real(SP) :: mean, sd

if(present(tag)) then
    mytag=tag
else
    mytag='mytag'
endif

!
!get the array size
!

!
!allocate space
!

!
!sort the src_cdf_ref and desti_cdf_ref
!

!-xyuan call sort(n1,src_cdf_ref)
!-xyuan call sort(n2,desti_cdf_ref)

!
!calculate mean and standard deviations
!
!    print*,'INSIDE TRANSFER'
!    print*,src_cdf_ref
    call moment(n1,src_cdf_ref,src_mean,dummy_adev,src_sd,dummy_var,src_skew,dummy_curt)
!    print*,desti_cdf_ref
    call moment(n2,desti_cdf_ref,desti_mean,dummy_adev,desti_sd,dummy_var,desti_skew,dummy_curt)

!
!initialize percentile vectors using Gringorten plotting position
!Gringorton gives exceedance probability, so 1-gringorton
!

do i = 1, n1
      quan1(i) = 1.-real(n1+1-i-0.44)/real(n1+0.12)
end do

do i = 1, n2
      quan2(i) = 1.-real(n2+1-i-0.44)/real(n2+0.12)
end do

!
!===============================================================
!find the quantile first
!

!
!if src is smaller than the smallest value in the reference
!
if(src < src_cdf_ref(1))then
!print*, 'finding quantile case 1', src, src_cdf_ref(1)

    if(extype .eq. 'norm')then
        quan = get_f_from_data_normal(src_mean,src_sd,src)
    else
!        print*,'-----------------',trim(mytag),'-------------------------'
!        print*, 'using get_f_from_data_weibul ', src, src_cdf_ref(1),src_mean,src_sd,src_skew
        quan = get_f_from_data_weibul(src_mean,src_sd,src_skew,src)
    endif

!if query falls above maximum value in referecne set
elseif(src > src_cdf_ref(n1))then
!print*, 'find quantile case 2', src, src_cdf_ref(n1)
    if(extype .eq. 'norm')then
        quan = get_F_from_data_normal(src_mean,src_sd,src)
    else
        quan = get_F_from_data_EVI(src_mean,src_sd,src)
    endif
else
!print*, 'finding quantile regular case'
    do i=1,n1
        if(src .le. src_cdf_ref(i))then
            ihigh=i
            exit
        endif
        ihigh=n1
    end do
    
    do i=n1,1,-1
        if(src .ge. src_cdf_ref(i))then
            ilow=i
            exit
        endif
        ilow=1
    end do
    
    if(src_cdf_ref(ihigh) .ne. src_cdf_ref(ilow))then
    X=(src_cdf_ref(ihigh)-src)/(src_cdf_ref(ihigh)-src_cdf_ref(ilow))
    quan = (1.-X)*quan1(ihigh)+X*quan1(ilow)
    else
    quan = (quan1(ihigh)+quan1(ilow))/2.
    endif
endif
!print*,'Value is:',src,' Quantile is: ',quan

!=============================================================================================
!after the quantile is obtained, we now need to find the values in desti_cdf_ref

if(quan < quan2(1))then
!print*, 'finding value case 1',quan, quan2(1)
    if(extype .eq. 'norm')then
        desti = get_data_from_F_normal(desti_mean,desti_sd,quan)
    else
        desti = get_data_from_F_weibul(desti_mean,desti_sd,desti_skew,quan)
    endif
elseif(quan > quan2(n2))then
!print*, 'finding value case 2',quan, quan2(n2)    
    if(extype .eq. 'norm')then
        desti = get_data_from_F_normal(desti_mean,desti_sd,quan)
    else
        desti = get_data_from_F_EVI(desti_mean,desti_sd,quan)
    endif
else
!print*,'finding value regular case'
    do i=1,n2
        if(quan .le. quan2(i))then
            ihigh=i
            exit
        endif
        ihigh=n2
    end do
    
    do i=n2,1,-1
        if(quan .ge. quan2(i))then
            ilow=i
            exit
        endif
        ilow=1
    end do

    if(ihigh .ne. ilow)then
    X=(quan2(ihigh)-quan)/(quan2(ihigh)-quan2(ilow))
    desti = (1.-X)*desti_cdf_ref(ihigh)+X*desti_cdf_ref(ilow)
    else
    desti = desti_cdf_ref(ihigh)
    endif
endif
!print*,'Quantile is:',quan,' Value is: ',desti
end subroutine


!#######################################################################################
!
!     function get_F_from_data_normal(mean, sd, x)
!
!########################################################################################

function get_F_from_data_normal(mean, sd, x)
! uses approximation from Handbook of Hydrology eq. 18.2.2 
use NRTYPE
implicit none

integer :: sign
real(SP) :: mean, sd, x
real(SP) :: z, F, get_F_from_data_normal

sign=1
z = (x-mean)/sd
if(z < 0.) then
    sign=-1
    z=sign*z
endif
  F = 1-0.5*exp(-1.*((83.*z+351.)*z+562.)/(165.+703./z))
!

if(F>1.0 .or. F<0.0)then
    print*,'Error in get_F_from_data_normal'
    print*,'Bad quantile: data=',x,' mean=',mean,' sd=',sd,' quant=',F
    stop
!    if(F>=1.0) F=0.99999
!    if(F<=0.0) F=0.00001
else
    if(sign == -1) F = 1-F
endif
get_F_from_data_normal  = F
end function get_F_from_data_normal


!#######################################################################################
!
!  function get_data_from_F_normal(mean,sd,F)
!
!########################################################################################

function get_data_from_F_normal(mean,sd,F)
use NRTYPE
implicit none
real(SP), intent(in) :: mean, sd, F
real(SP) :: get_data_from_F_normal, x,z
! /* uses approximation from Handbook of Hydrology eq. 18.2.3a */
  if(F > 1.0 .or. F < 0.0)then
    print*,'Error in get_data_from_F_normal'
    print*,'Bad quantile: mean=',mean,' sd=',sd,' quant=',F
    stop
  endif
  
! /*  z = ((F**0.135)-((1-F)**0.135))/0.1975;*/
  z = ((F**0.135)-((1.-F)**0.135))/0.1975;
  if(z > 5)then
    print*,'Pushing upper limit of applicability of equation'
    print*,'in get_data_from_F_normal best range is 0<z<5. z= ',z
  endif
  x = z*sd+mean
  get_data_from_F_normal=x
end function get_data_from_F_normal


!#######################################################################################
!
!  function get_F_from_data_EVI(mean, sd, x)
!
!########################################################################################

function get_F_from_data_EVI(mean, sd, x)
use NRTYPE
implicit none
!  Gumbel (EV Type I distribution) for maxima 
real(SP), intent(in) :: mean, sd, x
real(SP) :: a, b, F, get_F_from_data_EVI

 b = 3.14159/(sd*sqrt(6.))
 a = mean-0.5772/b
 F = exp(-1.*exp(-b*(x-a)))

!  added following adjustment so that outlier ensembles don't */
!  get totally booted out AWW-020903 */
if(F == 1.0) then
   F = F-.00238  !to adjust quantile if it is out of bounds 
   print*, 'Invoked TINY in get_F_from_data_EVI: data=',x,' mean=',mean,' sd=',sd,' quant=',F
elseif( F == 0.0) then
   F = F +.00238  !to adjust quantile if it is out of bounds 
   print*, 'Invoked TINY in get_F_from_data_EVI: data=',x,' mean=',mean,' sd=',sd,' quant=',F
endif

if(F>=1.0 .or. F<=0.0)then
   print*,'Error in get_F_from_data_EVI'
   print*,'Bad quantile: data=',x,' mean=',mean,' sd=',sd,' quant=',F
   stop
endif
get_F_from_data_EVI =F
end function  get_F_from_data_EVI

!#######################################################################################
!
!   function get_data_from_F_EVI(mean, sd, F)
!
!########################################################################################

function get_data_from_F_EVI(mean, sd, F)
use NRTYPE
implicit none
!   Gumbel (EV Type I distribution) for maxima */
real(SP), intent(in) :: mean, sd, F
real(SP) :: x, get_data_from_F_EVI
real(SP) :: a, b

  if(F>=1.0 .or. F<=0.0)then
    print*,'Error in get_data_from_F_EVI'
    print*,'Bad quantile: mean=',mean,' sd=',sd,' quant=',F
    stop
  endif
  
  b = 3.14159/(sd*sqrt(6.))
  a = mean-0.5772/b
  x = a-(1/b)*(log(-log(F)))
  get_data_from_F_EVI = x
end function get_data_from_F_EVI

!#######################################################################################
!
!  function get_F_from_data_weibul(mean, sd, skew, x)
!
!########################################################################################

function get_F_from_data_weibul(mean, sd, skew, x)
use NRTYPE
implicit none
real(SP),intent(in) :: mean, sd, skew, x
!   Weibull (EV Type III distribution) for minima, bounded at zero 
!   approximation for a (alpha) is eq. 11-32 from Kite (1977) 

real(SP) :: a,b,AA,BB,bound
real(SP) :: get_F_from_data_weibul, F

!  print*,'begin get_f_from_data_weibul'
!  print*,'INSIDE get_F_from_data_weibul, mean sd skew x',mean,sd,skew,x

  if(skew .gt. 8.214 .or. skew .lt. -1.)then
    print*,'Outside limit for table in get_F_from_data_weibul'
    print*,'best range is -1<skew<8.2 skew=',skew
  endif

  call weibul_params(skew,a,AA,BB)
  b = AA*sd+mean;
  bound = b-BB*sd;

! lower bound minimum of zero, but allow larger minima if data say so */
  if(bound < 0.) bound=0.
  if(bound > x ) bound=0.
if(x<0.)then
  F=0.00001
else
! F =1.-exp(-1.* ((x-bound)/(b-bound))**a)
  F =1.-exp(min(5., -1.* ((x-bound)/(b-bound))**a))  !-xyuan
endif
  

!  if(F >= 1.0 .or.  F <= 0.0) then
! based on AWW-042001: altered as follows to allow for 0 precip pass through
  if(F >= 1.0 .or.  F <= 0.0 ) then
!    print*, 'Error in get_F_from_data_weibul'
!    print*, 'Bad quantile: data=',x, ' mean=',mean, ' sd=',sd,' skew=',skew,' quant=',F
!    print*, 'a = ',a,' AA=',aa,' bb=',bb,' b=',b,' bound=',bound

   F=0.00001
   endif
   get_F_from_data_weibul = F
!   print*,'INSIDE get_F_from_data_weibul, F=',get_F_from_data_weibul
!   print*,'finish get_f_from_data_weibul'
end function get_F_from_data_weibul


!#######################################################################################
!
!   function get_data_from_F_weibul(mean,sd,skew,F)
!
!########################################################################################

function get_data_from_F_weibul(mean,sd,skew,F)
!  /* Weibull (EV Type III distribution) for minima, bounded at zero */
!  /* approximation for a (alpha) is eq. 11-32 from Kite (1977) */
use NRTYPE
implicit none
real(SP), intent(in) :: mean, sd, skew, F
real(SP) :: a, b, AA, BB, bound, x, get_data_from_F_weibul
  if(F >= 1.0 .or. F <= 0.0)then
    print*,'Error entering get_data_from_F_weibul'
    print*,'Bad quantile: mean=',mean,' sd=',sd,' skew=',skew,' quant=',F
    stop
  endif
  
  if(skew>8.2 .or. skew<-1)then
    print*,'Outside limit for table in get_data_from_F_weibul'
    print*,'best range is -1<skew<8.2 skew=',skew
  endif
  
!  a = 1/(0.2777757913+0.3132617714*skew+0.0575670910*pow(skew,2)-
!      0.0013038566*pow(skew,3)-0.0081523408*pow(skew,4)); 
  call weibul_params(skew,a,AA,BB)
  b = AA*sd+mean
  bound = b-BB*sd
!  /*  lower bound minimum of zero, but allow larger minima if data say so */
  if(bound<0) bound=0
  x = ((-log(1-F))**(1/a))*(b-bound) + bound
  get_data_from_F_weibul = x
end function get_data_from_F_weibul


!#######################################################################################
!
!   subroutine weibul_params(skew,a,aa,bb)
!
!########################################################################################

subroutine weibul_params(skew,a,aa,bb)
use NRTYPE
implicit none
!   returns alpha, Aalpha, and Balpha for the Weibull distrubution 
!   table taken from Statistical Methods in Hydrology, Haan, shown below 
real(SP) :: skew
real(SP), intent(out)  :: a
real(SP), intent(out)  :: aa
real(SP), intent(out)  :: bb

integer:: i
real(SP) :: X

real(SP),dimension(26) :: sk = (/-1.000, -0.971, -0.917, -0.867, -0.638, -0.254, 0.069, 0.359,&
		 0.631,  0.896,  1.160,  1.430,  1.708,  2.000, 2.309, 2.640, &
                 2.996,  3.382,  3.802,  4.262,  4.767,  5.323, 5.938, 6.619, &
 		 7.374,  8.214 /)
real(SP),dimension(26) :: inva = (/0.020, 0.030, 0.040, 0.050, 0.100, 0.200, 0.300, 0.400, 0.500, &
	          0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400, &
	          1.500, 1.600, 1.700, 1.800, 1.900, 2.000, 2.100, 2.200 /)
real(SP),dimension(26) :: Avec = (/0.446,  0.444,  0.442,  0.439,  0.425,  0.389,  0.346,  0.297, &
                  0.246,  0.193,  0.142,  0.092,  0.044,  0.000, -0.040, -0.077, &
                 -0.109, -0.136, -0.160, -0.180, -0.196, -0.208, -0.217, -0.224,&
		 -0.227, -0.229 /)
real(SP),dimension(26) :: Bvec = (/40.005, 26.987, 20.481, 16.576, 8.737, 4.755, 3.370, 2.634,&
		   2.159,  1.815,  1.549,  1.334, 1.154, 1.000, 0.867, 0.752, &
                   0.652,  0.563,  0.486,  0.418, 0.359, 0.308, 0.263, 0.224, &
		   0.190, 0.161 /)

if(skew>sk(26)) then
   skew=sk(26)
else if (skew<sk(1)) then
   skew=sk(1)
endif

do i=1,25
      if(skew <= sk(i+1) .and. skew >= sk(i))then
            if((sk(i+1)-sk(i)) == 0.) then
            X=1.
            else
            X=(sk(i+1)-skew)/(sk(i+1)-sk(i))
            a=1./(X*inva(i)+(1.-X)*inva(i+1))
            aa= X*Avec(i)+(1.-X)*Avec(i+1)
            bb= X*Bvec(i)+(1.-X)*Bvec(i+1)
            endif
      endif
end do
end subroutine  weibul_params   

end module
!-----------------------------------------------------------------------------------
! end of module cdftransfer
!===================================================================================

!===================================================================================
! The main program for regridding, bias correction and spi calculation
!-----------------------------------------------------------------------------------
program NMMEspi
use netcdf
IMPLICIT NONE
integer, parameter::clon=360,clat=181,dims=5,fcst=6,year=29,yclim=61,nparam=5,ens=47,nvar=2,nvar1=5
integer :: i,ii,j,jj,jjj,k,kk,kkk,ll,nkk
integer :: status,ncid,ncid1,varid(nvar1+3),varid1
integer :: start(dims),count(dims),dimids(dims)
integer lon_dimid, lat_dimid, fcst_dimid, ens_dimid, mod_dimid, year_dimid, bmon_dimid
integer cnt,tnm,ntag,ens1(nparam),spilead(4),yr1,nm1,day(12),nday
data ens1/6,10,11,10,10/,spilead/1,3,6,12/,day/31,28,31,30,31,30,31,31,30,31,30,31/
!data ens1/1,1,1,1,1/,spilead/1,3,6,12/,day/31,28,31,30,31,30,31,31,30,31,30,31/
character(len=*),parameter::model(nparam) = (/'COLA-RSMAS-CCSM3','GFDL-CM2p1-aer04','NASA-GMAO-062012','CMC1-CanCM3','CMC2-CanCM4'/), &
                            vcase(nparam) = (/'COLA-RSMAS-CCSM3','GFDL-CM2p1','NASA-GMAO','CMC1-CanCM3','CMC2-CanCM4'/), &
                            varname(nvar) = (/'PR','T2M'/), &
                            vardesc(nvar) = (/'total precipitation','2m temperature'/), &
                            varunit(nvar) = (/'mm/day','C'/), &
                            varcdf(nvar)  = (/'prec','tref'/), &
                            varname1(nvar1) = (/'SPI1','SPI3','SPI6','SPI12','T2ano'/), &
                            vardesc1(nvar1) = (/'1-month SPI','3-month SPI','6-month SPI','12-month SPI','monthly mean T2M anomaly'/), &
                            varunit1(nvar1) = (/'-','-','-','-','C'/)
character*2 cmon,cmon1,cens
character*4 cyr,cyr1,vtag,vtag1(2)
data vtag1/'prcp','norm'/
character*255 ext,ext0,ext1
real,allocatable::pr0(:,:,:,:),obsc0(:,:,:,:),var(:,:,:,:,:),var1(:,:,:,:,:),ovar3(:,:,:,:,:),ovar(:,:,:)
real,allocatable::lonw(:),latw(:),lcc(:,:)
real lonc(clon+1),latc(clat),tlon,tvar(fcst,ens,nvar)
real ww,we,ws,wn,tdis
real stdnorm_cdf(1001)
real pr(year,ens),pr1(ens),obsc(yclim)
real var2(31),tvar1
real alpha,beta,gamm,pzero,gamcdf,anvnrm,variance,lead(fcst)
integer wlon,wlat,yr,nm
real slon,slat,res
namelist /nmmepar/ wlon,wlat,slon,slat,res,yr,nm,ext,ext0

open(2,file='namelist')
read(2,nmmepar)
close(2)
allocate(lonw(wlon),latw(wlat),lcc(wlon,wlat))
write(cyr,'(i4)') yr
write(cmon,'(i2.2)') nm
call check( nf90_open(trim(ext)//'hindcast/landmask.nc', NF90_NOWRITE, ncid1) )
call check( nf90_inq_varid(ncid1, 'PR', varid1) )
call check( nf90_get_var(ncid1, varid1, lcc) )
call check( nf90_close(ncid1) )
do i=1,wlon
   lonw(i)=slon+res*(i-1)
end do
do j=1,wlat
   latw(j)=slat+res*(j-1)
end do
do i=1,clon+1
   lonc(i)=i-1
end do
do j=1,clat
   latc(j)=-90.+(j-1)
end do
do i=1,fcst
   lead(i)=i-0.5
end do
!====================================================================
! $1. regrid the raw forecast data
!====================================================================
allocate(var(wlon,wlat,fcst,ens,nvar),var1(clon+1,clat,fcst,ens,nvar))
nkk=0
do kk=1,nparam
   do jjj=1,ens1(kk)
      nkk=nkk+1
      write(cens,'(i2)') jjj
      do ll=1,nvar
         call check( nf90_open(trim(ext)//cyr//cmon//'/'//trim(varcdf(ll))//'/'//trim(model(kk))//'_'//trim(adjustl(cens))//'.nc', NF90_NOWRITE, ncid1) )
         call check( nf90_inq_varid(ncid1, varcdf(ll), varid1) )
         if(ll==2 .and. kk<=3)then
            start(1:4)=(/1,1,1,1/); count(1:4)=(/clon,clat,1,fcst/)
            call check( nf90_get_var(ncid1, varid1, var1(1:clon,:,:,nkk,ll), start(1:4), count(1:4)) )
         else
           start(1:3)=(/1,1,1/); count(1:3)=(/clon,clat,fcst/)
           call check( nf90_get_var(ncid1, varid1, var1(1:clon,:,:,nkk,ll), start(1:3), count(1:3)) )
         endif
         call check( nf90_close(ncid1) )
      end do
   end do
end do
var1(clon+1,:,:,:,:)=var1(1,:,:,:,:)
var=-9.99e+08
do i=1,wlon
   tlon=lonw(i)
   if(tlon<0.) tlon=tlon+360.
   do ii=1,clon
      if(tlon>=lonc(ii) .and. tlon<lonc(ii+1))then
         ww=(tlon-lonc(ii))/(lonc(ii+1)-lonc(ii))
         we=1.-ww
         exit
      endif
   end do
   do j=1,wlat
      if(lcc(i,j)>=0.5)then
         do jj=1,clat-1
            if(latw(j)>=latc(jj) .and. latw(j)<latc(jj+1))then
               ws=(latw(j)-latc(jj))/(latc(jj+1)-latc(jj))
               wn=1.-ws
               exit
            endif
         end do
         tvar=0.; tdis=0.
         if(var1(ii,jj,1,1,1)>-100.)then
            tdis=tdis+we*wn
            tvar=tvar+var1(ii,jj,:,:,:)*we*wn
         endif
         if(var1(ii+1,jj,1,1,1)>-100.)then
            tdis=tdis+ww*wn
            tvar=tvar+var1(ii+1,jj,:,:,:)*ww*wn
         endif
         if(var1(ii,jj+1,1,1,1)>-100.)then
            tdis=tdis+we*ws
            tvar=tvar+var1(ii,jj+1,:,:,:)*we*ws
         endif
         if(var1(ii+1,jj+1,1,1,1)>-100.)then
            tdis=tdis+ww*ws
            tvar=tvar+var1(ii+1,jj+1,:,:,:)*ww*ws
         endif
         if(tdis>0.)then
            var(i,j,:,:,:)=tvar/tdis
            var(i,j,:,:,2)=var(i,j,:,:,2)-273.16
         endif
      endif
   end do
end do
deallocate(var1)

!====================================================================
! $2. Bias correction for PR and T2M
!====================================================================
allocate(pr0(wlon,wlat,year,ens))
allocate(obsc0(wlon,wlat,fcst,yclim))
allocate(var1(wlon,wlat,fcst,ens,nvar))
write(ext1,'(a,a)') trim(ext),'hindcast/stdnorm_sample.dat'
call read_stdnorm(ext1,stdnorm_cdf)
call sort(1001,stdnorm_cdf)
do ntag=1,nvar
   vtag=vtag1(ntag)
!==============================
! read observation climatology
!==============================
   call check( nf90_open(trim(ext)//'hindcast/OBS.nc', NF90_NOWRITE, ncid) )
   if(vtag=='prcp')then
      call check( nf90_inq_varid(ncid, 'PR', varid1) )
   elseif(vtag=='norm')then
      call check( nf90_inq_varid(ncid, 'T2M', varid1) )
   endif
   do k=1,fcst
      tnm=nm+k-1
      if(tnm>12) tnm=tnm-12
      start(1:4)=(/1,1,1,tnm/);  count(1:4)=(/wlon,wlat,yclim,1/)
      call check( nf90_get_var(ncid, varid1, obsc0(:,:,k,:), start(1:4), count(1:4)) )
   end do
   call check( nf90_close(ncid) )
   do k=1,fcst
!==============================
! read hindcast
!==============================
      nkk=0
      do kk=1,nparam
         do kkk=1,ens1(kk)
            nkk=nkk+1
            write(cens,'(i2.2)') kkk
            call check( nf90_open(trim(ext)//'hindcast/'//trim(vcase(kk))//cens//'.nc', NF90_NOWRITE, ncid1) )
            if(vtag=='prcp')then
               call check( nf90_inq_varid(ncid1, 'PR', varid1) )
            elseif(vtag=='norm')then
               call check( nf90_inq_varid(ncid1, 'T2M', varid1) )
            endif
            start=(/1,1,k,1,nm/); count=(/wlon,wlat,1,year,1/)
            call check( nf90_get_var(ncid1, varid1, pr0(:,:,:,nkk), start, count) )
            call check( nf90_close(ncid1) )
         end do
      end do
      if(vtag=='norm')then
         pr0=pr0-273.16
      endif
!==============================
! bias correction
!==============================
      do i=1,wlon
         do j=1,wlat
            if(obsc0(i,j,k,1)>-100. .and. lcc(i,j)>=0.5)then
               obsc=obsc0(i,j,k,:)
               pr=pr0(i,j,:,:)
               pr1=var(i,j,k,:,ntag)
               cnt=0
               call BiasCorrect(yclim,year,nparam,ens,ens1,obsc,pr,pr1,var1(i,j,k,:,ntag),stdnorm_cdf,vtag)
            else
               var1(i,j,k,:,ntag)=-9.99e+08
            endif
         end do
      end do
   end do !end of fcst lead
end do !end of ntag
deallocate(pr0,obsc0,var)

!==============================
! output bias corrected results
!==============================
nkk=0
do kk=1,nparam
   do jjj=1,ens1(kk)
      nkk=nkk+1
      write(cens,'(i2)') jjj
      call check( nf90_create(trim(ext)//cyr//cmon//'/'//trim(model(kk))//'_'//trim(adjustl(cens))//'_BC.nc', nf90_netcdf4, ncid) )
      call check( nf90_def_dim(ncid, 'lon',  wlon, lon_dimid) )
      call check( nf90_def_dim(ncid, 'lat',  wlat, lat_dimid) )
      call check( nf90_def_dim(ncid, 'fcst', fcst, fcst_dimid) )
      dimids(1:3) = (/ lon_dimid, lat_dimid, fcst_dimid /)
      do i=1,nvar
         call check( nf90_def_var(ncid, varname(i), NF90_REAL,  dimids(1:3), varid(i), deflate_level = 1) )
         call check( nf90_put_att(ncid, varid(i), 'long_name',  vardesc(i)) )
         call check( nf90_put_att(ncid, varid(i), 'units',      varunit(i)) )
         call check( nf90_put_att(ncid, varid(i), '_FillValue',    -9.99e+08) )
         call check( nf90_put_att(ncid, varid(i), 'missing_value', -9.99e+08) )
      end do
      dimids(1) = lon_dimid
      call check( nf90_def_var(ncid, 'lon', NF90_REAL, dimids(1), varid(nvar+1)) )
      call check( nf90_put_att(ncid, varid(nvar+1), 'units', 'degrees_east') )
      dimids(1) = lat_dimid
      call check( nf90_def_var(ncid, 'lat', NF90_REAL, dimids(1), varid(nvar+2)) )
      call check( nf90_put_att(ncid, varid(nvar+2), 'units', 'degrees_north') )
      dimids(1) = fcst_dimid
      call check( nf90_def_var(ncid, 'fcst', NF90_REAL, dimids(1), varid(nvar+3)) )
      call check( nf90_put_att(ncid, varid(nvar+3), 'units', 'months') )
      call check( nf90_put_att(ncid, nf90_global, 'period',  cyr//cmon) )
      call check( nf90_put_att(ncid, nf90_global, 'source',  'bias corrected NMME realtime forecast') )
      call check( nf90_put_att(ncid, nf90_global, 'contact', 'Dr. Xing Yuan, xingy@princeton.edu, 6092583967') )
      call check( nf90_put_att(ncid, nf90_global, 'history', 'Created by Dr. Xing Yuan, Princeton, 2013-08-25') )
      call check( nf90_enddef(ncid) )
      do i=1,nvar
         call check( nf90_put_var(ncid, varid(i), var1(:,:,:,nkk,i)) )
      end do
      call check( nf90_put_var(ncid, varid(nvar+1), lonw) )
      call check( nf90_put_var(ncid, varid(nvar+2), latw) )
      call check( nf90_put_var(ncid, varid(nvar+3), lead) )
      call check( nf90_close(ncid) )
   end do
end do

!====================================================================
! $3. Calculate 1,3,6,12-month SPI drought indices
!====================================================================
allocate(var(wlon,wlat,fcst,ens,nvar1),ovar3(wlon,wlat,fcst,4,30),ovar(wlon,wlat,11))
call check( nf90_open(trim(ext)//'hindcast/PR-ma.nc', NF90_NOWRITE, ncid) )
call check( nf90_inq_varid(ncid, 'PR', varid1) )
do k=1,fcst
   tnm=nm+k-1
   if(tnm>12) tnm=tnm-12
   start=(/1,1,1,32,tnm/);  count=(/wlon,wlat,4,30,1/)
   call check( nf90_get_var(ncid, varid1, ovar3(:,:,k,:,:), start, count) )
end do
call check( nf90_close(ncid) )

do k=1,11
   yr1=yr; nm1=nm-k
   if(nm1<1)then
      nm1=nm1+12
      yr1=yr1-1
   endif
   write(cyr1,'(i4)') yr1
   write(cmon1,'(i2.2)') nm1
   print*,yr1,nm1
   call check( nf90_open(trim(ext0)//'3B42RT_BC_'//cyr1//cmon1//'_monthly_0.250deg.nc', NF90_NOWRITE, ncid) )
   call check( nf90_inq_varid(ncid, 'prec', varid1) )
   start(1:3)=(/1,1,1/);  count(1:3)=(/wlon,wlat,1/)
   call check( nf90_get_var(ncid, varid1, ovar(:,:,k), start(1:3), count(1:3)) )
   call check( nf90_close(ncid) )
   nday=day(nm1)
   if(mod(yr1,4)==0 .and. nm1==2) nday=nday+1
   ovar(:,:,k)=ovar(:,:,k)/nday
end do

var(:,:,:,:,1:4)=-9.99e+08
do i=1,wlon
   do j=1,wlat
      if(lcc(i,j)>=0.5 .and. ovar3(i,j,1,1,1)>-100.)then
         do k=1,fcst
            do ii=1,4
               var2(1:30)=ovar3(i,j,k,ii,:)
               if(variance(30,var2(1:30))>1.e-3)then
                  call gamfit(var2(1:30),30,alpha,beta,gamm,pzero)
                  do kk=1,ens
                     var2(31)=0.
                     if(k<spilead(ii))then
                        do kkk=1,spilead(ii)-k
                           !var2(31)=var2(31)+ovar(i,j,12-kkk)
                           var2(31)=var2(31)+ovar(i,j,kkk)
                        end do
                        do kkk=1,k
                           var2(31)=var2(31)+var1(i,j,k-kkk+1,kk,1)
                        end do 
                     else
                        do kkk=1,spilead(ii)
                           var2(31)=var2(31)+var1(i,j,k-kkk+1,kk,1)
                        end do
                     endif
                     var2(31)=var2(31)/spilead(ii)
                     if(var1(i,j,k,kk,1)>-100.)then
                        tvar1=max(1.e-15, gamcdf(beta,gamm,pzero,var2(31)))
                        var(i,j,k,kk,ii)=min(100.,anvnrm(tvar1))
                     endif
                  end do
               endif
            end do
         end do
      endif
   end do
end do

!====================================================================
! $4. Calculate monthly T2 anomaly
!====================================================================
call check( nf90_open(trim(ext)//'hindcast/OBSy.nc', NF90_NOWRITE, ncid) )
call check( nf90_inq_varid(ncid, 'T2M', varid1) )
do k=1,fcst
   tnm=nm+k-1
   if(tnm>12) tnm=tnm-12
   start(1:3)=(/1,1,tnm/);  count(1:3)=(/wlon,wlat,1/)
   call check( nf90_get_var(ncid, varid1, ovar(:,:,k), start(1:3), count(1:3)) )
end do
call check( nf90_close(ncid) )

var(:,:,:,:,nvar1)=-9.99e+08
do i=1,wlon
   do j=1,wlat
      if(lcc(i,j)>=0.5 .and. ovar(i,j,1)>-100.)then
         do k=1,fcst
            do kk=1,ens
               if(var1(i,j,k,kk,2)>-100.) var(i,j,k,kk,nvar1)=var1(i,j,k,kk,2)-ovar(i,j,k)
            end do
         end do
      endif
   end do
end do

!====================================================================
! $5. Output 1,3,6,12-month SPI and monthly T2 anomaly
!====================================================================
nkk=0
do kk=1,nparam
   do jjj=1,ens1(kk)
      nkk=nkk+1
      write(cens,'(i2)') jjj
      call check( nf90_create(trim(ext)//cyr//cmon//'/'//trim(model(kk))//'_'//trim(adjustl(cens))//'_sf.nc', nf90_netcdf4, ncid) )
      call check( nf90_def_dim(ncid, 'lon',  wlon, lon_dimid) )
      call check( nf90_def_dim(ncid, 'lat',  wlat, lat_dimid) )
      call check( nf90_def_dim(ncid, 'fcst', fcst, fcst_dimid) )
      dimids(1:3) = (/ lon_dimid, lat_dimid, fcst_dimid /)
      do i=1,nvar1
         call check( nf90_def_var(ncid, varname1(i), NF90_REAL,  dimids(1:3), varid(i), deflate_level = 1) )
         call check( nf90_put_att(ncid, varid(i), 'long_name',  varname1(i)) )
         call check( nf90_put_att(ncid, varid(i), 'units',      varunit1(i)) )
         call check( nf90_put_att(ncid, varid(i), '_FillValue',    -9.99e+08) )
         call check( nf90_put_att(ncid, varid(i), 'missing_value', -9.99e+08) )
      end do
      dimids(1) = lon_dimid
      call check( nf90_def_var(ncid, 'lon', NF90_REAL, dimids(1), varid(nvar1+1)) )
      call check( nf90_put_att(ncid, varid(nvar1+1), 'units', 'degrees_east') )
      dimids(1) = lat_dimid
      call check( nf90_def_var(ncid, 'lat', NF90_REAL, dimids(1), varid(nvar1+2)) )
      call check( nf90_put_att(ncid, varid(nvar1+2), 'units', 'degrees_north') )
      dimids(1) = fcst_dimid
      call check( nf90_def_var(ncid, 'fcst', NF90_REAL, dimids(1), varid(nvar1+3)) )
      call check( nf90_put_att(ncid, varid(nvar1+3), 'units', 'months') )
      call check( nf90_put_att(ncid, nf90_global, 'period',  cyr//cmon) )
      call check( nf90_put_att(ncid, nf90_global, 'source',  '1,3,6,12-month SPI and monthly T2 anomaly forecast') )
      call check( nf90_put_att(ncid, nf90_global, 'contact', 'Dr. Xing Yuan, xingy@princeton.edu, 6092583967') )
      call check( nf90_put_att(ncid, nf90_global, 'history', 'Created by Dr. Xing Yuan, Princeton, 2013-08-25') )
      call check( nf90_enddef(ncid) )
      do i=1,nvar1
         call check( nf90_put_var(ncid, varid(i), var(:,:,:,nkk,i)) )
      end do
      call check( nf90_put_var(ncid, varid(nvar1+1), lonw) )
      call check( nf90_put_var(ncid, varid(nvar1+2), latw) )
      call check( nf90_put_var(ncid, varid(nvar1+3), lead) )
      call check( nf90_close(ncid) )
   end do
end do

deallocate(lonw,latw,lcc,ovar3,ovar,var,var1)
!====================================================================
contains
  subroutine check(status)
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then
     print *, trim(nf90_strerror(status))
     stop "Stopped"
  end if
  end subroutine check
END

subroutine rdplot(dat,no,io,nrec,finame)
integer no,io,nrec
real  dat(no)
character(len=*)::  finame
open (io,status='old',file=trim(finame)&
      ,form='unformatted',access='direct',recl=no*1)
read (io,rec=nrec) (dat(i),i=1,no)
close(io)
end

subroutine read_stdnorm(ext,stdnorm_cdf)
implicit none
real::stdnorm_cdf(1001)
integer :: i
character*255 ext
!from methworld.wolfram.com
!the distribution of standard normal is
!
!D(x)= 1/2 [erf(x/sqrt(2))+1]
!
!The following data was generated with rnorm in R.

open(99,file=trim(ext),status='old')
do i=1,1001
    read(99,*),stdnorm_cdf(i)
end do
close(99)
end subroutine read_stdnorm

subroutine sort(n,a)
implicit none
integer i,j,n
real a(n),c
do i=1,n-1
   do j=i+1,n
      if(a(i)>a(j))then
         c=a(j)
         a(j)=a(i)
         a(i)=c
      endif
   end do
end do
end

SUBROUTINE moment(n,data,ave,adev,sdev,var,skew,curt)
USE nrtype
USE nrutil,only:nrerror
IMPLICIT NONE
REAL(SP), INTENT(OUT) :: ave,adev,sdev,var,skew,curt
INTEGER(I4B) :: n
REAL(SP), DIMENSION(n), INTENT(IN) :: data
REAL(SP) :: ep
REAL(SP), DIMENSION(size(data)) :: p,s
if (n <= 1) call nrerror('moment: n must be at least 2')
ave=sum(data(:))/n
s(:)=data(:)-ave
ep=sum(s(:))
adev=sum(abs(s(:)))/n
p(:)=s(:)*s(:)
var=sum(p(:))
p(:)=p(:)*s(:)
skew=sum(p(:))
p(:)=p(:)*s(:)
curt=sum(p(:))
var=(var-ep**2/n)/(n-1)
sdev=sqrt(var)
if (var /= 0.0) then
      skew=skew/(n*sdev**3)
      curt=curt/(n*var**2)-3.0_sp
else
      !call nrerror('moment: no skew or kurtosis when zero variance')
      !modified by Lifeng
      !print*,'moment: no skew or kurtosis when zero variance'
      skew=0
      curt=0
end if
END SUBROUTINE moment

subroutine BiasCorrect(yclim,n,nm,nparam,ens1,obsc,pr,pr1,var,stdnorm_cdf,vtag)
!=========================================================
! The bias correction program using quantile mapping.
!
! Author: Dr. Xing Yuan @ Princeton Univerity, 2011-10-02
! Revised: Xing Yuan, 2013-08-06
!
!=========================================================
use cdftransfer,only:transfer
implicit none
integer,intent(in)::nm,       &! number of model
                    nparam,   &! total number of ensemble
                    yclim,    &! length of observation
                    n          ! length of hindcast
integer,intent(in)::ens1(nm)   ! ensemble for each model
real,intent(in):: &
    stdnorm_cdf(1001),        &! standard normal distribution cdf
    obsc(yclim),              &! observation climatology
    pr(n,nparam),             &! hindcast
    pr1(nparam)                ! current forecast
real,intent(out)::var(nparam)  ! downscaled forecast
real,dimension(yclim)::obsc_cdf
real,dimension((n+1)*24,nm)::model_cdf
real :: quantile
integer i,j,k,is
character*4 vtag

obsc_cdf=obsc
call sort(yclim,obsc_cdf)
if(vtag=='prcp')then
   if(sum(obsc_cdf)<6.3)then
      var=0.
      goto 100
   endif
endif

is=0
do k=1,nm
   do i=1,ens1(k)
      is=is+1
      model_cdf((i-1)*(n+1)+1:i*(n+1)-1,k)=pr(:,is)
      model_cdf(i*(n+1),k)=pr1(is)
   end do
   call sort((n+1)*ens1(k),model_cdf(1:(n+1)*ens1(k),k))
end do

is=0
do k=1,nm
   do j=1,ens1(k)
      is=is+1
      call transfer(pr1(is),model_cdf(1:(n+1)*ens1(k),k),(n+1)*ens1(k),quantile,obsc_cdf,yclim,var(is),vtag,'pcptag1')
   end do
end do

100 continue
end

!=====================================================================
!
!  input prob; return z.
!
!=====================================================================
real function anvnrm(prob)
implicit none
real,parameter::c0=2.515517,c1=0.802853,c2=0.010328
real,parameter::d1=1.432788,d2=0.189269,d3=0.001308
real prob,sign,t
  
if (prob .gt. 0.5) then
   sign = 1.0
   prob = 1.0 - prob
else
   sign = -1.0
endif
  
if (prob .lt. 0.0) then
   write(*, *) 'Error in anvnrm(). Prob. not in [0,1.0]'
   anvnrm = 0.0
   return
endif

if (prob .eq. 0.0) then
   anvnrm = 1.0e37 * sign
   return
endif
  
t = sqrt(alog (1.0 / (prob * prob)))
anvnrm = (sign * (t - ((((c2 * t) + c1) * t) + c0) /  &
         ((((((d3 * t) + d2) * t) + d1) * t) + 1.0)))
return
end

!=====================================================================
!
!  Estimate incomplete gamma parameters.
!
!  Input:
!     datarr - data array
!     n - size of datarr
!
! Output:
!     alpha, beta, gamma - gamma paarameters
!     pzero - probability of zero.
!
!=====================================================================
subroutine gamfit(datarr, n, alpha, beta, gamm, pzero)
implicit none
integer i,n,nact
real alpha, beta, gamm, pzero
real datarr(n)
real sum,sumlog,av

if (n .le. 0) then
   write(0, *) 'Error in gamfit - empty data array'
   stop
endif

sum = 0.0
sumlog = 0.0
pzero = 0.0
nact = 0
  
!  compute sums
do i = 1, n
   if (datarr(i) .gt. 0.0) then
      sum = sum + datarr(i)
      sumlog = sumlog + alog (datarr(i))
      nact = nact + 1
   else
      pzero = pzero + 1
   endif
end do
pzero = pzero / n
if(nact .ne. 0) av = sum / nact
  
!  Bogus data array but do something reasonable
if(nact .eq. 1) then
   alpha = 0.0
   gamm = 1.0
   beta = av
   return
endif

!  They were all zeroes. 
if(pzero .eq. 1.0) then
   alpha = 0.0
   gamm = 1.0
   beta = av
   return
endif

!  Use MLE
alpha = alog (av) - sumlog / nact
if(abs(alpha)<1.e-15)then
   alpha = 0.0
   gamm = 1.0
   beta = av
   return
endif  !+xyuan
gamm = (1.0 + sqrt (1.0 + 4.0 * alpha / 3.0)) / (4.0 * alpha)
beta = av / gamm
  
return
end


!=====================================================================
!
!  Compute probability of a<=x using incomplete gamma parameters.
!
!  Input:
!      beta, gamma - gamma parameters
!      pzero - probability of zero.
!      x - value.
!
!  Return:
!      Probability  a<=x.
!
!=====================================================================

real function gamcdf(beta, gamm, pzero, x)
implicit none
real beta, gamm, pzero, x, gammap

if(x .le. 0.0) then
   gamcdf = pzero
else
!if(gamm<0.)then
!print*,'gammap(gamm, x / beta)',gamm, x / beta
!endif!xyuan
   gamcdf = pzero + (1.0 - pzero) * gammap(gamm, x / beta)
endif
return
end

!=========================================================================
!
!  Compute inverse gamma function i.e. return x given p where CDF(x) = p.
!
!  Input:
!      beta, gamma - gamma parameters
!      pzero - probability of zero.
!      prob - probability.
!
!  Return:
!      x as above.
!
!  Method:
!      We use a simple binary search to first bracket out initial
!      guess and then to refine our guesses until two guesses are within
!      tolerance (eps).  Is there a better way to do this?
!
!=========================================================================
real function gaminv(beta, gamm, pzero, prob)
implicit none
real,parameter::eps=1.0e-7
integer niter
real beta, gamm, pzero, prob, thigh, phigh, gamcdf, tlow, t, p
  
!  Check if prob < prob of zero
if (prob .le. pzero) then
   gaminv = 0.0
   return
endif
  
!     Otherwise adjust prob
prob = (prob - pzero) / (1.0 - pzero)
  
!  Make initial guess. Keep doubling estimate until prob is
!  bracketed.
thigh = 2.0*eps 
10 continue
   phigh = gamcdf (beta, gamm, pzero, thigh)
   if(phigh .ge. prob) goto 20
   thigh = thigh*2.0
   goto 10
20 continue
tlow = thigh / 2.0
  
!  Iterate to find root.
   niter = 0
30 continue
if((thigh - tlow) .le. eps) goto 40
niter = niter + 1
t = (tlow + thigh) / 2.0
p = gamcdf (beta, gamm, pzero, t)
      
if (p .lt. prob) then
   tlow = t
else
   thigh = t
endif
goto 30
40 continue
gaminv = (tlow + thigh) / 2.0
return
end

!=========================================================================
!
!  Functions for the incomplete gamma functions P and Q
!
!                  1     /x  -t a-1
!   P (a, x) = -------- |   e  t    dt,  a > 0
!              Gamma(x)/ 0
!
!   Q (a, x) = 1 - P (a, x)
!
!
!=========================================================================
!
! Evaluate P(a,x) by its series representation.  
!
real function gamser (a, x)
!  Maximum number of iterations, and bound on error.
implicit none
integer,parameter::maxitr=100
real,parameter::eps=3.0e-7
integer n
real gln, gammln, a, x, ap, sum, del

gln = gammln (a)
if (x .eq. 0.0) then
   gamser = 0.0
   return
endif
ap = a
sum = 1.0 / a
del = sum
  
do n = 1, maxitr
   ap = ap + 1.0
   del = del * (x / ap)
   sum = sum + del
   if (abs (del) .lt. eps * abs (sum)) exit
end do
gamser =  sum * exp (-x + a * alog (x) - gln)
return
end

!
!     Evaluate P(a,x) in its continued fraction representation.
!
real function gammcf (a, x)
implicit none
integer,parameter::maxitr=100
real,parameter::eps=3.0e-7
integer n
real  g, gln, gammln, a, x
real gold,a0,a1,b0,b1,fac,an,ana,anf
  
g = 0.0
gln = gammln (a)
gold = 0.0
a0 = 1.0
a1 = x
b0 = 0.0
b1 = 1.0
fac = 1.0
do n = 1, maxitr
   an = n
   ana = an - a
   a0 = (a1 + a0 * ana) * fac
   b0 = (b1 + b0 * ana) * fac
   anf = an * fac
   a1 = x * a0 + anf * a1
   b1 = x * b0 + anf * b1
   if (a1 .ne. 0.0) then
      fac = 1.0 / a1
      g = b1 * fac
      if (abs((g - gold) / g) .lt. eps) exit
      gold = g
   endif
end do 
gammcf =  g * exp (-x + a * alog (x) - gln)
return
end

!
!     Evaluate the incomplete gamma function P(a,x), choosing the most 
!     appropriate representation.
!
real function gammap (a, x)
implicit none
real a, x, gamser, gammcf
if (x .lt. a + 1.0) then
   gammap = gamser (a, x)
else
   gammap = 1.0 - gammcf (a, x)
endif
return
end

!
!     Evaluate the incomplete gamma function Q(a,x), choosing the most 
!   appropriate representation.
!
real function gammaq (a, x)
implicit none
real a, x, gamser, gammcf
if (x .lt. a + 1.0) then
   gammaq = 1.0 - gamser (a, x)
else
   gammaq = gammcf (a, x) 
endif
return
end

!
!     For those who don't have a ln(gamma) function.
!
real function gammln(xx)
implicit none
integer j
real xx,x,tmp,ser
real cof(6)
data cof /76.18009173, -86.50532033, 24.01409822, -1.231739516, &
          0.120858003e-2, -0.536382e-5/          
x = xx - 1.0
tmp = x + 5.5
!if(tmp<=0.)then
!print*,'xx,x,tmp',xx,x,tmp
!stop
!endif !xyuan
tmp = tmp - (x+0.5) * alog (tmp)
ser = 1.0
do j = 1, 5
   x = x + 1.0
   ser = ser + cof(j) / x
end do
gammln = -tmp + alog (2.50662827465 * ser)
return
end

function variance(n,tmp)
implicit none
integer :: n,i
real,dimension(n) :: tmp
real :: total, variance
real :: mean
mean=sum(tmp)/float(n)
total=0
do i=1,n
    total=total+(tmp(i)-mean)*(tmp(i)-mean)
end do
variance = total/float(n)
end function
