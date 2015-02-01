module complex_functions

  implicit none

  private

  public :: operator(>)
  public :: operator(<)
  public :: operator(>=)
  public :: operator(<=)
  public :: ccdim, ccmin, ccmax, ccabs, ccsign
  public :: cctanh, ccacos, ccasin, cclog10, ccatan2
  public :: cccosh, ccatan
  public :: ccmaxloc, ccminloc
  public :: ccmaxval, ccminval
  public :: ccrandom_number
  public :: cccpu_time
  public :: ccnint, ccceiling
  public :: extract_imaginary_part
  public :: o

  interface o
    module procedure real_o
    module procedure cmplx_o
  end interface

  interface extract_imaginary_part
    module procedure real_extract_imaginary_part
    module procedure cmplx_extract_imaginary_part
  end interface

  interface ccmaxval
    module procedure ccmaxval_c
    module procedure ccmaxval_c2
    module procedure ccmaxval_r
    module procedure ccmaxval_i
    module procedure ccmaxval_i1
    module procedure ccmaxval_ii
  end interface ccmaxval

  interface ccminval
    module procedure ccminval_c
    module procedure ccminval_c2
    module procedure ccminval_r
    module procedure ccminval_i
    module procedure ccminval_i1
    module procedure ccminval_ii
  end interface ccminval

  interface ccmaxloc
    module procedure ccmaxloc_c
    module procedure ccmaxloc_c2
    module procedure ccmaxloc_c_dim
    module procedure ccmaxloc_i
  end interface ccmaxloc

  interface ccminloc
    module procedure ccminloc_c
    module procedure ccminloc_c2
    module procedure ccminloc_c_dim
    module procedure ccminloc_i
    module procedure ccminloc_i_dim
  end interface ccminloc

  interface cccpu_time
    module procedure cccpu_time_r
    module procedure cccpu_time_c
  end interface cccpu_time

  interface ccceiling
    module procedure ccceiling_r
    module procedure ccceiling_r4
    module procedure ccceiling_c
    module procedure ccceiling_c4
  end interface ccceiling

  interface ccnint
    module procedure ccnint_r
    module procedure ccnint_c
    module procedure ccnint_r4
    module procedure ccnint_c4
  end interface ccnint

  interface ccdim
    module procedure ccdim_cc
  end interface ccdim

  interface ccrandom_number
    module procedure ccrandom_number_r
    module procedure ccrandom_number_r4
    module procedure ccrandom_number_c
    module procedure ccrandom_number_c4
  end interface ccrandom_number

  interface operator(>)
    module procedure gt_ca_c
    module procedure gt_ca_r
    module procedure gt_c_c
    module procedure gt_c4_c4
    module procedure gt_r_c
    module procedure gt_r4_c
    module procedure gt_r4_c4
    module procedure gt_c_r
    module procedure gt_c4_r4
  end interface

  interface operator(<)
    module procedure lt_ca_c
    module procedure lt_c_c
    module procedure lt_c4_c4
    module procedure lt_r_c
    module procedure lt_r4_c
    module procedure lt_r4_c4
    module procedure lt_c_r
  end interface

  interface operator(>=)
    module procedure ge_ca_c
    module procedure ge_c_c
    module procedure ge_c4_c4
    module procedure ge_c_r
    module procedure ge_c_i
    module procedure ge_r_c
    module procedure ge_r4_c
    module procedure ge_i_c
  end interface

  interface operator(<=)
    module procedure le_ca_c
    module procedure le_c_c
    module procedure le_c4_c4
    module procedure le_r_c
    module procedure le_r4_c4
    module procedure le_c_r
  end interface

  interface ccmax
    module procedure ccmax_cc
    module procedure ccmax_c4cd
    module procedure ccmax_c4c4
    module procedure ccmax_r4r4
    module procedure ccmax_r8r8
    module procedure ccmax_ccc
    module procedure ccmax_cccc
    module procedure ccmax_cr
    module procedure ccmax_rc
    module procedure ccmax_ii
    module procedure ccmax_ii1
    module procedure ccmax_iii
    module procedure ccmax_iiii
    module procedure ccmax_iiiii
    module procedure ccmax_iiiiii
    module procedure ccmax_rcc
    module procedure ccmax_i2_i2
    module procedure ccmax_i1_i1
    module procedure ccmax_i1_i
  end interface ccmax

  interface ccmin
    module procedure ccmin_cc
    module procedure ccmin_ccc
    module procedure ccmin_cccc
    module procedure ccmin_cccccc
    module procedure ccmin_cr
    module procedure ccmin_rc
    module procedure ccmin_rr
    module procedure ccmin_ii
    module procedure ccmin_iii
    module procedure ccmin_iiii
    module procedure ccmin_iiiiii
    module procedure ccmin_iiiiiiii
    module procedure ccmin_rcc
    module procedure ccmin_ii_kind_i1
  end interface ccmin

  interface ccacos
    module procedure ccacos_c
    module procedure ccacos_r
  end interface ccacos

  interface ccasin
    module procedure ccasin_c
    module procedure ccasin_r
  end interface ccasin

  interface ccatan2
    module procedure ccatan2_cc
    module procedure ccatan2_rr
  end interface ccatan2

  interface ccatan
    module procedure ccatan_cc
    module procedure ccatan_rr
    module procedure ccatan_cc4
    module procedure ccatan_rr4
  end interface ccatan

  interface ccabs
    module procedure ccabs_c
    module procedure ccabs_c4
    module procedure ccabs_i
    module procedure ccabs_r
    module procedure ccabs_r4
  end interface ccabs

  interface cctanh
    module procedure cctanh_c
    module procedure cctanh_r
  end interface cctanh

  interface cclog10
    module procedure cclog10_c
    module procedure cclog10_c4
    module procedure cclog10_r
    module procedure cclog10_r4
  end interface cclog10

  interface ccsign
    module procedure ccsign_cc
    module procedure ccsign_rc
    module procedure ccsign_ccv
    module procedure ccsign_ii
  end interface ccsign

contains

!================================== CCRANDOM_NUMBER_* ========================80
!
! Complex RANDOM_NUMBER() for various argument types.
!
!=============================================================================80

  subroutine ccrandom_number_r( harvest_r )
    use set_precision, only : r8
    real(r8), intent(out) :: harvest_r
    real(r8) :: harvest
   continue
    call random_number(harvest)
    harvest_r = harvest
  end subroutine ccrandom_number_r

  subroutine ccrandom_number_r4( harvest_r4 )
    use set_precision, only : r4
    real(r4), intent(out) :: harvest_r4
    real(r4) :: harvest
   continue
    call random_number(harvest)
    harvest_r4 = harvest
  end subroutine ccrandom_number_r4

  subroutine ccrandom_number_c( harvest_c )
    use set_precision, only : r8
    complex(r8), intent(out) :: harvest_c
    real(r8) :: harvest
   continue
    call random_number(harvest)
    harvest_c = harvest
  end subroutine ccrandom_number_c

  subroutine ccrandom_number_c4( harvest_c4 )
    use set_precision, only : r4
    complex(r4), intent(out) :: harvest_c4
    real(r4) :: harvest
   continue
    call random_number(harvest)
    harvest_c4 = harvest
  end subroutine ccrandom_number_c4

!===================================== LT_* ==================================80
!
!  Complex OPERATOR(<) for various combinations.
!
!=============================================================================80

  pure function lt_ca_c(a,b)
    use set_precision, only : r8
    complex(r8), dimension(:), intent(in) :: a
    complex(r8),               intent(in) :: b
    logical, dimension(size(a))                  :: lt_ca_c
   continue
    lt_ca_c(:) = (real(a(:),r8) < real(b,r8))
  end function lt_ca_c

  pure function lt_c_c(a,b)
    use set_precision, only : r8
    complex(r8), intent(in) :: a, b
    logical                        :: lt_c_c
   continue
    lt_c_c = ( real(a,r8) < real(b,r8) )
  end function lt_c_c

  pure function lt_c4_c4(c4p1,c4p2)
    use set_precision, only : r4
    complex(r4), intent(in) :: c4p1, c4p2
    logical                        :: lt_c4_c4
   continue
    lt_c4_c4 = ( real(c4p1,r4) < real(c4p2,r4) )
  end function lt_c4_c4

  pure function lt_r_c(a,b)
    use set_precision, only : r8
    real(r8),    intent(in) :: a
    complex(r8), intent(in) :: b
    logical                        :: lt_r_c
   continue
    lt_r_c = ( a < real(b,r8) )
  end function lt_r_c

  pure function lt_r4_c(a,b)
    use set_precision, only : r4, r8
    real(r4),    intent(in) :: a
    complex(r8), intent(in) :: b
    logical                        :: lt_r4_c
   continue
    lt_r4_c = ( a < real(b,r4) )
  end function lt_r4_c

  pure function lt_r4_c4(a,b)
    use set_precision, only : r4
    real(r4),    intent(in) :: a
    complex(r4), intent(in) :: b
    logical                        :: lt_r4_c4
   continue
    lt_r4_c4 = ( a < real(b,r4) )
  end function lt_r4_c4

  pure function lt_c_r(a,b)
    use set_precision, only : r8
    real(r8),    intent(in) :: b
    complex(r8), intent(in) :: a
    logical                        :: lt_c_r
   continue
    lt_c_r = ( real(a,r8) < b )
  end function lt_c_r

!===================================== GT_* ==================================80
!
!  Complex OPERATOR(>) for various combinations.
!
!=============================================================================80

  pure function gt_ca_c(a,b)
    use set_precision, only : r8
    complex(r8), dimension(:), intent(in) :: a
    complex(r8),               intent(in) :: b
    logical, dimension(size(a))                  :: gt_ca_c
   continue
    gt_ca_c(:) = ( real(a(:),r8) > real(b,r8) )
  end function gt_ca_c

  pure function gt_ca_r(a,b)
    use set_precision, only : r8
    complex(r8), dimension(:), intent(in) :: a
    real(r8),                  intent(in) :: b
    logical, dimension(size(a))                  :: gt_ca_r
   continue
    gt_ca_r(:) = ( real(a(:),r8) > b )
  end function gt_ca_r

  pure elemental function gt_c_c(a,b)
    use set_precision, only : r8
    complex(r8), intent(in) :: a, b
    logical                        :: gt_c_c
   continue
    gt_c_c = ( real(a,r8) > real(b,r8) )
  end function gt_c_c

  pure elemental function gt_c4_c4(a,b)
    use set_precision, only : r4
    complex(r4), intent(in) :: a, b
    logical                        :: gt_c4_c4
   continue
    gt_c4_c4 = ( real(a,r4) > real(b,r4) )
  end function gt_c4_c4

  pure elemental function gt_r_c(a,b)
    use set_precision, only : r8
    real(r8),    intent(in) :: a
    complex(r8), intent(in) :: b
    logical                        :: gt_r_c
   continue
    gt_r_c = ( a > real(b,r8) )
  end function gt_r_c

  pure elemental function gt_r4_c(a,b)
    use set_precision, only : r4, r8
    real(r4),    intent(in) :: a
    complex(r8), intent(in) :: b
    logical                        :: gt_r4_c
   continue
    gt_r4_c = ( a > real(b,r4) )
  end function gt_r4_c

  pure elemental function gt_r4_c4(a,b)
    use set_precision, only : r4
    real(r4),    intent(in) :: a
    complex(r4), intent(in) :: b
    logical                        :: gt_r4_c4
   continue
    gt_r4_c4 = ( a > real(b,r4) )
  end function gt_r4_c4

  pure elemental function gt_c_r(a,b)
    use set_precision, only : r8
    real(r8),    intent(in) :: b
    complex(r8), intent(in) :: a
    logical                        :: gt_c_r
   continue
    gt_c_r = ( real(a,r8) > b )
  end function gt_c_r

  pure elemental function gt_c4_r4(c4p,r4p)
    use set_precision, only : r4
    complex(r4), intent(in) :: c4p
    real(r4),    intent(in) :: r4p
    logical                        :: gt_c4_r4
   continue
    gt_c4_r4 = ( real(c4p,r4) > r4p )
  end function gt_c4_r4

!===================================== LE_* ==================================80
!
!  Complex OPERATOR(<=) for various combinations.
!
!=============================================================================80

  pure function le_ca_c(a,b)
    use set_precision, only : r8
    complex(r8), dimension(:), intent(in) :: a
    complex(r8),               intent(in) :: b
    logical, dimension(size(a))                  :: le_ca_c
   continue
    le_ca_c(:) = ( real(a(:),r8) <= real(b,r8) )
  end function le_ca_c

  pure elemental function le_c_c(a,b)
    use set_precision, only : r8
    complex(r8), intent(in) :: a, b
    logical                        :: le_c_c
   continue
    le_c_c = ( real(a,r8) <= real(b,r8) )
  end function le_c_c

  pure elemental function le_c4_c4(c4p1,c4p2)
    use set_precision, only : r4
    complex(r4), intent(in) :: c4p1, c4p2
    logical                        :: le_c4_c4
   continue
    le_c4_c4 = ( real(c4p1,r4) <= real(c4p2,r4) )
  end function le_c4_c4

  pure elemental function le_r_c(a,b)
    use set_precision, only : r8
    real(r8),    intent(in) :: a
    complex(r8), intent(in) :: b
    logical                        :: le_r_c
   continue
    le_r_c = ( a <= real(b,r8) )
  end function le_r_c

  pure elemental function le_r4_c4(a,b)
    use set_precision, only : r4
    real(r4),    intent(in) :: a
    complex(r4), intent(in) :: b
    logical                        :: le_r4_c4
   continue
    le_r4_c4 = ( a <= real(b,r4) )
  end function le_r4_c4

  pure elemental function le_c_r(a,b)
    use set_precision, only : r8
    real(r8),    intent(in) :: b
    complex(r8), intent(in) :: a
    logical                        :: le_c_r
   continue
    le_c_r = ( real(a,r8) <= b )
  end function le_c_r

!===================================== GE_* ==================================80
!
!  Complex OPERATOR(>=) for various combinations.
!
!=============================================================================80

  pure function ge_ca_c(a,b)
    use set_precision, only : r8
    complex(r8), dimension(:), intent(in) :: a
    complex(r8),               intent(in) :: b
    logical, dimension(size(a))                  :: ge_ca_c
   continue
    ge_ca_c(:) = ( real(a(:),r8) >= real(b,r8) )
  end function ge_ca_c

  pure elemental function ge_c_c(a,b)
    use set_precision, only : r8
    complex(r8), intent(in) :: a, b
    logical                        :: ge_c_c
   continue
    ge_c_c = ( real(a,r8) >= real(b,r8) )
  end function ge_c_c

  pure elemental function ge_c4_c4(c4p1,c4p2)
    use set_precision, only : r4
    complex(r4), intent(in) :: c4p1, c4p2
    logical                        :: ge_c4_c4
   continue
    ge_c4_c4 = ( real(c4p1,r4) >= real(c4p2,r4) )
  end function ge_c4_c4

  pure elemental function ge_c_r(a,b)
    use set_precision, only : r8
    complex(r8), intent(in) :: a
    real(r8),    intent(in) :: b
    logical                        :: ge_c_r
   continue
    ge_c_r = ( real(a,r8) >= b )
  end function ge_c_r

  pure elemental function ge_c_i(a,b)
    use set_precision, only : r8
    complex(r8), intent(in) :: a
    integer,            intent(in) :: b
    logical                        :: ge_c_i
   continue
    ge_c_i = ( real(a,r8) >= b )
  end function ge_c_i

  pure elemental function ge_r_c(a,b)
    use set_precision, only : r8
    real(r8),    intent(in) :: a
    complex(r8), intent(in) :: b
    logical                        :: ge_r_c
   continue
    ge_r_c = ( a >= real(b,r8) )
  end function ge_r_c

  pure elemental function ge_r4_c(a,b)
    use set_precision, only : r8, r4
    real(r4),    intent(in) :: a
    complex(r8), intent(in) :: b
    logical                        :: ge_r4_c
   continue
    ge_r4_c = ( a >= real(b,r8) )
  end function ge_r4_c

  pure elemental function ge_i_c(a,b)
    use set_precision, only : r8
    integer,            intent(in) :: a
    complex(r8), intent(in) :: b
    logical                        :: ge_i_c
   continue
    ge_i_c = ( a >= real(b,r8) )
  end function ge_i_c

!================================== CCDIM_* ==================================80
!
! Complex DIM() for various argument types.
!
!=============================================================================80

  elemental function ccdim_cc(a1,a2)
    use set_precision, only : r8
    complex(r8), intent(in) :: a1, a2
    complex(r8)             :: ccdim_cc
   continue
    ccdim_cc = ccmax_cr(a2-a1, 0.0_r8)
  end function ccdim_cc

!================================== CCMAX_* ==================================80
!
! Complex MAX() for various argument types.
!
!=============================================================================80

  elemental function ccmax_cc(a1,a2)
    use set_precision, only : r8
    complex(r8), intent(in) :: a1, a2
    complex(r8)             :: ccmax_cc
   continue
    if ( real(a1,r8) > real(a2,r8) ) then
      ccmax_cc = a1
    else
      ccmax_cc = a2
    endif
  end function ccmax_cc

  elemental function ccmax_c4cd(c4p,cdp)
    use set_precision, only : r4, r8
    complex(r4), intent(in) :: c4p
    complex(r8), intent(in) :: cdp
    complex(r4)             :: ccmax_c4cd
   continue
    if ( real(c4p,r4) > real(cdp,r4) ) then
      ccmax_c4cd = c4p
    else
      ccmax_c4cd = cmplx(real(cdp,r4),aimag(cdp),r4)
    endif
  end function ccmax_c4cd

  elemental function ccmax_c4c4(c4p1,c4p2)
    use set_precision, only : r4
    complex(r4), intent(in) :: c4p1, c4p2
    complex(r4)             :: ccmax_c4c4
   continue
    if ( real(c4p1,r4) > real(c4p2,r4) ) then
      ccmax_c4c4 = c4p1
    else
      ccmax_c4c4 = c4p2
    endif
  end function ccmax_c4c4

  elemental function ccmax_r4r4(r4p1,r4p2)
    use set_precision, only : r4
    real(r4), intent(in) :: r4p1, r4p2
    real(r4)             :: ccmax_r4r4
   continue
    ccmax_r4r4 = max(r4p1, r4p2)
  end function ccmax_r4r4

  elemental function ccmax_r8r8(r8p1,r8p2)
    use set_precision, only : r8
    real(r8), intent(in) :: r8p1, r8p2
    real(r8)             :: ccmax_r8r8
   continue
    ccmax_r8r8 = max(r8p1, r8p2)
  end function ccmax_r8r8

  elemental function ccmax_ccc(a1,a2,a3)
    use set_precision, only : r8
    complex(r8), intent(in) :: a1, a2, a3
    complex(r8)             :: ccmax_ccc
   continue
    if ( real(a1,r8) > real(a2,r8) ) then
      if ( real(a1,r8) > real(a3,r8) ) then
        ccmax_ccc = a1
      else
        ccmax_ccc = a3
      endif
    else
      if ( real(a2,r8) > real(a3,r8) ) then
        ccmax_ccc = a2
      else
        ccmax_ccc = a3
      endif
    endif
  end function ccmax_ccc

  elemental function ccmax_cccc(a1,a2,a3,a4)
    use set_precision, only : r8
    complex(r8), intent(in) :: a1, a2, a3, a4
    complex(r8)             :: ccmax_cccc
   continue
    ccmax_cccc = ccmax( ccmax( a1, a2 ), ccmax( a3, a4 ) )
  end function ccmax_cccc

  elemental function ccmax_cr(a1,a2)
    use set_precision, only : r8
    complex(r8), intent(in) :: a1
    real(r8),    intent(in) :: a2
    complex(r8)             :: ccmax_cr
   continue
    if ( real(a1,r8) > a2 ) then
      ccmax_cr = a1
    else
      ccmax_cr = cmplx(a2, 0.0_r8, r8)
    endif
  end function ccmax_cr

  elemental function ccmax_rc(a1,a2)
    use set_precision, only : r8
    real(r8),    intent(in) :: a1
    complex(r8), intent(in) :: a2
    complex(r8)             :: ccmax_rc
   continue
    if ( a1 > real(a2,r8) ) then
      ccmax_rc = cmplx(a1, 0.0_r8, r8)
    else
      ccmax_rc = a2
    endif
  end function ccmax_rc

  elemental function ccmax_rcc(a1,a2,a3)
    use set_precision, only : r8
    real(r8),    intent(in) :: a1
    complex(r8), intent(in) :: a2
    complex(r8), intent(in) :: a3
    complex(r8)             :: ccmax_rcc
   continue
    ccmax_rcc = ccmax( ccmax( a1, a2), a3)
  end function ccmax_rcc

  elemental function ccmax_ii(a1,a2)
    integer, intent(in) :: a1, a2
    integer             :: ccmax_ii
   continue
    if ( a1 > a2 ) then
      ccmax_ii = a1
    else
      ccmax_ii = a2
    endif
  end function ccmax_ii

  elemental function ccmax_ii1(a1,a2)
    use set_precision, only : system_i1
    integer,            intent(in) :: a1
    integer(system_i1), intent(in) :: a2
    integer                        :: ccmax_ii1
   continue
    if ( a1 > a2 ) then
      ccmax_ii1 = a1
    else
      ccmax_ii1 = a2
    endif
  end function ccmax_ii1

  elemental function ccmax_iii(a1,a2,a3)
    integer, intent(in) :: a1, a2, a3
    integer             :: ccmax_iii
   continue
    ccmax_iii = ccmax(ccmax(a1,a2),a3)
  end function ccmax_iii

  elemental function ccmax_iiii(a1,a2,a3,a4)
    integer, intent(in) :: a1, a2, a3, a4
    integer             :: ccmax_iiii
   continue
    ccmax_iiii = ccmax(ccmax(a1,a2),ccmax(a3,a4))
  end function ccmax_iiii

  elemental function ccmax_iiiii(a1,a2,a3,a4,a5)
    integer, intent(in) :: a1, a2, a3, a4, a5
    integer             :: ccmax_iiiii
   continue
    ccmax_iiiii = ccmax(ccmax(a1,a2),ccmax(a3,a4,a5))
  end function ccmax_iiiii

  elemental function ccmax_iiiiii(a1,a2,a3,a4,a5,a6)
    integer, intent(in) :: a1, a2, a3, a4, a5, a6
    integer             :: ccmax_iiiiii
   continue
    ccmax_iiiiii = ccmax(ccmax(a1,a2,a3),ccmax(a4,a5,a6))
  end function ccmax_iiiiii

!================================== CCMIN_* ==================================80
!
! Complex MIN() for various argument types.
!
!=============================================================================80

  elemental function ccmin_cc(a1,a2)
    use set_precision, only : r8
    complex(r8), intent(in)  :: a1, a2
    complex(r8)              :: ccmin_cc
   continue
    if ( real(a1,r8) < real(a2,r8) ) then
      ccmin_cc = a1
    else
      ccmin_cc = a2
    endif
  end function ccmin_cc

  elemental function ccmin_ccc(a1,a2,a3)
    use set_precision, only : r8
    complex(r8), intent(in)  :: a1, a2, a3
    complex(r8)              :: ccmin_ccc
   continue
    ccmin_ccc = ccmin(a1,ccmin(a2,a3))
  end function ccmin_ccc

  elemental function ccmin_cccc(a1,a2,a3,a4)
    use set_precision, only : r8
    complex(r8), intent(in)  :: a1, a2, a3, a4
    complex(r8)              :: ccmin_cccc
   continue
    ccmin_cccc = ccmin(ccmin(a1,a2),ccmin(a3,a4))
  end function ccmin_cccc

  elemental function ccmin_cccccc(a1,a2,a3,a4,a5,a6)
    use set_precision, only : r8
    complex(r8), intent(in)  :: a1, a2, a3, a4, a5, a6
    complex(r8)              :: ccmin_cccccc
   continue
    ccmin_cccccc = ccmin(ccmin(a1,a2),ccmin(a3,a4),ccmin(a5,a6))
  end function ccmin_cccccc

  elemental function ccmin_cr(a1,a2)
    use set_precision, only : r8
    complex(r8), intent(in) :: a1
    real(r8),    intent(in) :: a2
    complex(r8)             :: ccmin_cr
   continue
    if ( real(a1,r8) < a2 ) then
      ccmin_cr = a1
    else
      ccmin_cr = cmplx(a2, 0.0_r8, r8)
    endif
  end function ccmin_cr

  elemental function ccmin_rc(a1,a2)
    use set_precision, only : r8
    real(r8),    intent(in) :: a1
    complex(r8), intent(in) :: a2
    complex(r8)             :: ccmin_rc
   continue
    if ( a1 < real(a2,r8) ) then
      ccmin_rc = cmplx(a1, 0.0_r8, r8)
    else
      ccmin_rc = a2
    endif
  end function ccmin_rc

  elemental function ccmin_rr(a1,a2)
    use set_precision, only : r8
    real(r8), intent(in) :: a1
    real(r8), intent(in) :: a2
    complex(r8)          :: ccmin_rr
   continue
    ccmin_rr = min(a1,a2)
  end function ccmin_rr

  elemental function ccmin_rcc(a1,a2,a3)
    use set_precision, only : r8
    real(r8),    intent(in) :: a1
    complex(r8), intent(in) :: a2
    complex(r8), intent(in) :: a3
    complex(r8)             :: ccmin_rcc
   continue
    ccmin_rcc = ccmin( ccmin( a1, a2 ), a3 )
  end function ccmin_rcc

  elemental function ccmax_i2_i2(a1,a2)
    use set_precision, only : system_i2
    integer(system_i2), intent(in) :: a1, a2
    integer(system_i2)             :: ccmax_i2_i2
   continue
    if ( a1 > a2 ) then
      ccmax_i2_i2 = a1
    else
      ccmax_i2_i2 = a2
    endif
  end function ccmax_i2_i2

  elemental function ccmax_i1_i1(a1,a2)
    use set_precision, only : system_i1
    integer(system_i1), intent(in) :: a1, a2
    integer(system_i1)             :: ccmax_i1_i1
   continue
    if ( a1 > a2 ) then
      ccmax_i1_i1 = a1
    else
      ccmax_i1_i1 = a2
    endif
  end function ccmax_i1_i1

  elemental function ccmax_i1_i(a1,a2)
    use set_precision, only : system_i1
    integer(system_i1), intent(in) :: a1
    integer,            intent(in) :: a2
    integer                        :: ccmax_i1_i
   continue
    if ( a1 > a2 ) then
      ccmax_i1_i = a1
    else
      ccmax_i1_i = a2
    endif
  end function ccmax_i1_i

  elemental function ccmin_ii(a1,a2)
    integer, intent(in) :: a1, a2
    integer             :: ccmin_ii
   continue
    ccmin_ii = min(a1,a2)
  end function ccmin_ii

  elemental function ccmin_iii(a1,a2,a3)
    integer, intent(in) :: a1, a2, a3
    integer             :: ccmin_iii
   continue
    ccmin_iii = min(a1,a2,a3)
  end function ccmin_iii

  elemental function ccmin_iiii(a1,a2,a3,a4)
    integer, intent(in) :: a1, a2, a3, a4
    integer             :: ccmin_iiii
   continue
    ccmin_iiii = min(a1,a2,a3,a4)
  end function ccmin_iiii

  elemental function ccmin_iiiiii(a1,a2,a3,a4,a5,a6)
    integer, intent(in) :: a1, a2, a3, a4, a5, a6
    integer             :: ccmin_iiiiii
   continue
    ccmin_iiiiii = min(a1,a2,a3,a4,a5,a6)
 end function ccmin_iiiiii

  elemental function ccmin_iiiiiiii(a1,a2,a3,a4,a5,a6,a7,a8)
    integer, intent(in) :: a1, a2, a3, a4, a5, a6, a7, a8
    integer             :: ccmin_iiiiiiii
   continue
    ccmin_iiiiiiii = min(a1,a2,a3,a4,a5,a6,a7,a8)
 end function ccmin_iiiiiiii

  elemental function ccmin_ii_kind_i1(a1,a2)
    use set_precision, only : system_i1
    integer(system_i1), intent(in) :: a1, a2
    integer(system_i1)             :: ccmin_ii_kind_i1
   continue
    ccmin_ii_kind_i1 = min(a1,a2)
  end function ccmin_ii_kind_i1

!================================ CCCPU_TIME_* ===============================80
!
! Complex CPU_TIME() for various argument types.
!
!=============================================================================80

  subroutine cccpu_time_r( time_r )
    real, intent(out) :: time_r !FIXME: needs r8
   continue
    call cpu_time( time_r )
  end subroutine cccpu_time_r

  subroutine cccpu_time_c( time_c )
    complex, intent(out) :: time_c !FIXME: needs r8
    real :: time_r !FIXME: needs r8
   continue
    call cpu_time( time_r )
    time_c = time_r !FIXME: needs r8
  end subroutine cccpu_time_c

!================================ CCCEILING_* ================================80
!
! Complex CEILING() for various argument types.
!
!=============================================================================80

  elemental function ccceiling_r(a)
    use set_precision, only : r8
    real(r8), intent(in) :: a
    integer                     :: ccceiling_r
   continue
    ccceiling_r = ceiling( a )
  end function ccceiling_r

  elemental function ccceiling_r4(a)
    use set_precision, only : r4
    real(r4), intent(in) :: a
    integer                     :: ccceiling_r4
   continue
    ccceiling_r4 = ceiling( a )
  end function ccceiling_r4

  elemental function ccceiling_c(a)
    use set_precision, only : r8
    complex(r8), intent(in) :: a
    integer                        :: ccceiling_c
   continue
    ccceiling_c = ceiling( real(a,r8) )
  end function ccceiling_c

  elemental function ccceiling_c4(a)
    use set_precision, only : r4
    complex(r4), intent(in) :: a
    integer                        :: ccceiling_c4
   continue
    ccceiling_c4 = ceiling( real(a,r4) )
  end function ccceiling_c4

!================================= CCNINT_* ==================================80
!
! Complex NINT() for various argument types.
!
!=============================================================================80

  elemental function ccnint_r(a)
    use set_precision, only : r8
    real(r8), intent(in) :: a
    integer                     :: ccnint_r
   continue
    ccnint_r = nint( a )
  end function ccnint_r

  elemental function ccnint_c(a)
    use set_precision, only : r8
    complex(r8), intent(in) :: a
    integer                        :: ccnint_c
   continue
    ccnint_c = nint( real(a,r8) )
  end function ccnint_c

  elemental function ccnint_r4(a)
    use set_precision, only : r4
    real(r4), intent(in) :: a
    integer                     :: ccnint_r4
   continue
    ccnint_r4 = nint( a )
  end function ccnint_r4

  elemental function ccnint_c4(a)
    use set_precision, only : r4
    complex(r4), intent(in) :: a
    integer                        :: ccnint_c4
   continue
    ccnint_c4 = nint( real(a,r4) )
  end function ccnint_c4

!================================== CCABS_* ==================================80
!
! Complex ABS() for various argument types.
!
!=============================================================================80

  pure elemental function ccabs_c(a)
    use set_precision, only : r8
    complex(r8), intent(in) :: a
    complex(r8)             :: ccabs_c
   continue
    if (real(a,r8) < 0.0_r8) then
      ccabs_c = -a
    else
      ccabs_c = a
    end if
  end function ccabs_c

  pure elemental function ccabs_c4(c4p)
    use set_precision, only : r4
    complex(r4), intent(in) :: c4p
    complex(r4)             :: ccabs_c4
   continue
    if ( real(c4p,r4) > 0.0_r4 ) then
      ccabs_c4 =  c4p
    else
      ccabs_c4 = -c4p
    endif
  end function ccabs_c4

  pure elemental function ccabs_r(a)
    use set_precision, only : r8
    real(r8), intent(in) :: a
    real(r8)             :: ccabs_r
   continue
    ccabs_r = abs( a )
  end function ccabs_r

  pure elemental function ccabs_r4(a)
    use set_precision, only : r4
    real(r4), intent(in) :: a
    real(r4)             :: ccabs_r4
   continue
    ccabs_r4 = abs( a )
  end function ccabs_r4

  pure elemental function ccabs_i(a)
    integer, intent(in) :: a
    integer             :: ccabs_i
   continue
    ccabs_i = abs( a )
  end function ccabs_i

!================================== CCACOS_* =================================80
!
! Complex ACOS() for various argument types.
!
!=============================================================================80

  pure function ccacos_r(a)
    use set_precision, only : r8
    real(r8), intent(in) :: a
    real(r8)             :: ccacos_r
   continue
    ccacos_r = acos(a)
  end function ccacos_r

  pure function ccacos_c(a)
    use set_precision, only : r8
    complex(r8), intent(in) :: a
    complex(r8)             :: ccacos_c
    real(r8) :: x2, y2, xp1, xm1, part1, part2, apiece, bpiece
   continue
!   from http://astronomy.swin.edu.au/~pbourke/oldstuff/mathlib/ccl.c
    if ( aimag(a) < epsilon(1.0_r8) ) then
      ccacos_c = cmplx(acos(real(a,r8)), 0.0_r8, r8)
    else
      x2 = real(a,r8) * real(a,r8)
      y2 = aimag(a) * aimag(a)
      xp1 = x2 + 2.0_r8 * real(a,r8) + 1.0_r8
      xm1 = x2 - 2.0_r8 * real(a,r8) + 1.0_r8
      part1 = 0.5_r8 * sqrt(xp1 + y2)
      part2 = 0.5_r8 * sqrt(xm1 + y2)
      apiece = part1 + part2
      bpiece = part1 - part2
      ccacos_c = cmplx(acos(bpiece),                                           &
                       log(apiece+sqrt(apiece*apiece-1.0_r8)),          &
                       r8)
    endif
  end function ccacos_c

!================================== CCASIN_* =================================80
!
! Complex ASIN() for various argument types.
!
!=============================================================================80

  function ccasin_r(a)
    use set_precision, only : r8
    real(r8), intent(in) :: a
    real(r8)             :: ccasin_r
   continue
    ccasin_r = asin(a)
  end function ccasin_r

! FIXME: not completely general - see ccacos_c
  function ccasin_c(a)
    use set_precision, only : r8
    complex(r8), intent(in) :: a
    complex(r8)             :: ccasin_c
   continue
    ccasin_c = cmplx( asin(real(a, r8)),                                &
                      aimag(a)/sqrt(1.0_r8 - real(a, r8)**2),    &
                      r8 )
  end function ccasin_c

!================================== CCATAN2_* ================================80
!
! Complex ATAN2() for various argument types.
!
!=============================================================================80

  pure function ccatan2_rr(a,b)
    use set_precision, only : r8
    real(r8), intent(in) :: a, b
    real(r8)             :: ccatan2_rr
   continue
    ccatan2_rr = atan2(a,b)
  end function ccatan2_rr

! FIXME: not completely general - see ccacos_c
  pure function ccatan2_cc(a,b)
    use set_precision, only : r8
    complex(r8), intent(in) :: a,b
    complex(r8) :: ccatan2_cc
   continue
    ccatan2_cc = cmplx( atan2(real(a,r8),real(b,r8)),            &
                        (real(b,r8)*aimag(a)-real(a,r8)*aimag(b))&
                       /(real(a, r8)**2+real(b,r8)**2), r8)
  end function ccatan2_cc

!================================== CCATAN_* =================================80
!
! Complex ATAN() for various argument types.
!
!=============================================================================80

  function ccatan_rr(a)
    use set_precision, only : r8
    real(r8), intent(in) :: a
    real(r8)             :: ccatan_rr
   continue
    ccatan_rr = atan(a)
  end function ccatan_rr

  function ccatan_rr4(a)
    use set_precision, only : r4
    real(r4), intent(in) :: a
    real(r4)             :: ccatan_rr4
   continue
    ccatan_rr4 = atan(a)
  end function ccatan_rr4

! FIXME: not completely general - see ccacos_c
  function ccatan_cc(a)
    use set_precision, only : r8
    complex(r8), intent(in) :: a
    complex(r8)             :: ccatan_cc
   continue
    ccatan_cc = cmplx(atan(real(a)),aimag(a)/(1.+real(a)**2),r8)
  end function ccatan_cc

! FIXME: not completely general - see ccacos_c
  function ccatan_cc4(a)
    use set_precision, only : r4
    complex(r4), intent(in) :: a
    complex(r4)             :: ccatan_cc4
   continue
    ccatan_cc4 = cmplx(atan(real(a)),aimag(a)/(1.+real(a)**2),r4)
  end function ccatan_cc4

!================================== CCSIGN_* =================================80
!
! Complex SIGN() for various argument types.
!
!=============================================================================80

  pure function ccsign_ccv(a,b)
    use set_precision, only : r8
    complex(r8),               intent(in) :: a
    complex(r8), dimension(:), intent(in) :: b
    complex(r8), dimension(size(b,1))     :: ccsign_ccv
    real(r8), dimension(size(b,1)) :: sgnb
    integer :: i
   continue
    sgnb(:) = 1.0_r8
    do i=1,size(b,1)
      if (real(b(i),r8) < 0.0_r8) sgnb(i) = -1.0_r8
    enddo
    ccsign_ccv = ccabs(a)*sgnb
  end function ccsign_ccv

  pure function ccsign_cc(a,b)
    use set_precision, only : r8
    complex(r8), intent(in) :: a, b
    complex(r8)             :: ccsign_cc
    real(r8) :: sgnb
   continue
    sgnb = 1.0_r8
    if (real(b,r8) < 0.0_r8) sgnb = -1.0_r8
    ccsign_cc = ccabs(a)*sgnb
  end function ccsign_cc

  pure function ccsign_ii(a,b)
    integer, intent(in) :: a, b
    integer             :: ccsign_ii
    integer :: sgnb
   continue
    sgnb = 1
    if ( b < 0 ) sgnb = -1
    ccsign_ii = abs(a)*sgnb
  end function ccsign_ii

  pure function ccsign_rc(a,b)
    use set_precision, only : r8
    real(r8),    intent(in) :: a
    complex(r8), intent(in) :: b
    real(r8)                :: ccsign_rc
    real(r8) :: sgnb
   continue
    sgnb = 1.0_r8
    if (real(b,r8) < 0.0_r8) sgnb = -1.0_r8
    ccsign_rc = abs(a)*sgnb
  end function ccsign_rc

!================================== CCTANH_* =================================80
!
! Complex TANH() for various argument types.
!
!=============================================================================80

  pure function cctanh_r(a)
    use set_precision, only : r8
    real(r8), intent(in) :: a
    complex(r8)          :: cctanh_r
   continue
    cctanh_r = cmplx( tanh(a), 0.0_r8, r8 )
  end function cctanh_r

  pure function cctanh_c(a)
    use set_precision, only : r8
    complex(r8), intent(in) :: a
    complex(r8)             :: cctanh_c
    complex(r8) :: eplus, eminus
   continue
    if (real(a,r8) > 50.0_r8) then ! magic number
      cctanh_c = 1.0_r8
    else
      eplus  = exp(a)
      eminus = exp(-a)
      cctanh_c = (eplus - eminus)/(eplus + eminus)
    end if
  end function cctanh_c

!================================== CCLOG10_* ================================80
!
! Complex LOG10() from various argument types.
!
!=============================================================================80

  elemental function cclog10_c(a)
    use set_precision, only : r8
    complex(r8), intent(in) :: a
    complex(r8)             :: cclog10_c
   continue
    cclog10_c = log(a) / log((10.0_r8, 0.0_r8))
  end function cclog10_c

  elemental function cclog10_c4(a)
    use set_precision, only : r4
    complex(r4), intent(in) :: a
    complex(r4)             :: cclog10_c4
   continue
    cclog10_c4 = log(a) / log((10.0_r4, 0.0_r4))
  end function cclog10_c4

  elemental function cclog10_r(a)
    use set_precision, only : r8
    real(r8), intent(in) :: a
    real(r8)             :: cclog10_r
   continue
    cclog10_r = log10(a)
  end function cclog10_r

  elemental function cclog10_r4(a)
    use set_precision, only : r4
    real(r4), intent(in) :: a
    real(r4)             :: cclog10_r4
   continue
    cclog10_r4 = log10(a)
  end function cclog10_r4

!================================== CCCOSH_* =================================80
!
! Complex COSH() for various argument types.
!
!=============================================================================80

  function cccosh(a)
    use set_precision, only : r8
    complex(r8), intent(in) :: a
    complex(r8)             :: cccosh
   continue
    cccosh = cmplx( cosh(real(a, r8)),                                  &
                    aimag(a)*sinh(real(a, r8)),                         &
                    r8 )
  end function cccosh

!================================== CCMAXVAL_* ===============================80
!
! Complex MAXVAL() for various argument types.
!
!=============================================================================80

  function ccmaxval_c(a)
    use set_precision, only : r8
    complex(r8), dimension(:), intent(in) :: a
    complex(r8)                           :: ccmaxval_c
    integer, dimension(1) :: max_location
   continue
    max_location = ccmaxloc(a)
    ccmaxval_c = a(max_location(1))
  end function ccmaxval_c

  function ccmaxval_c2(a)
    use set_precision, only : r8
    complex(r8), dimension(:,:), intent(in) :: a
    complex(r8)                             :: ccmaxval_c2
    integer, dimension(2) :: max_location
   continue
    max_location = ccmaxloc(a)
    ccmaxval_c2 = a(max_location(1),max_location(2))
  end function ccmaxval_c2

  function ccmaxval_r(a)
    use set_precision, only : r8
    real(r8), dimension(:), intent(in) :: a
    real(r8)                           :: ccmaxval_r
   continue
    ccmaxval_r = maxval(a)
  end function ccmaxval_r

  function ccmaxval_i(a)
    integer, dimension(:), intent(in) :: a
    integer                           :: ccmaxval_i
   continue
    ccmaxval_i = maxval(a)
  end function ccmaxval_i

  function ccmaxval_i1(a)
    use set_precision, only : system_i1
    integer(system_i1), dimension(:), intent(in) :: a
    integer(system_i1)                           :: ccmaxval_i1
   continue
    ccmaxval_i1 = maxval(a)
  end function ccmaxval_i1

  function ccmaxval_ii(a)
    integer, dimension(:,:), intent(in) :: a
    integer                             :: ccmaxval_ii
   continue
    ccmaxval_ii = maxval(a)
  end function ccmaxval_ii

!================================== CCMINVAL_* ===============================80
!
! Complex MINVAL() for various argument types.
!
!=============================================================================80

  pure function ccminval_c(a)
    use set_precision, only : r8
    complex(r8), dimension(:), intent(in) :: a
    complex(r8)                           :: ccminval_c
    integer, dimension(1) :: min_location
   continue
    min_location = ccminloc(a)
    ccminval_c = a(min_location(1))
  end function ccminval_c

  pure function ccminval_c2(a)
    use set_precision, only : r8
    complex(r8), dimension(:,:), intent(in) :: a
    complex(r8)                             :: ccminval_c2
    integer, dimension(2) :: min_location
   continue
    min_location = ccminloc(a)
    ccminval_c2 = a(min_location(1),min_location(2))
  end function ccminval_c2

  pure function ccminval_r(a)
    use set_precision, only : r8
    real(r8), dimension(:), intent(in) :: a
    real(r8)                           :: ccminval_r
   continue
    ccminval_r = minval(a)
  end function ccminval_r

  pure function ccminval_i(a)
    integer, dimension(:), intent(in) :: a
    integer                           :: ccminval_i
   continue
    ccminval_i = minval(a)
  end function ccminval_i

  pure function ccminval_i1(a)
    use set_precision, only : system_i1
    integer(system_i1), dimension(:), intent(in) :: a
    integer(system_i1)                           :: ccminval_i1
   continue
    ccminval_i1 = minval(a)
  end function ccminval_i1

  pure function ccminval_ii(a)
    integer, dimension(:,:), intent(in) :: a
    integer                             :: ccminval_ii
   continue
    ccminval_ii = minval(a)
  end function ccminval_ii

!================================== CCMAXLOC_* ===============================80
!
! Complex MAXLOC() for various argument types.
!
!=============================================================================80

  function ccmaxloc_c(a, mask)
    use set_precision, only : r8
    complex(r8), dimension(:),      intent(in) :: a
    logical, dimension(size(a)), optional, intent(in) :: mask
    integer                                           :: ccmaxloc_c(1)
   continue
    mask_conditional : if ( present(mask) ) then
      ccmaxloc_c = maxloc( real( a, r8 ) , mask )
    else mask_conditional
      ccmaxloc_c = maxloc( real( a, r8 ) )
    end if mask_conditional
  end function ccmaxloc_c

  function ccmaxloc_c2(a, mask)
    use set_precision, only : r8
    complex(r8), dimension(:,:),                intent(in) :: a
    logical, dimension(size(a,1),size(a,2)), optional, intent(in) :: mask
    integer                                                    :: ccmaxloc_c2(2)
   continue
    mask_conditional : if ( present(mask) ) then
      ccmaxloc_c2 = maxloc( real( a, r8 ) , mask )
    else mask_conditional
      ccmaxloc_c2 = maxloc( real( a, r8 ) )
    end if mask_conditional
  end function ccmaxloc_c2

  function ccmaxloc_c_dim(a, dim)
    use set_precision, only : r8
    complex(r8), dimension(:), intent(in) :: a
    integer,                          intent(in) :: dim
    integer                                      :: ccmaxloc_c_dim
   continue
    ccmaxloc_c_dim = maxloc( real( a, r8 ), dim )
  end function ccmaxloc_c_dim

  function ccmaxloc_i(a, mask)
    integer, dimension(:),                 intent(in) :: a
    logical, dimension(size(a)), optional, intent(in) :: mask
    integer                                           :: ccmaxloc_i(1)
   continue
    mask_conditional : if ( present(mask) ) then
      ccmaxloc_i = maxloc( a, mask )
    else mask_conditional
      ccmaxloc_i = maxloc( a )
    end if mask_conditional
  end function ccmaxloc_i

!================================== CCMINLOC_* ===============================80
!
! Complex MINLOC() for various argument types.
!
!=============================================================================80

  pure function ccminloc_c(a, mask)
    use set_precision, only : r8
    complex(r8), dimension(:),      intent(in) :: a
    logical, dimension(size(a)), optional, intent(in) :: mask
    integer                                           :: ccminloc_c(1)
   continue
    mask_conditional : if ( present(mask) ) then
      ccminloc_c = minloc( real( a, r8 ) , mask )
    else mask_conditional
      ccminloc_c = minloc( real( a, r8 ) )
    end if mask_conditional
  end function ccminloc_c

  pure function ccminloc_c2(a, mask)
    use set_precision, only : r8
    complex(r8), dimension(:,:),                intent(in) :: a
    logical, dimension(size(a,1),size(a,2)), optional, intent(in) :: mask
    integer                                                    :: ccminloc_c2(2)
   continue
    mask_conditional : if ( present(mask) ) then
      ccminloc_c2 = minloc( real( a, r8 ) , mask )
    else mask_conditional
      ccminloc_c2 = minloc( real( a, r8 ) )
    end if mask_conditional
  end function ccminloc_c2

  pure function ccminloc_c_dim(a, dim)
    use set_precision, only : r8
    complex(r8), dimension(:), intent(in) :: a
    integer,                          intent(in) :: dim
    integer :: ccminloc_c_dim
   continue
    ccminloc_c_dim = minloc( real( a, r8 ), dim )
  end function ccminloc_c_dim

  pure function ccminloc_i_dim(a, dim)
    integer, dimension(:), intent(in) :: a
    integer,               intent(in) :: dim
    integer                           :: ccminloc_i_dim
   continue
    ccminloc_i_dim = minloc( a, dim )
  end function ccminloc_i_dim

  pure function ccminloc_i(a, mask)
    integer, dimension(:),                 intent(in) :: a
    logical, dimension(size(a)), optional, intent(in) :: mask
    integer                                           :: ccminloc_i(1)
   continue
    mask_conditional : if ( present(mask) ) then
      ccminloc_i = minloc( a, mask )
    else mask_conditional
      ccminloc_i = minloc( a )
    end if mask_conditional
  end function ccminloc_i

!========================== REAL_EXTRACT_IMAGINARY_PART  =====================80
!
! Allow compilation of complex code before complexification.
!
!=============================================================================80

  subroutine real_extract_imaginary_part(arg_in,arg_out)

    use set_precision, only : dp

    real(dp), intent(in)  :: arg_in
    real(dp), intent(out) :: arg_out

    continue

    arg_out = arg_in

  end subroutine real_extract_imaginary_part

!================================= CMPLX_EXTRACT_IMAGINARY_PART  =============80
!
! Imaginary part of arg_in returned in real part of arg_out.
!
!=============================================================================80

  subroutine cmplx_extract_imaginary_part(arg_in,arg_out)

    use set_precision, only : dp

    complex(dp), intent(in)  :: arg_in
    complex(dp), intent(out) :: arg_out

    continue

    arg_out = aimag(arg_in)

  end subroutine cmplx_extract_imaginary_part

!========================== REAL_O  ==========================================80
!
! Real part of arg_in.
!
!=============================================================================80

  pure function real_o( arg_in )

    use set_precision, only : dp

    real(dp), intent(in)  :: arg_in
    real(dp)              :: real_o

    continue

    real_o = real(arg_in,dp)

  end function real_o

!================================= CMPLX_O  ==================================80
!
! Real part of arg_in.
!
!=============================================================================80

  pure function cmplx_o( arg_in )

    use set_precision, only : dp

    complex(dp), intent(in)  :: arg_in
    real(dp)                 :: cmplx_o

    continue

    cmplx_o = real(arg_in,dp)

  end function cmplx_o

end module complex_functions
