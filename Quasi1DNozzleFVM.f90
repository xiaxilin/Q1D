
!This code will utilize the Finite Volume Method to 
!solve the quasi-1D nozzle using the Euler equations.
!Runge-Kutta time integration is used.
!There are six flux options found in Set_Inputs:
!Central difference with J-S-T damping
!The following use MUSCL interpolation of primitive variables
!Steger-Warming FVS
!Van Leer FVS (Needs to be fixed for shock case)
!Liou-Steffen (AUSM) FVS
!Roe's (Roe-Pike) FDS
!HLLC FDS (To be added)

!============================== set_precision ================================80
!
! Provides IEEE 754 compliant kinds
!
!=============================================================================80
module set_precision

  implicit none

  private

  public :: sngl, dbl, quad, dp

  integer, parameter :: sngl = selected_real_kind( 6, 37)
  integer, parameter :: dbl  = selected_real_kind(15, 307)
  integer, parameter :: quad = selected_real_kind(33, 4931)
  integer, parameter :: dp   = dbl

end module set_precision

!================================ constants ==================================80
!
! Provides useful constants
!
!=============================================================================80
module constants

  use set_precision, only : dp

  implicit none

  private

  public :: zero, one, two, three, four, five, six
  public :: sixth, fifth, fourth, third, half
  public :: small, large

  real(dp), parameter :: zero   = 0.0_dp
  real(dp), parameter :: one    = 1.0_dp
  real(dp), parameter :: two    = 2.0_dp
  real(dp), parameter :: three  = 3.0_dp
  real(dp), parameter :: four   = 4.0_dp
  real(dp), parameter :: five   = 5.0_dp
  real(dp), parameter :: six    = 6.0_dp

  real(dp), parameter :: sixth  = one/six
  real(dp), parameter :: fifth  = 0.2_dp
  real(dp), parameter :: fourth = 0.25_dp
  real(dp), parameter :: third  = one/three
  real(dp), parameter :: half   = 0.5_dp

  real(dp), parameter :: small  = tiny(one)
  real(dp), parameter :: large  = huge(one)

end module constants


!=============================================================================80
module set_inputs

  use set_precision, only : dp

  implicit none

  integer  :: i, j, k, n, z, RK
  integer  :: Num1stOrderIter, Flux_type, RKorder, nout, nmax, imax
  real(dp) :: R, gamma, Mref, conv, Toint, Point, Pback, kfour, ktwo,	       &
              kappa, eps, CFL, pi

! Sonic point fix parameter
  real(dp), parameter :: lambdaeps = 0.1_dp

end module set_inputs


!=============================================================================80
module solution_variables

  use set_precision, only : dp
  use set_inputs

  implicit none

  real(dp) :: T, P, psi, dt, VinMax, Rrho, Rrhou, Rrhoet
  real(dp) :: L1Rrho, L1Rrhou, L1Rrhoet, L1Rrhoinit, L1Rrhouinit, L1Rrhoetinit
  real(dp), allocatable, dimension(:)   :: a, nu, mach
  real(dp), allocatable, dimension(:,:) :: V, U, Unew, U0, Resid

  contains

  subroutine allocate_solution_variables(faces)

    implicit none

    integer, intent(in) :: faces
    
    allocate( a(faces+1), nu(faces+1), mach(faces+1) )
    allocate( V(3,faces+1), U(3,faces+1), Unew(3,faces+1), U0(3,faces+1) )
    allocate( Resid(3,faces+1) )

  end subroutine allocate_solution_variables

end module solution_variables

!=============================================================================80
module geometry

  use Set_precision, only : dp
  use Set_Inputs

  Implicit None

  Save

  real(dp) :: dx
  real(dp), allocatable, dimension(:) :: Areaface, Areacent, dAdx
  real(dp), allocatable, dimension(:) :: x
  real(dp), parameter                 :: NozzleLength = 2.0_dp

  contains

  subroutine allocate_geometry_variables(faces)

    implicit none

    integer, intent(in)::faces

    allocate( Areaface(faces+1), Areacent(faces+1), dAdx(faces+1), x(faces+1) )

  end subroutine allocate_geometry_variables

end module geometry

!=============================================================================80
module flux_constants

  use set_precision, only : dp
  use set_inputs

  implicit none

  save

  real(dp), allocatable, dimension(:,:) :: F,S,d
  real(dp), allocatable, dimension(:,:) :: Umin,Uplus,Vmin,Vplus

  contains

  subroutine allocate_flux_constants(faces)

    implicit none

    integer, intent(in)::faces

    allocate( F(3,faces), S(3,faces), d(3,faces) )
    allocate( Umin(3,faces+1), Uplus(3,faces+1) )
    allocate( Vmin(3,faces+1), Vplus(3,faces+1) )

  end subroutine allocate_flux_constants

end module flux_constants

!=============================================================================80
subroutine set_constants
  use set_precision, only : dp
  use constants,     only : one
  use set_inputs,    only : pi

  implicit none

  pi = acos(-one)

end subroutine set_constants

!=============================================================================80
subroutine read_input

  use set_precision, only : dp
  use set_inputs

  implicit none

  open(11,file="Quasi1DNozzleInput.txt",status='old')
  rewind(11)
  read(11, *) imax
  read(11, *) nmax
  read(11, *) nout
  read(11, *) RKorder
  read(11, *) Flux_type
  read(11, *) eps
  read(11, *) kappa
  read(11, *) Num1stOrderIter
  read(11, *) CFL
  read(11, *) conv
  read(11, *) ktwo
  read(11, *) kfour
  read(11, *) Mref
  read(11, *) Point
  read(11, *) Toint
  read(11, *) Pback
  read(11, *) gamma 
  read(11, *) R
  close(11)
  print *,'*********************************************'
  print *,'     Quasi-1D Nozzle : RELEASE 3.0, 2011     '
  print *,'*********************************************'
  print *,'Number of cells: ',imax
  print *,'Flux type:'
  if (Flux_type == 1) print*,'Central Difference w/ JST Damping'
  if (Flux_type == 2) print*,'Steger-Warming FVS'
  if (Flux_type == 3) print*,'Van Leer FVS'
  if (Flux_type == 4) print*,'AUSM'
  if (Flux_type == 5) print*,'Roe`s FDS'
  print *,'Runge-Kutta Order: ',RKorder
  print *,'CFL number: ',CFL
  print *,'MUSCL (0=No, 1=Yes): ',eps

end subroutine read_input

!=============================================================================80
subroutine primitive_to_conserved(faces, gamma, qp, Q)

  use set_precision, only : dp
  use constants,     only : half, one

  implicit none

  integer,  intent(in) :: faces
  real(dp), intent(in) :: gamma
  real(dp), dimension(3,faces+1), intent(in)    :: qp
  real(dp), dimension(3,faces+1), intent(inout) :: Q

  Q(1,:) = qp(1,:)
  Q(2,:) = qp(1,:)*qp(2,:)
  Q(3,:) = (qp(3,:)/(gamma-one))+half*(qp(1,:)*qp(2,:)*qp(2,:))

end subroutine primitive_to_conserved

!=============================================================================80
subroutine conserved_to_primitive(faces,gamma,Q,qp)

  use set_precision, only : dp
  use constants,     only : half, one

  implicit none

  integer,  intent(in) :: faces
  real(dp), intent(in) :: gamma
  real(dp), dimension(3,faces+1), intent(in)    :: Q
  real(dp), dimension(3,faces+1), intent(inout) :: qp

  qp(1,:) = Q(1,:)
  qp(2,:) = Q(2,:)/Q(1,:)
  qp(3,:) = (gamma-one)*(Q(3,:)-half*(Q(2,:)*Q(2,:)/Q(1,:)))

end subroutine conserved_to_primitive

!=============================================================================80
subroutine muscl_primitive_var(V, iterations, FirstOrderIters,                 &
                               epsin, kappa, faces, Vplus, Vmin)

  use set_precision, only : dp
  use constants,     only : zero, fourth, one

  implicit none

  integer,  intent(in) :: iterations,faces,FirstOrderIters
  real(dp), intent(in) :: epsin
  real(dp), intent(in) :: kappa
  real(dp), dimension(3,faces+1), intent(in) :: V
  real(dp), dimension(3,faces+1), intent(out) :: Vplus,Vmin

  integer  :: i
  real(dp) :: eps
  real(dp), dimension(3) :: Rmin, Rplus, Psimin, Psiplus

  continue

  if (iterations <= FirstOrderIters) then
    eps = zero
  else
    eps = epsin
  end if

  i=1
    Vplus(:,i) = V(:,i)+(fourth*eps)*((one+kappa)*(V(:,i+1)-V(:,i)))
    Vmin(:,i)  = zero

  do i=2,faces
    Rmin(:)   = (V(:,i)-V(:,i-1))/(V(:,i+1)-V(:,i)+0.000001_dp)
    Rplus(:)  = (V(:,i+1)-V(:,i))/(V(:,i)-V(:,i-1)+0.000001_dp)

    Psimin(:)  = one
!max(zero, (Rmin(:)+Rmin(:)**2)/(one+Rmin(:)**2))
    Psiplus(:) = one
!max(zero, (Rplus(:)+Rplus(:)**2)/(one+Rplus(:)**2))

    Vmin(:,i)  = V(:,i) - (fourth*eps)                                         &
               * ((one-kappa)*Psimin(:)  * (V(:,i+1) - V(:,i))                 &
               +  (one+kappa)*Psiplus(:) * (V(:,i)   - V(:,i-1)))
    Vplus(:,i) = V(:,i) + (fourth*eps)                                         &
               * ((one-kappa)*Psiplus(:) * (V(:,i)   - V(:,i-1))               &
               +  (one+kappa)*Psimin(:)  * (V(:,i+1) - V(:,i)))
  end do

  i=faces+1
    Vmin(:,i)  = V(:,i)-(fourth*eps)*((one+kappa)*(V(:,i)-V(:,i-1)))
    Vplus(:,i) = zero

! Apply floors to density and pressure
  Vmin(1,:)  = max(Vmin(1,:), 0.0001_dp)
  Vmin(3,:)  = max(Vmin(3,:), 500.0_dp)

  Vplus(1,:) = max(Vplus(1,:), 0.0001_dp)
  Vplus(3,:) = max(Vplus(3,:), 500.0_dp)

end subroutine muscl_primitive_var

!=============================================================================80
subroutine set_geometry

  use set_inputs
  use geometry

  implicit none

  call set_constants

  dx = NozzleLength/real(imax-1,dp)
  do i=1,imax+1
    x(i) = -0.5_dp*NozzleLength+dx*real(i-2,dp)
  end do

  do i=1,imax
    Areaface(i) = 0.2_dp+0.4_dp*(1.0_dp+sin(pi*(-1.0_dp+real(i-1,dp)*          &
                  dx-0.5_dp)))
    Areacent(i) = 0.2_dp+0.4_dp*(1.0_dp+sin(pi*(-1.0_dp+real(i-2,dp)*          &
                  dx+0.5_dp*dx-0.5_dp)))
    dAdx(i) = 0.4_dp*pi*cos(pi*(-1.0_dp+real(i-2,dp)*dx+0.5_dp*dx-0.5_dp))
  end do

  dAdx(1) = 0.0_dp
  Areacent(1) = Areacent(2)
  Areacent(imax+1) = Areacent(imax)

end subroutine set_geometry

!=============================================================================80
subroutine set_initial_values

  use set_inputs
  use solution_variables
  use flux_constants

  implicit none

!Set Initial Conditions
  psi = 1.0_dp+0.5_dp*(gamma-1.0_dp)*(Mref*Mref)
  T = Toint/psi
  P = Point/(psi**(gamma/(gamma-1.0_dp)))
!Initialize Primitive Variable array [rho U p]'
  do i=1,imax+1
    V(1,i) = P/(R*T)
    V(2,i) = Mref*sqrt(gamma*R*T)
    V(3,i) = P
    Vmin(:,i) = 0.0_dp
    Vplus(:,i) = 0.0_dp
  end do
!Initialize Conserved Variable array [rho rho*U rho*e0]'
  call primitive_to_conserved(imax,gamma,V,U)

!Initialize Unew array
  Unew = U

!Initialize Speed of Sound, mach, and smoothing array
  do i =1,imax+1
    a(i) = sqrt(gamma*V(3,i)/V(1,i))  !speed of sound a
    mach(i) = V(2,i)/a(i)
    nu(i) = 0.0_dp
  end do

end subroutine set_initial_values

!=============================================================================80
subroutine central_diff_with_jst_damping(faces,ktwo,kfour,a,qp,Q,F)

  use set_precision, only : dp
  use constants,     only : zero, half, two, three

  implicit none

  integer,  intent(in) :: faces
  real(dp), intent(in) :: ktwo,kfour
  real(dp), dimension(faces+1),   intent(in)  :: a
  real(dp), dimension(3,faces+1), intent(in)  :: qp,Q
  real(dp), dimension(3,faces),   intent(out) :: F

  integer :: i
  real(dp) :: epstwo,epsfour,lambda
  real(dp), dimension(faces+1) :: nu
  real(dp), dimension(3,faces) :: d

  do i=2,faces
    nu(i) = abs((qp(3,i-1) - two*qp(3,i) + qp(3,i+1))                          &
              / (qp(3,i-1) + two*qp(3,i) + qp(3,i+1)))
  end do
  nu(1)       = two*nu(2) - nu(3)
  nu(faces+1) = zero

!Calculate smoothing terms and dissipation vector
  i=1
    lambda  = half*(abs(qp(2,i+1))+a(i+1) + abs(qp(2,i))+a(i))
    epstwo  = ktwo*max(nu(i), nu(i+1), nu(i+2))
    epsfour = max(zero, kfour-epstwo)
    d(:,i)  = -lambda* (epstwo*(Q(:,i+1) - Q(:,i))                             &
              -epsfour*(Q(:,i+2) - three*Q(:,i+1) + two*Q(:,i)))
                
  do i=2,faces-1
    lambda  = half*(abs(qp(2,i+1))+a(i+1) + abs(qp(2,i))+a(i))
    epstwo  = ktwo*max(nu(i-1), nu(i), nu(i+1), nu(i+2))
    epsfour = max(zero, kfour-epstwo)
    d(:,i)  = -lambda*(epstwo*(Q(:,i+1)-Q(:,i))                                &
              -epsfour*(Q(:,i+2) - three*Q(:,i+1) + three*Q(:,i) - Q(:,i-1)))
  end do

  i=faces
    lambda  = half*(abs(qp(2,i+1))+a(i+1) + abs(qp(2,i))+a(i))
    epstwo  = ktwo*max(nu(i-1), nu(i), nu(i+1))
    epsfour = max(zero, kfour-epstwo)
    d(:,i)  = -lambda*(epstwo*(Q(:,i+1) - Q(:,i))                              &
              -epsfour*(-two*Q(:,i+1) + three*Q(:,i) - Q(:,i-1)))

!Update fluxes with smoothing terms  
  do i=1,faces
    F(1,i) = d(1,i) + half * (Q(2,i) + Q(2,i+1))
    F(2,i) = d(2,i) + half * ((Q(2,i)**2 / Q(1,i)) + qp(3,i) +                 &
                              (Q(2,i+1)**2 / Q(1,i+1)) +qp(3,i+1))
    F(3,i) = d(3,i) + half * ((Q(2,i)/Q(1,i))*(Q(3,i)+qp(3,i))+                &
                              (Q(2,i+1)/Q(1,i+1))*(Q(3,i+1)+qp(3,i+1)))
  end do

end subroutine central_diff_with_jst_damping

!=============================================================================80
subroutine steger_warming_fvs(faces,gamma,qplus,qmin,F)

  use set_precision, only : dp
  use constants,     only : half, one, two, three

  Implicit None

  integer,  intent(in) :: faces
  real(dp), intent(in) :: gamma
  real(dp), dimension(3,faces+1), intent(in)  :: qplus, qmin
  real(dp), dimension(3,faces),   intent(out) :: F

  integer  :: i
  real(dp) :: al, ar
  real(dp), dimension(3) :: lambdal, lambdar, lambdaplus, lambdamin
  real(dp), dimension(3,faces+1) :: Fiplus, Fimin

!Calculate Left (+) Fluxes
  do i=1,faces
    al = sqrt(gamma*qplus(3,i)/qplus(1,i))
    lambdal(1) = qplus(2,i)
    lambdal(2) = qplus(2,i)+al
    lambdal(3) = qplus(2,i)-al

    lambdaplus(:) = half*(lambdal(:)+abs(lambdal(:)))
    
    Fiplus(1,i) = two*(gamma-one)*lambdaplus(1) + lambdaplus(2) + lambdaplus(3)
    Fiplus(2,i) = two*(gamma-one)*                                             &
                  lambdaplus(1)*lambdal(1) +                                   &
                  lambdaplus(2)*lambdal(2) +                                   &
                  lambdaplus(3)*lambdal(3)
    Fiplus(3,i) = (gamma-one)*lambdaplus(1)*lambdal(1)**2                      &
                +half*(lambdaplus(2)*lambdal(2)**2+lambdaplus(3)*lambdal(3)**2)&
                + ((three-gamma)/(two*(gamma-one)))*                           &
                  (lambdaplus(2)+lambdaplus(3))*al**2
                    
    Fiplus(:,i) = (half*qplus(1,i)/gamma)*Fiplus(:,i)
  end do
!Calculate Right (-) Fluxes
  do i=2,faces+1
    ar = sqrt(gamma*qmin(3,i)/qmin(1,i))
    lambdar(1) = qmin(2,i)
    lambdar(2) = qmin(2,i)+ar
    lambdar(3) = qmin(2,i)-ar

    lambdamin(:) = half*(lambdar(:)-abs(lambdar(:)))
    
    Fimin(1,i) = two*(gamma-one)*lambdamin(1)+lambdamin(2)+lambdamin(3)
    Fimin(2,i) = two*(gamma-one)*                                              &
                 lambdamin(1)*lambdar(1) +                                     &
                 lambdamin(2)*lambdar(2) +                                     &
                 lambdamin(3)*lambdar(3)
    Fimin(3,i) = (gamma-one)*lambdamin(1)*lambdar(1)**2                        &
               + half*(lambdamin(2)*lambdar(2)**2 + lambdamin(3)*lambdar(3)**2)&
               + ((three-gamma)/(two*(gamma-one)))*                            &
                 (lambdamin(2)+lambdamin(3))*ar**2

    Fimin(:,i) = (half*qmin(1,i)/gamma)*Fimin(:,i)
  end do

!Calculate Interface Fluxes
  do i=1,faces
    F(:,i) = Fiplus(:,i)+Fimin(:,i+1)
  end do

end subroutine steger_warming_fvs

!=============================================================================80
subroutine van_leer_fvs(faces,gamma,qplus,qmin,F)

  use set_precision, only : dp
  use constants,     only : zero, fourth, half, one, two

  implicit none

  integer,  intent(in) :: faces
  real(dp), intent(in) :: gamma

  real(dp), dimension(3,faces+1), intent(in)    :: qplus, qmin
  real(dp), dimension(3,faces),   intent(inout) :: F

  integer  :: i
  real(dp) :: rho, a, M, Mfloor, fa, fb, switch
  real(dp), dimension(3) :: Fisub, Fiss
  real(dp), dimension(3,faces+1) :: Fiplus, Fimin

  continue

!Calculate Left (+) Fluxes
  do i = 1, faces
    rho = qplus(1,i)
    a = sqrt(gamma*qplus(3,i)/rho)
    M = qplus(2,i)/a

!Left sub(sonic) flux
    fa = fourth*rho*a*(M+one)**2
    fb = a*((gamma-one)*M + two)

    Fisub(1) = fa
    Fisub(2) = fa*fb / gamma
    Fisub(3) = fa*fb*fb / (two*(gamma**2-one))

!Floor the mach number for left supersonic flows
    Mfloor = half*(M + abs(M))
!Left supersonic flux
    Fiss(1) = rho*a*Mfloor
    Fiss(2) = (Mfloor/M)*rho*a**2*(Mfloor**2+one/gamma)
    Fiss(3) = rho*a**3*Mfloor*(half*Mfloor**2+one/(gamma-one))

!Combine sub and supersonic fluxes
    switch = real(max(0, 1-int(abs(M))),dp)
    Fiplus(:,i) = (one-switch)*Fiss(:) + switch*Fisub(:)
  end do

!Calculate Right (-) Fluxes
  do i=2,faces+1
    rho = qmin(1,i)
    a = sqrt(gamma*qmin(3,i)/rho)
    M = qmin(2,i)/a

!Right sub(sonic) flux
    fa = -fourth*rho*a*(M-one)**2
    fb = a*((gamma-one)*M - two)

    Fisub(1) = fa
    Fisub(2) = fa*fb / gamma
    Fisub(3) = fa*fb*fb / (two*(gamma**2-one))

!Floor the mach number
    Mfloor = half*(M - abs(M))

    Fiss(1) = rho*a*Mfloor
    Fiss(2) = (Mfloor/M)*rho*a**2*(Mfloor**2+one/gamma)
    Fiss(3) = rho*a**3*Mfloor*(half*Mfloor**2+one/(gamma-one))

!Combine sub and supersonic fluxes
    switch = real(max(0, 1-int(abs(M))),dp)
    Fimin(:,i) = (one-switch)*Fiss(:) + switch*Fisub(:) 
  end do
!Calculate Interface Fluxes
  do i=1,faces
    F(:,i) = Fiplus(:,i) + Fimin(:,i+1)
  end do

end subroutine van_leer_fvs

!=============================================================================80
subroutine ausm_fds(faces,gamma,qpplus,qpmin,Qplus,Qmin,F)

  use set_precision, only : dp

  implicit none

  integer,  intent(in) :: faces
  real(dp), intent(in) :: gamma

  real(dp), dimension(3,faces+1), intent(in)  :: qpplus, qpmin, Qplus, Qmin
  real(dp), dimension(3,faces),   intent(out) :: F

  integer  :: i
  real(dp) :: PL, PR, HTL, HTR, al, ar, Ml, Mr, Mlold, Mrold

  continue

  do i = 1,faces
!Calculate left (+) state
    al = sqrt(gamma*qpplus(3,i)/qpplus(1,i))
    Ml = qpplus(2,i)/al
    PL = qpplus(3,i)
    HTL = (Qplus(3,i)+PL)/Qplus(1,i)
    
    if (abs(Ml)<1.0_dp) then
      PL = PL*0.5_dp*(1.0_dp+Ml)
      Ml = 0.25_dp*(Ml+1.0_dp)**2
    else
      Mlold = Ml
      Ml = 0.5_dp*(Ml+abs(Ml))
      PL = PL*Ml/Mlold	
    end if
!Calculate right (-) state
    ar = sqrt(gamma*qpmin(3,i+1)/qpmin(1,i+1))
    Mr = qpmin(2,i+1)/ar
    PR = qpmin(3,i+1)
    HTR = (Qmin(3,i+1)+PR)/Qmin(1,i+1)

    if (abs(Mr)<1.0_dp) then
      PR = PR*0.5_dp*(1.0_dp-Mr)
      Mr = -0.25_dp*(Mr-1.0_dp)**2
    else
      Mrold = Mr
      Mr = 0.5_dp*(Mr-abs(Mr))
      PR = PR*Mr/Mrold
    end if
!Combine
    F(1,i) = 0.5_dp*(Ml+Mr)*(Qmin(1,i+1)*ar+Qplus(1,i)*al)-                    &
             0.5_dp*abs(Ml+Mr)*(Qmin(1,i+1)*ar-Qplus(1,i)*al)
    F(2,i) = 0.5_dp*(Ml+Mr)*(Qmin(2,i+1)*ar+Qplus(2,i)*al)-                    &
             0.5_dp*abs(Ml+Mr)*(Qmin(2,i+1)*ar-Qplus(2,i)*al)+(PL+PR)
    F(3,i) = 0.5_dp*(Ml+Mr)*(Qmin(1,i+1)*HTR*ar+Qplus(1,i)*HTL*al)-            &
             0.5_dp*abs(Ml+Mr)*(Qmin(1,i+1)*HTR*ar-Qplus(1,i)*HTL*al)
  end do

end subroutine ausm_fds

!=============================================================================80
subroutine roes_fds(faces,gamma,lambdaeps,qpplus,qpmin,Qplus,Qmin,F)

  use set_precision, only : dp
  use constants,     only : half, one, two, four

  implicit none

  integer,  intent(in) :: faces
  real(dp), intent(in) :: gamma, lambdaeps

  real(dp), dimension(3,faces+1), intent(in)  :: qpplus, qpmin, Qplus, Qmin
  real(dp), dimension(3,faces),   intent(out) :: F

  integer :: i,z
  real(dp) :: PL, PR, Rhalf, RoeAvgrho, RoeAvgu, RoeAvght, RoeAvga
  real(dp),dimension(3) :: FL, FR, r1RoeAvg, r2RoeAvg, r3RoeAvg, dw, lambdaRoe

  continue

  do i = 1,faces
    PL = qpplus(3,i)       
    PR = qpmin(3,i+1)

    FL(:) = (/Qplus(2,i),Qplus(2,i)**2/Qplus(1,i)+PL,   &
             (Qplus(2,i)/Qplus(1,i))*(Qplus(3,i)+PL)/)
    FR(:) = (/Qmin(2,i+1),Qmin(2,i+1)**2/Qmin(1,i+1)+PR,   &
             (Qmin(2,i+1)/Qmin(1,i+1))*(Qmin(3,i+1)+PR)/)

    Rhalf = sqrt(Qmin(1,i+1)/Qplus(1,i))                !Roe interface variable
    RoeAvgrho = Rhalf*Qplus(1,i)                                  !Avg Density
    RoeAvgu   = (Rhalf*Qmin(2,i+1)/Qmin(1,i+1)+                                &
                 Qplus(2,i)/Qplus(1,i))/(Rhalf+1)                 !Avg Velocity
    RoeAvght  = (Rhalf*(Qmin(3,i+1)+PR)/Qmin(1,i+1)+                           &
                (Qplus(3,i)+PL)/Qplus(1,i))/(Rhalf+1)             !Avg Enthalpy
    RoeAvga   = sqrt((gamma-1.0_dp)*(RoeAvght-0.5_dp*RoeAvgu**2)) !Avg a

    lambdaRoe(1) = RoeAvgu
    lambdaRoe(2) = RoeAvgu+RoeAvga
    lambdaRoe(3) = RoeAvgu-RoeAvga

!Entropy fix
    do z=1,3
      if ( abs(lambdaRoe(z)) <= two*lambdaeps*RoeAvga ) then
         lambdaRoe(z) = (lambdaRoe(z)**2)/(four*lambdaeps*RoeAvga) +         &
                        lambdaeps*RoeAvga;
      end if
    end do

    dw(1) = (Qmin(1,i+1)-Qplus(1,i))-(PR-PL)/RoeAvga**2
    dw(2) = (Qmin(2,i+1)/Qmin(1,i+1)-Qplus(2,i)/Qplus(1,i)) +                  &
            (PR-PL)/(RoeAvgrho*RoeAvga)
    dw(3) = (Qmin(2,i+1)/Qmin(1,i+1)-Qplus(2,i)/Qplus(1,i)) -                  &
            (PR-PL)/(RoeAvgrho*RoeAvga)

    r1RoeAvg(:) = (/one, RoeAvgu, half*RoeAvgu**2/)
    r2RoeAvg(:) = (half*RoeAvgrho/RoeAvga) *                                  &
                  (/one, RoeAvgu+RoeAvga, RoeAvght+RoeAvgu*RoeAvga/)
    r3RoeAvg(:) = (-half*RoeAvgrho/RoeAvga) *                                 &
                  (/one, RoeAvgu-RoeAvga, RoeAvght-RoeAvgu*RoeAvga/)

!Calculate Interface Fluxes
    F(:,i) = half*((FL(:)+FR(:))                                               &
           - (abs(lambdaRoe(1))*dw(1)*r1RoeAvg(:)                              &
           +  abs(lambdaRoe(2))*dw(2)*r2RoeAvg(:)                              &
           +  abs(lambdaRoe(3))*dw(3)*r3RoeAvg(:)))
  end do

end subroutine roes_fds
!=============================================================================80
!                         Begin Main Program Quasi1DNozzleFVM
!=============================================================================80
program quasi1dnozzlefvm

  use set_precision, only : dp
  use constants,     only : zero, one, two
  use set_inputs
  use solution_variables
  use geometry
  use flux_constants

  implicit none

  real(dp) :: entropy
  real(dp), dimension(3) :: ent

  continue

  call read_input()
!  call read_grid()

  call allocate_solution_variables(imax)
  call allocate_geometry_variables(imax)
  call allocate_flux_constants(imax)
  call set_initial_values
  call set_geometry

!Begin loop over max number of iterations
  do n=1,nmax

!Set Time Step

!!!! call set_time_step(dt,V,a)
    dt = 1000.0_dp
    do i=2,imax
      dt = min(((CFL*dx)/(abs(V(2,i))+a(i))),dt)
    end do

    do RK=1,RKorder

      select case (Flux_Type)
      case (1)
        call central_diff_with_jst_damping(imax,ktwo,kfour,a,V,U,F)
      case (2)
        call muscl_primitive_var(V,n,Num1stOrderIter,eps,kappa,imax,Vplus,Vmin)
        call steger_warming_fvs(imax,gamma,Vplus,Vmin,F)
      case (3)
        call muscl_primitive_var(V,n,Num1stOrderIter,eps,kappa,imax,Vplus,Vmin)
        call van_leer_fvs(imax,gamma,Vplus,Vmin,F)
      case (4)
        call muscl_primitive_var(V,n,Num1stOrderIter,eps,kappa,imax,Vplus,Vmin)
        call primitive_to_conserved(imax,gamma,Vplus,Uplus)
        call primitive_to_conserved(imax,gamma,Vmin,Umin)
        call ausm_fds(imax,gamma,Vplus,Vmin,Uplus,Umin,F)
      case (5)
        call muscl_primitive_var(V,n,Num1stOrderIter,eps,kappa,imax,Vplus,Vmin)
        call primitive_to_conserved(imax,gamma,Vplus,Uplus)
        call primitive_to_conserved(imax,gamma,Vmin,Umin)
        call roes_fds(imax,gamma,lambdaeps,Vplus,Vmin,Uplus,Umin,F)
      end select

!Source terms
      S = zero
      do i = 1,imax
        S(2,i) = V(3,i)*dAdx(i)
      end do

!M-Step Runge-Kutta Explicit Iteration
      if (RK.eq.1) then
        U0(:,:) = U(:,:)
      end if

      do i=2,imax
        Resid(:,i) = (one/(dx*Areacent(i))) *                                  &
               ( dx*S(:,i) - (F(:,i)*Areaface(i)) + (F(:,i-1)*Areaface(i-1)) )
        Unew(:,i) = U0(:,i)                                                    &
                  + dt*(one/(real(RKorder-RK,dp)+one))*Resid(:,i)
      end do

!Calculate new primitive variable vectors
      call conserved_to_primitive(imax,gamma,Unew,V)
      V(1,:) = max(V(1,:), 0.0001_dp)
      V(3,:) = max(V(3,:), 500.0_dp)

!Extrapolate primitive variables to first ghost cell
      VinMax = sqrt(two*gamma*R*Toint/(gamma-one))-one

!Extrapolated velocity limited by inflow conditions
      V(2,1) = max(-VinMax, min(two*V(2,2)-V(2,3), VinMax))	
      psi = Toint/(Toint-((gamma-one)*(V(2,1)*V(2,1))/(two*gamma*R)))
      V(1,1) = Point/(R*Toint*psi**(one/(gamma-one)))
      V(3,1) = Point/psi**(gamma/(gamma-one))
  
!Extrapolate interior values to end ghost cell
      V(1,imax+1) = max(two*V(1,imax)-V(1,imax-1), 0.0001_dp)
      V(2,imax+1) = two*V(2,imax)-V(2,imax-1)
      V(3,imax+1) = max(two*V(3,imax)-V(3,imax-1), 500.0_dp)

      if (Pback >= 0.00001_dp) then
        V(3,imax+1) = Pback
      end if

      do i=1,imax+1
        a(i) = sqrt(gamma*V(3,i)/V(1,i))
        mach(i) = V(2,i)/a(i)
      end do
  
!Updated conserved variables with limited primitive variables
      call primitive_to_conserved(imax,gamma,V,U)

    end do
!Calculate Residuals
    if (n.eq.1 .or. mod(n,nout).eq.0) then
      Rrho   = zero
      Rrhou  = zero
      Rrhoet = zero
      do i=2,imax
        Rrho   = Rrho+((U(1,i)-U0(1,i))/dt)**2
        Rrhou  = Rrhou+((U(2,i)-U0(2,i))/dt)**2
        Rrhoet = Rrhoet+((U(3,i)-U0(3,i))/dt)**2
      end do
      L1Rrho   = sqrt(amax1(Rrho, 1.e-20_dp))   / real(imax-1,dp)
      L1Rrhou  = sqrt(amax1(Rrhou, 1.e-20_dp))  / real(imax-1,dp)
      L1Rrhoet = sqrt(amax1(Rrhoet, 1.e-20_dp)) / real(imax-1,dp)

      if(n.eq.1) then
        L1Rrhoinit   = L1Rrho
        L1Rrhouinit  = L1Rrhou
        L1Rrhoetinit = L1Rrhoet
      end if
      L1Rrho   = L1Rrho   / L1Rrhoinit
      L1Rrhou  = L1Rrhou  / L1Rrhouinit
      L1Rrhoet = L1Rrhoet / L1Rrhoetinit
      write(*,300) n, L1Rrho, L1Rrhou, L1Rrhoet
300   format(1X,i8,2(e15.6),3(e15.6),4(e15.6))
      if (L1Rrho<=conv .and. L1Rrhou<=conv .and. L1Rrhoet<=conv) then
        goto 100
      end if
    end if

  end do

  write(*,*)'Solution failed to converge'
!  read(*,'()')
  open(40,file='q1Dnozzle.dat',status='unknown')
  write(40,*) 'TITLE = "Quasi-1D Nozzle Solution"'
  write(40,*) 'variables="x(m)""Area(m^2)""rho(kg/m^3)""u(m/s)"&
       &"Press(N/m^2)""Mach""U1""U2""U3""S1""S2""S3""entropy"'
  write(40,*) 'ZONE I=',imax+1
  write(40,*) 'DATAPACKING=POINT'
  write(40,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE &
       &DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
  do i=1,imax+1
    entropy = log(V(3,i)/V(1,i)**gamma)
    ent(1) = -entropy/(gamma-1.0_dp) + (gamma+1.0_dp)/(gamma-1.0_dp)-        &
         U(3,i)/V(1,i)/V(3,i)
    ent(2) = U(2,i)/V(3,i)
    ent(3) = -V(1,i)/V(3,i)
    write(40,*) x(i), Areacent(i), V(1,i), V(2,i), V(3,i),               &
         V(2,i)/a(i), U(1,i), U(2,i), U(3,i), ent(1), ent(2), ent(3), entropy
  end do
  stop
100 continue
  write(*,*)'Solution converged'
  open(40,file='q1Dnozzle.dat',status='unknown')
  write(40,*) 'TITLE = "Quasi-1D Nozzle Solution"'
  write(40,*) 'variables="x(m)""Area(m^2)""rho(kg/m^3)""u(m/s)"&
       &"Press(N/m^2)""Mach""U1""U2""U3""S1""S2""S3""entropy"'
  write(40,*) 'ZONE I=',imax+1
  write(40,*) 'DATAPACKING=POINT'
  write(40,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE &
       &DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
  do i=1,imax+1
    entropy = log(V(3,i)/V(1,i)**gamma)
    ent(1) = -entropy/(gamma-1.0_dp) + (gamma+1.0_dp)/(gamma-1.0_dp)-        &
         U(3,i)/V(1,i)/V(3,i)
    ent(2) = U(2,i)/V(3,i)
    ent(3) = -V(1,i)/V(3,i)
    write(40,*) x(i), Areacent(i), V(1,i), V(2,i), V(3,i),               &
         V(2,i)/a(i), U(1,i), U(2,i), U(3,i), ent(1), ent(2), ent(3), entropy
  end do
  stop

end program quasi1dnozzlefvm
