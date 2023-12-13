
module SHAKELie
   implicit none

   ! definition of integrator options
   type  :: SHAKELie_options
      ! system is constrained?
      integer :: constrained = 0
      ! use stabilized index two formulation
      integer :: stab2 = 0
      ! mass matrix is constant?
      integer :: const_mass_matrix = 0
      ! mass matrix is a diagonal matrix? If == 1 then SHAKELie_diag_M
      ! is used in order to calculate the mass matrix rather than
      ! SHAKELie_M
      integer :: diag_mass_matrix = 0
      ! iteration matrix is banded? (in contrained case: an ordered iteration matrix)
      integer :: banded_iteration_matrix = 0
      ! if iteration matrix is banded: number of subdiagonals
      integer :: nr_subdiag = 0
      ! if iteration matrix is banded: number of superdiagonals
      integer :: nr_superdiag = 0
      ! vector, that is used to order the iteration matrix in order to obtain a banded structure in the constrained case
      integer, dimension(:), allocatable  :: jour
      ! recalculate the iteration matrix in every Newton step
      integer :: recalc_iteration_matrix = 0
      ! use numerical approximation for Ct and Kt
      integer :: use_num_Ct = 1
      integer :: use_num_Kt = 1
      ! omit Ct and Kt resp. in the iteration matrix
      integer :: no_Ct = 0
      integer :: no_Kt = 0
      ! variables for error control of newton method (absolute and relativ tolerance)
      real(8)  :: atol = 1.0e-10_8
      real(8)  :: rtol = 1.0e-8_8
      integer :: imax   = 5
      ! integration span
      real(8) :: t0 = 0.0_8
      real(8) :: te = 1.0_8
      integer :: nsteps = 100
   end type SHAKELie_options

   ! definition of integrator statistics
   type  :: SHAKELie_statistics
      ! current number of newton steps TODO: private
      integer  :: newt_steps_curr = 0
      ! number of newton steps
      integer  :: newt_steps_sum = 0
      ! maximum number of newton steps
      integer  :: newt_steps_max = 0
      ! average number of newton steps
      real(8)  :: newt_steps_avg = 0.0_8
      ! number of calls
      integer  :: ngcalls = 0
      integer  :: nBcalls = 0
      ! integration time
      real(8)  :: time = 0.0_8
   end type SHAKELie_statistics

   ! definition of abstract problem type
   type, abstract :: SHAKELie_problem
      ! variables that define the state of the integrator
      real(8)                             :: t = 0.0_8
      integer                             :: sizeq = 0
      real(8), dimension(:), allocatable  :: q
      integer                             :: sizev = 0
      real(8), dimension(:), allocatable  :: v
      !real(8), dimension(:), allocatable  :: vd
      real(8), dimension(:), allocatable  :: Dq ! this is just for output
      real(8), dimension(:), allocatable  :: Dq_old
      integer                             :: sizel = 0 ! number of constraints, irrelevant for this%opts%constrained == 0
      real(8), dimension(:), allocatable  :: l ! Lagrange multipliers, only needed in the constrained case
      real(8), dimension(:), allocatable  :: eta ! Auxiliar variables eta, only needed in the stabilized index-2 case
      ! integrator options
      type(SHAKELie_options)                :: opts
      ! solver statistics
      type(SHAKELie_statistics)             :: SHAKELie_stats
      ! constant mass matrix (if opts%const_mass_matrix == 1)
      real(8), dimension(:,:), allocatable   :: SHAKELie_const_M
      ! constant diagonal mass matrix (if opts%diag_mass_matrix == 1, also)
      real(8), dimension(:), allocatable     :: SHAKELie_const_diag_M
      ! internal variables:
      ! SHAKELie uses the impulse instead of velocity
      real(8), dimension(:), allocatable  :: p
   ! definition of deferred procedures
   contains
         ! function $M$ in $M(q) \dot v = \hat v^T M v - f(q,v,t)$
      procedure(SHAKELie_M),              deferred :: SHAKELie_M
         ! function $M$ in $M(q) \dot v = \hat v^T M V - f(q,v,t)$, diagonal case
      procedure(SHAKELie_diag_M),         deferred :: SHAKELie_diag_M
         ! function $g$ in $M(q) \dot v = \hat v^T M V - f(q,v,t)$
      procedure(SHAKELie_f),              deferred :: SHAKELie_f
         ! operation in the lie space $q_n * \exp(h\cdot\widetilde{\Delta q_n})$
      procedure(SHAKELie_qlpexphDqtilde), deferred :: SHAKELie_qlpexphDqtilde
         ! operation in the lie algebra $inversetilde([\tilde{v},\tilde{w}])$, may be a dummy, if system is unconstrained
      procedure(SHAKELie_itlbtvtw),       deferred :: SHAKELie_itlbtvtw
         ! tangent damping matrix $C_t$
      procedure(SHAKELie_Ct),             deferred :: SHAKELie_Ct
         ! tangent stiffness matrix $K_t$
      procedure(SHAKELie_Kt),             deferred :: SHAKELie_Kt
         ! tangent stiffness matrix $K_t$ in the constrained case
         ! (may depend on the Lagrange multipliers)
      procedure(SHAKELie_Kt_lambda),      deferred :: SHAKELie_Kt_lambda
         ! tangent operator of the exponential map
      procedure(SHAKELie_Tg),             deferred :: SHAKELie_Tg
         ! transposed inverse of the tangent operator of the exponential map
      procedure(SHAKELie_Tg_inv_T),       deferred :: SHAKELie_Tg_inv_T
         ! derivative of transposed inverse of the tangent operator of the exponential map
      procedure(SHAKELie_d_Tg_inv_T),     deferred :: SHAKELie_d_Tg_inv_T
         ! subroutine for output, is called after every integration step
      procedure(SHAKELie_outputFunction), deferred :: SHAKELie_outputFunction
         ! subroutine to initialize the problem
      procedure(SHAKELie_init),           deferred :: SHAKELie_init
         ! function $\Phi(q)$
      procedure(SHAKELie_phi),            deferred :: SHAKELie_phi
         ! function $B(q)$, where $B(q)v = d/dq \Phi(q) * (DL_q(e) v)$
      procedure(SHAKELie_b),              deferred :: SHAKELie_b
         ! function $F(q,v) = d/(dq) (B(q)v) * (DL_q(e) v)$
      procedure(SHAKELie_Z),              deferred :: SHAKELie_Z
         ! function Z, where $Zw = Z(q(Dq)) (v(Dq + BT(q)eta), T w)$
      procedure(SHAKELie_matZ),           deferred :: SHAKELie_matZ
         ! subroutine to calculate correct initial values for vd and a
      procedure                  :: SHAKELie_calcInitial => SHAKELie_calcInitial ! TODO: private
         ! subroutine to calculate correct initial values for vd and a
      procedure                  :: SHAKELie_calcInitialConstrained => SHAKELie_calcInitialConstrained
         ! subroutine to integrate one time step
      procedure                  :: SHAKELie_solveTimeStep => SHAKELie_solveTimeStep ! TODO: private
         ! subroutine to integrate one time step in a constrained system
      procedure                  :: SHAKELie_solveConstrainedTimeStep => SHAKELie_solveConstrainedTimeStep
         ! subroutine to integrate one time step in a constrained system
      procedure                  :: SHAKELie_solveConstrainedTimeStep_stab2 => SHAKELie_solveConstrainedTimeStep_stab2
         ! subroutine for integration
      procedure                  :: SHAKELie_integrate => SHAKELie_integrate
         ! subroutine for numerical approximation of Ct
      procedure                  :: SHAKELie_num_Ct => SHAKELie_num_Ct
         ! subroutine for numerical approximation of Kt
      procedure                  :: SHAKELie_num_Kt => SHAKELie_num_Kt
         ! subroutine for numerical approximation of Kt in the constrained case
      procedure                  :: SHAKELie_num_Kt_lambda => SHAKELie_num_Kt_lambda
         ! In order to use the numerical approximations of Ct and Kt
         ! respectively, the procedures SHAKELie_Ct and SHAKELie_Kt as well as
         ! their constrained pendants gen_CT_lambda and SHAKELie_Kt_lambde
         ! have to be referenced to SHAKELie_num_Ct and SHAKELie_num_Kt resp.
      procedure                  :: SHAKELie_print_stats => SHAKELie_print_stats
         ! Clean up the SHAKELie_problem object (only internal variables
         ! are reset, SHAKELie_alpha_m, SHAKELie_alpha_f, SHAKELie_beta and
         ! SHAKELie_gamma are NOT reset)
      procedure                  :: SHAKELie_cleanup => SHAKELie_cleanup
   end type SHAKELie_problem

   ! abstract interface of the problem
   abstract interface
      pure function SHAKELie_M(this, q) result(rslt)
         import                                          :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),   intent(in)            :: this
         real(8), dimension(:), intent(in)               :: q
         ! result
         real(8), dimension(this%sizev,this%sizev)       :: rslt
      end function SHAKELie_M

      pure function SHAKELie_diag_M(this, q) result(rslt)
         import                                          :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),   intent(in)            :: this
         real(8), dimension(:), intent(in)               :: q
         ! result
         real(8), dimension(this%sizev)                  :: rslt
      end function SHAKELie_diag_M

      pure function SHAKELie_f(this, q, v, t) result(rslt)
         import                              :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8), dimension(:), intent(in)   :: v
         real(8),               intent(in)   :: t
         ! result
         real(8), dimension(this%sizev)      :: rslt
      end function SHAKELie_f

      pure function SHAKELie_qlpexphDqtilde(this, q, h, dq) result(rslt)
         import                              :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8),               intent(in)   :: h
         real(8), dimension(:), intent(in)   :: dq
         ! result
         real(8), dimension(this%sizeq)      :: rslt
      end function SHAKELie_qlpexphDqtilde

      pure function SHAKELie_itlbtvtw(this, v, w) result(rslt)
         import                              :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         real(8), dimension(:), intent(in)   :: w
         ! result
         real(8), dimension(this%sizev)      :: rslt
      end function SHAKELie_itlbtvtw

      pure function SHAKELie_Ct(this, q, v, t) result(rslt)
         import                                    :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),   intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function SHAKELie_Ct

      pure function SHAKELie_Kt(this, q, v, vd, t) result(rslt)
         import                                    :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),   intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8), dimension(:), intent(in)         :: vd
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function SHAKELie_Kt

      pure function SHAKELie_Kt_lambda(this, q, v, vd, l, t) result(rslt)
         import                                    :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),   intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8), dimension(:), intent(in)         :: vd
         real(8), dimension(:), intent(in)         :: l
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function SHAKELie_Kt_lambda

      pure function SHAKELie_Tg(this, h, dq) result(rslt)
         import                                    :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),intent(in)         :: this
         real(8),               intent(in)         :: h
         real(8), dimension(:), intent(in)         :: dq
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function SHAKELie_Tg

      pure function SHAKELie_Tg_inv_T(this, dq) result(rslt)
         import                                    :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),   intent(in)      :: this
         real(8), dimension(:), intent(in)         :: dq
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function SHAKELie_Tg_inv_T

      pure function SHAKELie_d_Tg_inv_T(this, v, w) result(rslt)
         import                                    :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),   intent(in)      :: this
         real(8), dimension(:), intent(in)         :: v, w
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function SHAKELie_d_Tg_inv_T

      pure function SHAKELie_norm(this, v) result(rslt)
         import                  :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         ! result
         real(8)                             :: rslt
      end function SHAKELie_norm

      subroutine SHAKELie_outputFunction(this,info)
         import                           :: SHAKELie_problem
         ! input
         class(SHAKELie_problem), intent(in)  :: this
         integer,             intent(in)  :: info
      end subroutine SHAKELie_outputFunction

      subroutine SHAKELie_init(this)
         import                              :: SHAKELie_problem
         ! input/output
         class(SHAKELie_problem), intent(inout)  :: this
      end subroutine SHAKELie_init

      pure function SHAKELie_phi(this,q) result(rslt)
         import                                       :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         ! result
         real(8), dimension(this%sizel)               :: rslt
      end function SHAKELie_phi

      pure function SHAKELie_B(this,q) result(rslt)
         import                                       :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         ! result
         real(8), dimension(this%sizel,this%sizev)    :: rslt
      end function SHAKELie_B

      pure function SHAKELie_Z(this,q,v) result(rslt)
         import                                       :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         real(8), dimension(:),            intent(in) :: v
         ! result
         real(8), dimension(this%sizel)               :: rslt
      end function SHAKELie_Z

      pure function SHAKELie_matZ(this,q,v,T) result(rslt)
         import                                       :: SHAKELie_problem
         ! input
         class(SHAKELie_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         real(8), dimension(:),            intent(in) :: v
         real(8), dimension(:,:),          intent(in) :: T
         ! result
         real(8), dimension(this%sizel, this%sizev)   :: rslt
      end function SHAKELie_matZ

   end interface

   contains

   ! subroutine for pertubing initial values
   subroutine SHAKELie_calcInitial(this, h)
      implicit none
      class(SHAKELie_problem),  intent(inout)  :: this  ! problem object
      real(8),                intent(in   )  :: h     ! step size
      !
      ! Nothing to do
   end subroutine SHAKELie_calcInitial

   subroutine SHAKELie_calcInitialConstrained(this, h) !TODO TODO
      implicit none
      class(SHAKELie_problem),  intent(inout)  :: this  ! problem object
      real(8),                intent(in   )  :: h     ! step size
      ! internal variables
      integer                                                              :: i
      real(8), dimension(this%sizev+this%sizel)                            :: vdl
      integer, dimension(this%sizev+this%sizel)                            :: ipiv0   ! pivot vector for dgesv for LU factorization of MBB0
      integer                                                              :: info    ! info flag for dgesv
      real(8), dimension(this%sizev+this%sizel, this%sizev+this%sizel)     :: MBB0    ! Matrix
      real(8), dimension(this%sizel, this%sizev)                           :: B0      ! $B(q_0)$
      !

      ERROR STOP 'constrained not yet supported'
#if 0
      associate ( t  => this%t,            &
                  q  => this%q,            &
                  v  => this%v,            &
                  vd => this%vd,           &
                  l  => this%l,            &
                  sv => this%sizev,        &
                  sl => this%sizel         )
         ! calulate $B(q_0)$
         B0 = this%SHAKELie_b(q)
         ! count calls
         this%SHAKELie_stats%nBcalls = this%SHAKELie_stats%nBcalls + 1

         ! calculate $MBB0$
         if (this%opts%diag_mass_matrix == 1) then
            ! Set the mass matrix part to zero beforehand
            MBB0(1:sv, 1:sv) = 0.0_8
            if (this%opts%const_mass_matrix == 1) then
               MBB0(1:sv, 1) = this%SHAKELie_const_diag_M
            else
               MBB0(1:sv, 1) = this%SHAKELie_diag_M(q)
            end if
            do concurrent (i=2:sv)
               MBB0(i,i) = MBB0(i,1)
               MBB0(i,1) = 0.0_8
            end do
         else
            if (this%opts%const_mass_matrix == 1) then
               MBB0(1:sv, 1:sv) = this%SHAKELie_const_M
            else
               MBB0(1:sv, 1:sv) = this%SHAKELie_M(q)
            end if
         end if
         MBB0(sv+1:sv+sl, 1:sv) =           B0
         MBB0(1:sv, sv+1:sv+sl) = transpose(B0)
         MBB0(sv+1:sv+sl, sv+1:sv+sl) = 0.0_8

         ! we need to solve a linear equation, first calulate the rhs
         vdl(1:sv)       = -this%SHAKELie_g(q, v, t)
         vdl(sv+1:sv+sl) = -this%SHAKELie_Z(q, v)
         ! count calls
         this%SHAKELie_stats%ngcalls = this%SHAKELie_stats%ngcalls + 1
         ! then solve the system
         call dgesv(          &! solve the System A*X=B and save the result in B
                     sv+sl,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                     1,       &! number of right hand sides (=size(B,2))
                     MBB0,    &! matrix A
                     sv+sl,   &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                     ipiv0,   &! integer pivot vector; it is not needed
                     vdl,     &! matrix B
                     sv+sl,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                     info)     ! integer information flag
         ! Now vdl actually contains $\dot v(t_0)$ and $\lambda(t_0)
         if (info .ne. 0)  print*, "dgesv sagt info=", info ! TODO

         ! apply the calculated values
         vd = vdl(   1:   sv)
         l  = vdl(sv+1:sv+sl)

      end associate
#endif
   end subroutine SHAKELie_calcInitialConstrained

   ! subroutine for integrating one time step
   subroutine SHAKELie_solveConstrainedTimeStep(this,t1)
      implicit none
      class(SHAKELie_problem), intent(inout)  :: this  ! problem object
      real(8)            , intent(in   )  :: t1    ! $t_{n+1}$

      ! internal integer variables
      integer                                      :: i     ! for iteration
      integer, dimension(this%sizev+this%sizel)    :: ipiv  ! needed for dgesv
      integer                                      :: info  ! needed for dgesv

      ! intenal real variables
      real(8)                                                           :: h       ! step size
      real(8), dimension(this%sizev)                                    :: vd1     ! $\dot v_{n+1}$
      !real(8), dimension(this%sizev)                                    :: a1      ! $a_{n+1}$ not needed 2015-06-23
      real(8), dimension(this%sizev)                                    :: v1      ! $v_{n+1}$
      real(8), dimension(this%sizev)                                    :: Dq      ! $\Delta q_n$
      real(8), dimension(this%sizeq)                                    :: q1      ! $q_{n+1}}
      real(8), dimension(this%sizev+this%sizel)                         :: res     ! $res(\dots)$
      real(8), dimension(this%sizev+this%sizel,this%sizev+this%sizel)   :: St      ! $S_t$
      real(8), dimension(:,:), allocatable                              :: Stb     ! permuted $S_t$ in banded format
      !real(8)                                                           :: beta_d  ! $\beta'$ not needed 2015-06-23
      !real(8)                                                           :: gamma_d ! $\gamma'$ not needed 2015-06-23
      real(8), dimension(this%sizev+this%sizel)                         :: Dxl     ! $(\Delta x, \Delta \lambda)^T$

      ! internal logical
      logical                                                           :: converged
      ERROR STOP "not yet implemeted"
#if 0
      ! associate construct for better readability
      associate (alf => this%SHAKELie_alpha_f,    &
                 alm => this%SHAKELie_alpha_m,    &
                 bet => this%SHAKELie_beta,       &
                 gam => this%SHAKELie_gamma,      &
                 sv  => this%sizev,           &
                 sl  => this%sizel,           &
                 v   => this%v,               &
                 vd  => this%vd,              &
                 a   => this%a,               &
                 q   => this%q,               &
                 t   => this%t,               &
                 l   => this%l)
         ! calculation of step size $h$
         h = t1 - t

#ifdef ZEROINIT
         Dq = 0.0_8
         v1 = (1 - gam/bet)*v + h*(1 - gam/(2*bet))*a
         vd1 = (1-alm)/(bet*(1-alf))*( (Dq - v)/h - 0.5*a)  +  (a - alf*vd)/(1-alf)
         l = 0.0_8
#else
         ! initialization of the new values $\dot v_{n+1}$, $a_{n+1}$ and $v_{n+1}$ and the value $\Delta q_n$ and $\lambda_{n+1}$
         !vd1 = 0.0_8 ! vd ! 0.0_8 DEBUG
         !a1  = (alf*vd - alm*a)/(1.0_8 - alm)
         !v1  = v + h*(1.0_8 - gam)*a + gam*h*a1
         !Dq  = v + (0.5_8 - bet)*h*a + bet*h*a1
         !!l   = l ! DEBUG
         !l   = 0.0_8
         Dq  = v + 0.5*h*a
         v1  = v +     h*a
         vd1 = (a - alf*vd)/(1-alf)
#endif


         ! loop of the newton method
         converged = .false.
         do i=1,this%opts%imax
            !! DEBUG
            !print *, t, i
            !print *, Dq
            !! GUBED
            !print *, 'Newton: ', i, this%opts%imax DEBUG
            ! calculate the new value $q_{n+1}$
            q1  = this%SHAKELie_qlpexphDqtilde(q,h,Dq)

            ! caluclate the residue $res$
            res(1:sv) = this%SHAKELie_g(q1,v1,t1) + matmul(transpose(this%SHAKELie_B(q1)),l)
            ! count calls
            this%SHAKELie_stats%ngcalls = this%SHAKELie_stats%ngcalls + 1
            this%SHAKELie_stats%nBcalls = this%SHAKELie_stats%nBcalls + 1
            !
            if (this%opts%const_mass_matrix == 1) then
               if (this%opts%diag_mass_matrix == 1) then
                  res(1:sv) = res(1:sv) + this%SHAKELie_const_diag_M*vd1
               else
                  res(1:sv) = res(1:sv) + matmul(this%SHAKELie_const_M,vd1)
               end if
            else
               if (this%opts%diag_mass_matrix == 1) then
                  res(1:sv) = res(1:sv) + this%SHAKELie_diag_M(q1)*vd1
               else
                  res(1:sv) = res(1:sv) + matmul(this%SHAKELie_M(q1),vd1)
               end if
            end if
            ! scale
            res(1:sv) = h*res(1:sv)
            res(sv+1:sv+sl) = this%SHAKELie_phi(q1)/h

            !! check if norm of the residue is sufficiently small to stop the iteration
            !if ((this%SHAKELie_norm(res(1:sv))       < this%opts%tolr  )  .and. &
            !    (this%SHAKELie_norm(res(sv+1:sv+sl)) < this%opts%tolphi)) then
            !   exit
            !end if

            ! solve the linear System $S_t \cdot \Delta xl = -res$
            Dxl = -res ! the result will be saved in the right hand side's spot

            if (i == 1 .or. this%opts%recalc_iteration_matrix == 1) then
               ! calculate iteration matrix
               if (this%opts%const_mass_matrix == 1) then
                  if (this%opts%diag_mass_matrix == 1) then
                     St(1:sv,1:sv) = 0.0_8
                     forall (i=1:sv)
                        St(    i   ,     i   ) = this%SHAKELie_const_diag_M(i) * (1-alm)/(bet*(1-alf))
                     end forall
                  else
                     St(   1:sv   ,   1:sv   ) = this%SHAKELie_const_M * (1-alm)/(bet*(1-alf))
                  end if
               else
                  if (this%opts%diag_mass_matrix == 1) then
                     St(1:sv,1:sv) = 0.0_8
                     St(   1:sv   ,   1   ) = this%SHAKELie_diag_M(q1) * (1-alm)/(bet*(1-alf))
                     forall (i=2:sv)
                        St(i,i) = St(i,1)
                        St(i,1) = 0.0_8
                     end forall
                  else
                     St(   1:sv   ,   1:sv   ) = this%SHAKELie_M(q1) * (1-alm)/(bet*(1-alf))
                  end if
               end if
               St(   1:sv   ,sv+1:sv+sl) = transpose(this%SHAKELie_B(q1)) ! TODO: dont calculate twice
               St(sv+1:sv+sl,   1:sv   ) = matmul(this%SHAKELie_B(q1), this%SHAKELie_Tg(h,Dq))
               St(sv+1:sv+sl,sv+1:sv+sl) = 0.0_8
               ! count calls
               this%SHAKELie_stats%nBcalls = this%SHAKELie_stats%nBcalls + 2

               if (this%opts%no_Ct == 0) then
                  if (this%opts%use_num_Ct == 1) then
                     St(1:sv,1:sv) = St(1:sv,1:sv) + this%SHAKELie_num_Ct(q1,v1,t1) * h*gam/bet
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%SHAKELie_stats%ngcalls = this%SHAKELie_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%SHAKELie_stats%ngcalls = this%SHAKELie_stats%ngcalls + this%sizev + 1
                     end if
                  else
                     St(1:sv,1:sv) = St(1:sv,1:sv) + this%SHAKELie_Ct(q1,v1,t1) * h*gam/bet
                  end if
               end if

               if (this%opts%no_Kt == 0) then
                  if (this%opts%use_num_Kt == 1) then
                     St(1:sv,1:sv) = St(1:sv,1:sv) + matmul(this%SHAKELie_num_Kt_lambda(q1,v1,vd1,l,t1),this%SHAKELie_Tg(h,Dq)) * h*h ! TODO: dgemm oder sowas
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%SHAKELie_stats%ngcalls = this%SHAKELie_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                        this%SHAKELie_stats%nBcalls = this%SHAKELie_stats%nBcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%SHAKELie_stats%ngcalls = this%SHAKELie_stats%ngcalls + this%sizev + 1
                        this%SHAKELie_stats%nBcalls = this%SHAKELie_stats%nBcalls + this%sizev + 1
                     end if
                  else
                     St(1:sv,1:sv) = St(1:sv,1:sv) + matmul(this%SHAKELie_Kt_lambda(q1,v1,vd1,l,t1),this%SHAKELie_Tg(h,Dq)) * h*h ! TODO: dgemm oder sowas
                  end if
               end if

               !! DEBUG !
               !print *, "t=", t
               !print *, "cond=", mycond(St)
               !print *, "maxstep=", this%SHAKELie_stats%newt_steps_max
               !!print *,
               !call print_vector(q,'qn')
               !call print_vector(v,'vn')
               !call print_vector(vd,'vdn')
               !call print_vector(a,'an')
               !call print_vector(l,'ln')
               !call print_vector(res,'resn')
               !call print_matrix(St,'Stn')
               !errorstop "mopp"
               !if (t > 2.5) then
               !   call print_matrix(St,'St')
               !   !call print_vector_int(this%opts%jour,'jour')
               !   stop
               !end if
               !! GUBED !

               !!DEBUG
               ! call print_matrix(St, 'Stb')
               !!GUBED

                        !!print *, St ! DEBUG
               ! Actually solve the linear System $S_t \cdot \Delta xl = -res$
               if (this%opts%banded_iteration_matrix == 1) then
#ifdef DEBUG_PRINT_ITERATION_MATRIX_AT
                  if (t >= DEBUG_PRINT_ITERATION_MATRIX_AT) then
                     call print_matrix(St,'St')
                     call print_vector_int(this%opts%jour,'jour')
                     stop "Successfully printed iteration matrix"
                  end if
#endif
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag)
                     if (.not. allocated(Stb)) then
                        allocate(Stb(lo*2+hi+1, sv+sl))
                     end if
                     !Stb = 0.0_8 ! DEBUG
                     forall (i=1:sv+sl)
                        Stb((max(1, i-hi)-i+lo+1+hi):(min(sv+sl,i+lo)-i+lo+1+hi),i) &
                            = St(this%opts%jour(max(1, i-hi):min(sv+sl, i + lo)), this%opts%jour(i))
                     end forall
                     Dxl = Dxl(this%opts%jour)
                     call dgbsv(           &! solve the system A*X=B and save the result in B
                                 sv+sl,    &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,       &!
                                 hi,       &!
                                 1,        &! number of right hand sides (=size(B,2))
                                 Stb,      &! matrix A
                                 2*lo+hi+1,&! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv,     &! integer pivot vector; it is not needed
                                 Dxl,      &! matrix B
                                 sv+sl,    &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)      ! integer information flag
                     Dxl(this%opts%jour) = Dxl
                     !call print_vector(Dxl,'Dxl_b') ! DEBUG
                     !stop ! DEBUG
                  end associate
               else
                  call dgesv(       &! solve the System A*X=B and save the result in B
                              sv+sl,&! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                              1,    &! number of right hand sides (=size(B,2))
                              St,   &! matrix A
                              sv+sl,&! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv, &! integer pivot vector; it is not needed
                              Dxl,  &! matrix B
                              sv+sl,&! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                              info)  ! integer information flag
                  !call print_vector(Dxl,'Dxl') ! DEBUG
                  !stop ! DEBUG
                  ! TODO: Errors abfragen
               end if

               if (info.ne.0) print *, "dgesv sagt info=", info

               ! update $\dot v_{n+1}$, $v_{n+1}$ and $\Delta q_n$
            else ! (St was not (re)calculated)
               ! Actually solve the linear System $S_t \cdot \Delta xl = -res$
               if (this%opts%banded_iteration_matrix == 1) then
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag, &
                             uStb => ubound(Stb,1))
                     Dxl = Dxl(this%opts%jour) !TODO
                     call dgbtrs(               &
                                 'No transpose',&! solve the System A*X=B and save the result in B
                                 sv+sl,         &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,            &!
                                 hi,            &!
                                 1,             &! number of right hand sides (=size(B,2))
                                 Stb,           &! matrix A
                                 2*lo+hi+1,     &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv,          &! integer pivot vector; it is not needed
                                 Dxl,           &! matrix B
                                 sv+sl,         &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)           ! integer information flag
                     Dxl(this%opts%jour) = Dxl
                  end associate
               else
                  call dgetrs(                &
                              'No transpose', & ! solve the System A*X=B and save the result in B
                              sv + sl,        & ! number of linear equations (=size(A,1))  ! Vorsicht: double precision muss real(8) sein
                              1,              & ! number of right hand sides (=size(B,2))
                              St,             & ! matrix A
                              sv + sl,        & ! leading dimesnion of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv,           & ! integer pivot vector; holds the pivoting that was calculated in the dgesv call
                              Dxl,            & ! matrix B
                              sv + sl,        & ! leading dimension of B, in thie case is equal to the number of linear equations (=size(A,1))
                              info)             ! integer information flag
               end if
               ! DEBUG
               if (info.ne.0) print*, 'dgetrs sagt info=', info
               ! GUBED
            end if

            !! DEBUG
            !call print_vector(Dxl,'Dxln')
            !stop
            !! GUBED

            ! Check for convergence
            if ( ( sum((Dxl(   1:   sv) / (this%opts%atol + this%opts%rtol*abs(Dq )))**2 )  &
                  +sum((Dxl(sv+1:sv+sl) / (this%opts%atol + this%opts%rtol*abs(h*l)))**2 ) )&
                   / (sv+sl)   <=  1.0_8 ) then
               converged = .true.
            end if

            ! update $\dot v_{n+1}$, $v_{n+1}$, $\Delta q_n$ and $\lambda_{n+1}$
#define REL *1
            Dq  = Dq  + Dxl(1:sv)                              REL
            v1  = v1  + gam/bet * Dxl(1:sv)                    REL
            vd1 = vd1 + (1-alm)/(bet*(1-alf)*h) * Dxl(1:sv)    REL
            l   = l   + Dxl(sv+1:sv+sl)/h                      REL
#undef REL

            ! update solver stats
            this%SHAKELie_stats%newt_steps_curr = i

            ! exit loop, if converged
            if ( converged ) exit

         end do ! end newton method

         ! DEBUG
         if (this%SHAKELie_stats%newt_steps_curr >= this%opts%imax) then
            print *, ' ' ! DEBUG
            print *, "Newton did not converge at t=", t
            !! DEBUG !
            call this%SHAKELie_print_stats
            print *, "Iteration matrix:"
            call print_matrix(St,'St')
            errorstop
            !! GUBED !
         end if


         ! The integration step is done; now the new values have to be
         ! applied to the problem object, and the new value for a can
         ! be calculated
         a  = ((1-alf)*vd1 + alf*vd - alm*a) / (1-alm)
         t  = t1
         q  = q1
         v  = v1
         vd = vd1

      end associate
#endif
   end subroutine SHAKELie_solveConstrainedTimeStep

   ! subroutine for integrating one time step
   subroutine SHAKELie_solveConstrainedTimeStep_stab2(this,t1)
      implicit none
      class(SHAKELie_problem), intent(inout)  :: this  ! problem object
      real(8)            , intent(in   )  :: t1    ! $t_{n+1}$

      ! internal integer variables
      integer                                    :: i     ! for iteration
      integer, dimension(this%sizev+2*this%sizel)  :: ipiv  ! needed for dgesv
      integer                                    :: info  ! needed for dgesv

      ! intenal real variables
      real(8)                                                          :: h       ! step size
      real(8), dimension(this%sizeq)                                   :: q1      ! $q_{k+1}}
      real(8), dimension(this%sizev+2*this%sizel)                      :: res     ! $res(\dots)$
      real(8), dimension(this%sizel,this%sizev)                        :: Bq
      real(8), dimension(this%sizev)                                   :: Dqs
      real(8), dimension(this%sizev)                                   :: M_Dq
      real(8), dimension(this%sizev)                                   :: M_Dq_old
      real(8), dimension(this%sizev+2*this%sizel,this%sizev+2*this%sizel)  :: St      ! $S_t$
      real(8), dimension(:,:), allocatable                             :: Stb     ! permuted $S_t$ in banded format
      real(8), dimension(this%sizev+2*this%sizel)                        :: Dxi   ! $(\Delta x, \Delta \lambda, \Delta \eta)^T$

      ! internal logical
      logical                                                              :: converged


      ! associate construct for better readability
      associate (sv  => this%sizev,           &
                 sl  => this%sizel,           &
                 v   => this%v,               &
                 Dq  => this%Dq,              &
                 Dq_old  => this%Dq_old,      &
                 q   => this%q,               &
                 t   => this%t,               &
                 l   => this%l,               &
                 eta => this%eta               )
         ! calculation of step size $h$
         h = t1 - t

         ! calculate B(q_n)
         Bq = this%SHAKELie_B(q)
         ! count calls
         this%SHAKELie_stats%nBcalls = this%SHAKELie_stats%nBcalls + 1

         ! initialization of the new values
         !l = l
         eta = 0 ! So Dqs is consistent without this:
         Dq  = Dq_old
         Dqs = Dq ! + matmul(Bq, eta)


         ! loop of the newton method
         converged = .false.
         do i=1,this%opts%imax

            ! calculate the new value $q_{k}$
            q1  = this%SHAKELie_qlpexphDqtilde(q,h,Dqs)

            ! Calculate M_Dq
            if (this%opts%const_mass_matrix == 1) then
               if (this%opts%diag_mass_matrix == 1) then
                  M_Dq = this%SHAKELie_const_diag_M*Dq
               else
                  M_Dq = matmul(this%SHAKELie_const_M,Dq)
               end if
            else
               if (this%opts%diag_mass_matrix == 1) then
                  M_Dq = this%SHAKELie_diag_M(q)*Dq
               else
                  M_Dq = matmul(this%SHAKELie_M(q),Dq)
               end if
            end if

            ! Calculate M_Dq_old
            if (this%opts%const_mass_matrix == 1) then
               if (this%opts%diag_mass_matrix == 1) then
                  M_Dq_old = this%SHAKELie_const_diag_M*Dq_old
               else
                  M_Dq_old = matmul(this%SHAKELie_const_M,Dq_old)
               end if
            else
               if (this%opts%diag_mass_matrix == 1) then
                  M_Dq_old = this%SHAKELie_diag_M(q)*Dq_old
               else
                  M_Dq_old = matmul(this%SHAKELie_M(q),Dq_old)
               end if
            end if

            ! calculate the residue $res_h$
            res(1:sv) = matmul(this%SHAKELie_Tg_inv_T(h*Dq_old),M_Dq_old) &
               - matmul(this%SHAKELie_Tg_inv_T(-h*Dq),M_Dq) - h*this%SHAKELie_f(q,v,t) - matmul(transpose(Bq),h*l)
            res(sv+1:sv+sl) = this%SHAKELie_phi(q1) / h
            res(sv+sl+1:sv+2*sl) = matmul(Bq, (Dq_old + Dq)/2)
            ! count calls
            this%SHAKELie_stats%ngcalls = this%SHAKELie_stats%ngcalls + 1

            ! solve the linear System $S_t \cdot \Delta xl = -res$
            Dxi = -res ! the result will be saved in the right hand side's spot

            if (i == 1 .or. this%opts%recalc_iteration_matrix == 1) then
#ifdef STNUM
               St = this%St_num(Dq,h,l,eta,t1)
#else
               ! calculate iteration matrix
               St(1:sv,1:sv) = -this%SHAKELie_Tg_inv_T(-h*Dq)
               if (this%opts%const_mass_matrix == 1) then
                  if (this%opts%diag_mass_matrix == 1) then
                     forall (i=1:sv)
                        St(:,i) = St(:,i)*this%SHAKELie_const_diag_M(i)
                     end forall
                  else
                     St(1:sv,1:sv) = matmul(St(1:sv,1:sv), this%SHAKELie_const_M)
                  end if
               else
                  ERROR STOP "Nonconstant mass matrix not supported yet"
               end if
               St(1:sv,1:sv) = St(1:sv,1:sv) + h*this%SHAKELie_d_Tg_inv_T(-h*Dq, M_Dq)

               St(      1:sv      ,    sv+1:sv+sl   ) = -transpose(Bq)
               St(sv   +1:sv+sl   ,       1:sv      ) = matmul(this%SHAKELie_B(q1), this%SHAKELie_Tg(h,Dqs))
               St(sv   +1:sv+sl   ,    sv+1:sv+sl   ) = 0.0_8

               St(      1:sv+sl   , sv+sl+1:sv+sl+sl) = 0.0_8
               St(sv+sl+1:sv+sl+sl, sv+sl+1:sv+sl+sl) = matmul(Bq,transpose(Bq))/2
               St(sv+sl+1:sv+sl+sl,       1:sv      ) = Bq/2
               St(sv+sl+1:sv+sl+sl, sv+   1:sv+sl   ) = 0.0_8

               ! count calls
               this%SHAKELie_stats%nBcalls = this%SHAKELie_stats%nBcalls + 1

! endif STNUM
#endif

               ! Actually solve the linear System $S_t \cdot \Delta xi = -res$
               if (this%opts%banded_iteration_matrix == 1) then
#ifdef DEBUG_PRINT_ITERATION_MATRIX_AT
                  if (t >= DEBUG_PRINT_ITERATION_MATRIX_AT) then
                     call print_matrix(St,'St')
                     call print_vector_int(this%opts%jour,'jour')
                     stop "Successfully printed iteration matrix"
                  end if
#endif
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag)
                     if (.not. allocated(Stb)) then
                        allocate(Stb(lo*2+hi+1, sv+2*sl))
                     end if
                     !Stb = 0.0_8 ! DEBUG
                     forall (i=1:sv+2*sl)
                        Stb((max(1, i-hi)-i+lo+1+hi):(min(sv+2*sl,i+lo)-i+lo+1+hi),i) &
                            = St(this%opts%jour(max(1, i-hi):min(sv+2*sl, i + lo)), this%opts%jour(i))
                     end forall
                     Dxi = Dxi(this%opts%jour)
                     call dgbsv(           &! solve the system A*X=B and save the result in B
                                 sv+2*sl,  &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,       &!
                                 hi,       &!
                                 1,        &! number of right hand sides (=size(B,2))
                                 Stb,      &! matrix A
                                 2*lo+hi+1,&! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv,     &! integer pivot vector; it is not needed
                                 Dxi,      &! matrix B
                                 sv+2*sl,  &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)      ! integer information flag
                     Dxi(this%opts%jour) = Dxi
                     !call print_vector(Dxl,'Dxl_b') ! DEBUG
                     !stop ! DEBUG
                  end associate
               else
                  call dgesv(          &! solve the System A*X=B and save the result in B
                              sv+2*sl, &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                              1,       &! number of right hand sides (=size(B,2))
                              St,      &! matrix A
                              sv+2*sl, &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv,    &! integer pivot vector; it is not needed
                              Dxi,     &! matrix B
                              sv+2*sl, &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                              info)     ! integer information flag
                  ! TODO: Errors abfragen
               end if

               if (info.ne.0) print *, "dgesv sagt info=", info

            else ! (St was not (re)calculated)
               ! Actually solve the linear System $S_t \cdot \Delta xl = -res$
               if (this%opts%banded_iteration_matrix == 1) then
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag, &
                             uStb => ubound(Stb,1))
                     Dxi = Dxi(this%opts%jour) !TODO
                     call dgbtrs(               &
                                 'No transpose',&! solve the System A*X=B and save the result in B
                                 sv+2*sl,       &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,            &!
                                 hi,            &!
                                 1,             &! number of right hand sides (=size(B,2))
                                 Stb,           &! matrix A
                                 2*lo+hi+1,     &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv,          &! integer pivot vector; it is not needed
                                 Dxi,           &! matrix B
                                 sv+2*sl,       &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)           ! integer information flag
                     Dxi(this%opts%jour) = Dxi
                  end associate
               else
                  call dgetrs(                &
                              'No transpose', & ! solve the System A*X=B and save the result in B
                              sv + 2*sl,      & ! number of linear equations (=size(A,1))  ! Vorsicht: double precision muss real(8) sein
                              1,              & ! number of right hand sides (=size(B,2))
                              St,             & ! matrix A
                              sv + 2*sl,      & ! leading dimesnion of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv,           & ! integer pivot vector; holds the pivoting that was calculated in the dgesv call
                              Dxi,            & ! matrix B
                              sv + 2*sl,      & ! leading dimension of B, in thie case is equal to the number of linear equations (=size(A,1))
                              info)             ! integer information flag
               end if
               ! DEBUG
               if (info.ne.0) print*, 'dgetrs sagt info=', info
               ! GUBED
            end if

            ! Check for convergence
            if ( ( sum((Dxi(      1:      sv) / (this%opts%atol + this%opts%rtol*abs(Dqs)))**2 )   &
                  +sum((Dxi(sv+   1:   sv+sl) / (this%opts%atol + this%opts%rtol*abs(h*l)))**2 )   &
                  +sum((Dxi(sv+sl+1:sv+sl+sl) / (this%opts%atol + this%opts%rtol*abs(eta)))**2 ) ) &
                   / (sv+sl+sl)   <=  1.0_8 ) then
               converged = .true.
            end if

            ! update
! DEBUG
#define REL *1
            Dqs   = Dqs  + Dxi(      1:sv      )     REL
            l     = l    + Dxi(sv   +1:sv+sl   )/h   REL
            eta   = eta  + Dxi(sv+sl+1:sv+sl+sl)     REL
#undef REL
            Dq = Dqs + matmul(transpose(Bq),eta)

            ! update solver stats
            this%SHAKELie_stats%newt_steps_curr = i

            ! exit loop, if converged
            if ( converged ) exit

         end do ! end newton method

         ! DEBUG
         if (this%SHAKELie_stats%newt_steps_curr >= this%opts%imax) then
            print *, ' ' ! DEBUG
            print *, "Newton did not converge at t=", t
            !! DEBUG !
            call this%SHAKELie_print_stats
            print *, "Iteration matrix:"
            !call print_matrix(St,'St')
            error stop
            !! GUBED !
         end if

         ! The calculation of the new q is done
         t      = t1
         q      = this%SHAKELie_qlpexphDqtilde(q, h, Dqs)
         v      = (Dq_old + Dq)/2
         Dq_old = Dq
      end associate
   end subroutine SHAKELie_solveConstrainedTimeStep_stab2

   ! subroutine for integrating one time step
   subroutine SHAKELie_solveTimeStep(this,t1)
      implicit none
      class(SHAKELie_problem), intent(inout)      :: this  ! problem object
      real(8),               intent(in   )      :: t1    ! $t_{n+1}$

      ! internal integer variables
      integer                                   :: i     ! for iteration
      integer, dimension(this%sizev)            :: ipiv  ! needed for dgesv
      integer                                   :: info  ! needed for dgesv

      ! internal real variables
      real(8)                                   :: h       ! step size
      real(8), dimension(this%sizev)            :: res     ! $res$
      real(8), dimension(this%sizev,this%sizev) :: St      ! $S_t$
      real(8), dimension(:,:), allocatable      :: Stb     ! $S_t$ in banded format
      real(8), dimension(this%sizev)            :: Dq      ! $v_{1/2}
      real(8), dimension(this%sizev)            :: Dx      ! $\Delta x$
      real(8), dimension(this%sizev)            :: M_Dq

      ! internal logical
      logical                                   :: converged

      ERROR STOP "Nonconstrained not supported yet"
#if 0
      ! associate construct for better readability
      associate (p      => this%p,               &
                 sv     => this%sizev,           &
                 v      => this%v,               &
                 q      => this%q,               &
                 t      => this%t                )
         ! calculation of step size $h$
         h = t1 - t

         ! Initial guess
         ! DEBUG
         Dq = 0.0_8
         !Dq = v

         ! loop of the newton method
         converged = .false.
         do i=1,this%opts%imax
            ! caluclate the residue $res$
            if (this%opts%const_mass_matrix == 1) then
               if (this%opts%diag_mass_matrix == 1) then
                  M_Dq = this%SHAKELie_const_diag_M*Dq
               else
                  M_Dq = matmul(this%SHAKELie_const_M,Dq)
               end if
            else
               if (this%opts%diag_mass_matrix == 1) then
                  M_Dq = this%SHAKELie_diag_M(q)*Dq
               else
                  M_Dq = matmul(this%SHAKELie_M(q),Dq)
               end if
            end if
            ! Apply inverse tangent operator and add f
            res = -p  + matmul(this%SHAKELie_Tg_inv_T(-h*Dq),M_Dq) + h/2*this%SHAKELie_f(q,v,t)
            ! count calls
            this%SHAKELie_stats%ngcalls = this%SHAKELie_stats%ngcalls + 1


            Dx = -res ! the result will be saved in the right hand side's spot

            if (i == 1  .or. this%opts%recalc_iteration_matrix == 1) then
               ! calculate iteration matrix
               St = this%SHAKELie_Tg_inv_T(-h*Dq)
               if (this%opts%const_mass_matrix == 1) then
                  if (this%opts%diag_mass_matrix == 1) then
                     forall (i=1:sv)
                        St(:,i) = St(:,i)*this%SHAKELie_const_diag_M(i)
                     end forall
                  else
                     St = matmul(St, this%SHAKELie_const_M)
                  end if
               else
                  ERROR STOP "Nonconstant mass matrix not supported yet"
               end if
               St = St - h*this%SHAKELie_d_Tg_inv_T(-h*Dq, M_Dq)

               !! DEBUG !
               !print *, 'i=', i
               !print *, 't=', t
               !call print_vector(q,'q')
               !call print_vector(v,'v')
               !call print_vector(vd,'vd')
               !call print_vector(a,'a')
               !!call print_matrix(St,'Stn')
               !call print_vector(res, 'res')
               !!call print_vector(this%SHAKELie_const_diag_M,'Mn')
               !!print *, "t=", t
               !!if (t > 2.0_8) then
               !!   call print_matrix(St,'St')
               !!   print *, "cond=", mycond(St)
               !!   exit
               !!endif
               !!print *,
               !!stop
               !! GUBED !

               ! solve the linear System $S_t \cdot \Delta x = -res$
               if (this%opts%banded_iteration_matrix == 1) then
#ifdef DEBUG_PRINT_ITERATION_MATRIX_AT
                  if (t >= DEBUG_PRINT_ITERATION_MATRIX_AT) then
                     call print_matrix(St,'St')
                     call print_vector_int(this%opts%jour,'jour')
                     stop "Successfully printed iteration matrix"
                  end if
#endif
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag, &
                             uStb => ubound(Stb,1))
                     if (.not. allocated(Stb)) then
                        allocate(Stb(lo*2+hi+1,this%sizev))
                     end if
                     forall (i=1:this%sizev)
                        Stb((max(1, i-hi)-i+lo+1+hi):(min(this%sizev,i+lo)-i+lo+1+hi),i) &
                            = St(max(1, i-hi):min(this%sizev, i + lo),i)
                     end forall
                     ! TODO: Conditional jump or move depends on unititialised value(s) in the next line?
                     call dgbsv(       &! solve the system A*X=B and save the result in B
                                 sv,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,   &!
                                 hi,   &!
                                 1,    &! number of right hand sides (=size(B,2))
                                 Stb,  &! matrix A
                                 2*lo+hi+1, &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv, &! integer pivot vector; it is not needed
                                 Dx,   &! matrix B
                                 sv,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)  ! integer information flag
                  end associate
               else
                  call dgesv(       &! solve the System A*X=B and save the result in B
                              sv,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                              1,    &! number of right hand sides (=size(B,2))
                              St,   &! matrix A
                              sv,   &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv, &! integer pivot vector; it is not needed
                              Dx,   &! matrix B
                              sv,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                              info)  ! integer information flag
               end if
               ! TODO: Errors abfragen
               !!! DEBUG
               !   print *, 'St = ...'
               !   print *, '[', St(1,1), ',', St(1,2), ';'
               !   print *, ' ', St(2,1), ',', St(2,2), '];'
               !   print *, 'mres = ...'
               !   print *, '[', -res(1), ',', -res(2), '];'
               !   print *, 'Dx = ...'
               !   print *, '[', Dx(1), ',', Dx(2), '];'
               !   print *, ''
               if (info.ne.0) print *, 'dgesv sagt info=', info
               !! GEBUD
            else
               ! solve the linear System $S_t \cdot \Delta x = -res$
               if (this%opts%banded_iteration_matrix == 1) then
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag  )
                     call dgbtrs(      &
                                 'No transpose', &! solve the System A*X=B and save the result in B
                                 sv,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,   &!
                                 hi,   &!
                                 1,    &! number of right hand sides (=size(B,2))
                                 Stb,  &! matrix A
                                 2*lo+hi+1, &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv, &! integer pivot vector; it is not needed
                                 Dx,   &! matrix B
                                 sv,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)  ! integer information flag
                  end associate
               else
                  call dgetrs(      &
                              'No transpose', &! solve the System A*X=B and save the result in B
                              sv,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                              1,    &! number of right hand sides (=size(B,2))
                              St,   &! matrix A
                              sv,   &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv, &!
                              Dx,   &! matrix B
                              sv,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                              info)  ! integer information flag
               end if
               ! DEBUG
               if (info.ne.0) print *, 'dgesv sagt info=', info
               ! GEBUD
            end if




            ! Check for convergence
            if ( sum((Dx / (this%opts%atol + this%opts%rtol*abs(Dq)))**2 )  &
                  / (sv)   <=  1.0_8 ) then
               converged = .true.
            end if

! DEBUG!!
#define REL *1
            ! update $\dot v_{n+1}$, $v_{n+1}$ and $\Delta q_n$
            Dq  = Dq  + Dx                       REL
#undef REL
! GUBED!!

            ! update solver stats
            this%SHAKELie_stats%newt_steps_curr = i

            ! exit loop, if converged
            if ( converged ) exit

         end do ! end newton method

         !! DEBUG
         if (this%SHAKELie_stats%newt_steps_curr >= this%opts%imax) then
            print *, ' ' ! DEBUG
            print *, "Newton did not converge at t=", t
            !!! DEBUG !
            !call this%SHAKELie_print_stats
            !print *, "Iteration matrix:"
            !call print_matrix(St,'St')
            errorstop
            !!! GUBED !
         end if

         ! Calculate new t, q, p and v
         t = t1
         q = this%SHAKELie_qlpexphDqtilde(q,h,Dq)
         if (this%opts%const_mass_matrix == 1) then
            if (this%opts%diag_mass_matrix == 1) then
               M_Dq = this%SHAKELie_const_diag_M*Dq
            else
               M_Dq = matmul(this%SHAKELie_const_M,Dq)
            end if
         else
            if (this%opts%diag_mass_matrix == 1) then
               M_Dq = this%SHAKELie_diag_M(q)*Dq
            else
               M_Dq = matmul(this%SHAKELie_M(q),Dq)
            end if
         end if
         p = matmul(this%SHAKELie_Tg_inv_T(h*Dq),M_Dq) - h/2*this%SHAKELie_f(q,Dq,t) ! TODO: Hier müsste ein v stehen, kein Dq, ist aber egal, wenn das Modell ungedämpft ist
         ! count calls
         this%SHAKELie_stats%ngcalls = this%SHAKELie_stats%ngcalls + 1

         ! Calculate velocity from impulse
         if (this%opts%const_mass_matrix == 1) then
            if (this%opts%diag_mass_matrix == 1) then
               forall (i=1:sv)
                  v(i) = p(i)/this%SHAKELie_const_diag_M(i)
               end forall
            else
               ERROR STOP "Nondiagonal mass matrix not supported yet"
            end if
         else
            ERROR STOP "Nonconstant mass matrix not supported yet"
         end if

      end associate
#endif

   end subroutine SHAKELie_solveTimeStep


!   pure function cum_sum_omit_last(vec) result(rslt)
!      ! input
!      real(8), intent(in)  :: vec(:)
!      ! result
!      real(8)              :: rslt(size(vec)-1)
!      ! internal
!      integer              :: i
!      !
!      rslt(1) = vec(1)
!      do i=2,size(vec)-1
!         rslt(i) = rslt(i-1) + vec(i)
!      end do
!   end function cum_sum_omit_last


   ! subroutine to integrate
   subroutine SHAKELie_integrate(this)
      implicit none
      ! input/output
      class(SHAKELie_problem),  intent(inout)  :: this
      ! internal
      class(SHAKELie_problem), allocatable     :: oldstep
      integer  :: n  ! needed for iteration
      real(8)  :: h  ! step size $h$
      real(8)  :: t1 ! next time $t_{n+1}$
      real(4)  :: times(2), time ! for dtime

      ! problem initialization
      call this%SHAKELie_init()

      ! TODO: More error checking
      if (this%opts%constrained == 1             .and. &
          this%opts%banded_iteration_matrix == 1 .and. &
          .not. allocated(this%opts%jour)               ) then
         print *, '[ERROR] opts%jour not allocated, aborting'
         errorstop
      end if
      ! TODO: t0 und t vergleichen (müssen gleich sein)

      ! initialize output function
      call this%SHAKELie_outputFunction(0)

      ! Calculate step size $h$
      h = (this%opts%te - this%opts%t0)/this%opts%nsteps

      ! Set stats of solver to zero
      this%SHAKELie_stats%newt_steps_curr = 0
      this%SHAKELie_stats%newt_steps_sum = 0
      this%SHAKELie_stats%newt_steps_max = 0
      this%SHAKELie_stats%newt_steps_avg = 0
      this%SHAKELie_stats%ngcalls = 0
      this%SHAKELie_stats%nBcalls = 0

      ! if mass matrix is constant, calculate it
      if (this%opts%const_mass_matrix == 1) then
         if (this%opts%diag_mass_matrix == 1) then
            if (.not. allocated(this%SHAKELie_const_diag_M)) then
               allocate(this%SHAKELie_const_diag_M(this%sizev))
            end if
            this%SHAKELie_const_diag_M = this%SHAKELie_diag_M(this%q)
         else
            if (.not. allocated(this%SHAKELie_const_M)) then
               allocate(this%SHAKELie_const_M(this%sizev,this%sizev))
            end if
            this%SHAKELie_const_M = this%SHAKELie_M(this%q)
         end if
      end if

      ! start stopwatch
      time = dtime(times)

      ! allocate space for Dq and Dq_old
      allocate(this%Dq(this%sizev))
      allocate(this%Dq_old(this%sizev))

      if (this%opts%constrained == 0) then

         ! calculate initial values
         call this%SHAKELie_calcInitial(h)

         ! output for the first time
         call this%SHAKELie_outputFunction(1)

         ! integration loop
         do n=1,this%opts%nsteps
            ! Calculate the next time $t_{n+1}$
            t1 = this%opts%t0 + n*h
            ! solve time step
            call this%SHAKELie_solveTimeStep(t1)
            ! update solver stats
            this%SHAKELie_stats%newt_steps_sum = this%SHAKELie_stats%newt_steps_sum &
               + this%SHAKELie_stats%newt_steps_curr
            this%SHAKELie_stats%newt_steps_avg = real(this%SHAKELie_stats%newt_steps_sum,8)/n
            this%SHAKELie_stats%newt_steps_max = max( &
               this%SHAKELie_stats%newt_steps_max,    &
               this%SHAKELie_stats%newt_steps_curr  )
            ! output normally
            call this%SHAKELie_outputFunction(1)
         end do


      else

         ! Allocate object to remember last step
         allocate(oldstep, source=this)

         ! First step has to be done separately
         !TODO DEBUG
         this%Dq_old = this%v
         this%l = 0.0_8
         ! GUBED

         ! integration loop, one extra step to calculate the last v and lambda
         do n=1,this%opts%nsteps+1
            ! Remember old values
            oldstep%t = this%t
            oldstep%q = this%q

            ! Calculate the next time $t_{n+1}$
            t1 = this%opts%t0 + n*h

            if (this%opts%stab2 == 1) then
               ! solve time step with the stabilized index-2 system
               call this%SHAKELie_solveConstrainedTimeStep_stab2(t1)
            else
               ERROR STOP "SHAKELie does not support index-3 integration"
               ! solve time step
               call this%SHAKELie_solveConstrainedTimeStep(t1)
            end if

            ! update solver stats
            this%SHAKELie_stats%newt_steps_avg =           &
               (this%SHAKELie_stats%newt_steps_avg*(n-1) + &
                this%SHAKELie_stats%newt_steps_curr)/n
            this%SHAKELie_stats%newt_steps_max = max( &
               this%SHAKELie_stats%newt_steps_max,    &
               this%SHAKELie_stats%newt_steps_curr  )

            !

            ! Set correct Lagrange multiplier and velocity of old step and output
            oldstep%v = this%v
            oldstep%l = this%l
            call oldstep%SHAKELie_outputFunction(1)
         end do

         ! Deallocate variable oldstep
         deallocate(oldstep)

      end if

      ! stop stopwatch
      this%SHAKELie_stats%time = dtime(times) - time

      ! output to terminate
      call this%SHAKELie_outputFunction(99)

   end subroutine SHAKELie_integrate

   pure function SHAKELie_num_Ct(this, q, v, t) result(rslt)
      ! input
      class(SHAKELie_problem),   intent(in)         :: this
      real(8), dimension(:), intent(in)         :: q
      real(8), dimension(:), intent(in)         :: v
      real(8),               intent(in)         :: t
      ! result
      real(8), dimension(this%sizev,this%sizev) :: rslt
      ! internal
      integer                                   :: i
      integer                                   :: j
      integer                                   :: diags
      real(8), dimension(this%sizev)            :: w
      real(8), dimension(this%sizev)            :: g0
      !
      ERROR STOP 'No Ct needed'
#if 0
      !
      ! Ct is the derivative of r wrt v. But v only appears in g, so
      ! it is sufficient to calculate the derivative of g.
      !
      ! Calculate g ! TODO: man könnte g0 durchschleifen oder so
      g0 = this%SHAKELie_g(q,v,t)
#ifdef NOBANDEDDER
      if (.false.) then
#else
      ! if the iteration matrix is banded, Ct must be banded, as well
      if (this%opts%banded_iteration_matrix == 1) then
#endif
         ! the number of total off-diagonals
         diags = this%opts%nr_subdiag+this%opts%nr_superdiag
         ! most of the result will be zero
         rslt = 0.0_8
         ! loop over the first diags columns of Ct
         do i=1,min(diags+1,this%sizev)
            ! Construct vector w
            w = 0.0_8
            do j=i, this%sizev, 1+diags
               w(j) = max( abs(v(j))*1.0e-8_8,  1.0e-12_8)
            end do
            ! calculate finite differences
            rslt(:,i) = this%SHAKELie_g(q,v+w,t) - g0
            ! move the parts of the finite difference, that belog to
            ! an other column of Ct, also divide by h=w(k)
            do j=i, this%sizev, 1+diags
               rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), j) =      &
                     rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), i)  &
                     / w(j)
               if (i .ne. j) then
                  rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), i) = 0.0_8
               end if
            end do
         end do
      else
         ! set w to zero
         w = 0.0_8
         ! loop over the columns of Ct
         do i=1,this%sizev
            if (.not. i == 1) w(i-1) = 0.0_8
            w(i) = max( abs(v(i))*1.0e-8_8,  1.0e-12_8)
            rslt(:,i) = (this%SHAKELie_g(q,v+w,t) - g0)/w(i)
         end do
      end if
#endif
   end function SHAKELie_num_Ct

   pure function SHAKELie_num_Kt(this, q, v, vd, t) result(rslt)
      ! input
      class(SHAKELie_problem),   intent(in)         :: this
      real(8), dimension(:), intent(in)         :: q
      real(8), dimension(:), intent(in)         :: v
      real(8), dimension(:), intent(in)         :: vd
      real(8),               intent(in)         :: t
      ! result
      real(8), dimension(this%sizev,this%sizev) :: rslt
      ! internal
      integer                                   :: i
      integer                                   :: j
      integer                                   :: diags
      real(8), dimension(this%sizev)            :: w
      real(8), dimension(this%sizev)            :: g0
      !
      ERROR STOP 'No Kt needed'
#if 0
      !
      ! Kt is the derivative of r wrt q. WE ASSUME THAT THE MASS MATRIX
      ! M DOES NOT DEPEND ON Q, OTHERWISE THIS APPROXIMATION WILL BE
      ! FALSE! TODO! TODO!
      ! Thus, it is sufficient to calculate the derivative of g.
      !
      ! Calculate g ! TODO: man könnte g0 durchschleifen oder so
      g0 = this%SHAKELie_g(q,v,t)
#ifdef NOBANDEDDER
      if (.false.) then
#else
      ! if the iteration matrix is banded, Kt must be banded, as well
      if (this%opts%banded_iteration_matrix == 1) then
#endif
         ! most of Kt is zero
         rslt = 0.0_8
         ! the number of total off-diagonals
         diags = this%opts%nr_subdiag+this%opts%nr_superdiag
         ! loop over the first diags columns of Kt
         do i=1,min(diags+1,this%sizev)
            ! Construct vector w
            w = 0.0_8
            do j=i, this%sizev, 1+diags
               w(j) = max( abs(q(j))*1.0e-8_8,  1.0e-12_8)
            end do
            ! calculate finite differences
            rslt(:,i) = this%SHAKELie_g(this%SHAKELie_qlpexphDqtilde(q,1.0_8, w),v,t) - g0
            ! move the parts of the finite difference, that belong to
            ! an other column of Kt, also divide by h=w(k)
            do j=i, this%sizev, 1+diags
               rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), j) =     &
                     rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), i) &
                       / w(j)
               if (i .ne. j) then
                  rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), i) = 0.0_8
               end if
            end do
         end do
      else
         ! set w to zero
         w = 0.0_8
         ! loop over the columns of Ct
         do i=1,this%sizev
            if (.not. i == 1) w(i-1) = 0.0_8
            w(i) = max( abs(q(i))*1.0e-8_8,  1.0e-12_8)
            rslt(:,i) = (this%SHAKELie_g(this%SHAKELie_qlpexphDqtilde(q,1.0_8, w),v,t) - g0) / w(i)
            !w    = 0.0_8
            !w(i) = h
            !rslt(:,i) = (this%SHAKELie_g(q+this%SHAKELie_tilde(w),v,t) - g0)/h
         end do
      end if
#endif
   end function SHAKELie_num_Kt

! DEBUG
   pure function SHAKELie_num_B(this, q) result(rslt)
      ! input
      class(SHAKELie_problem),   intent(in)         :: this
      real(8), dimension(:), intent(in)         :: q
      ! result
      real(8), dimension(this%sizel,this%sizev) :: rslt
      ! internal
      integer                                   :: i
      real(8), dimension(this%sizev)            :: w
      real(8), dimension(this%sizev)            :: Phi0
      !
      !
      ! Calculate Phi
      Phi0 = this%SHAKELie_Phi(q)

      ! set w to zero
      w = 0.0_8
      ! loop over the columns of Ct
      do i=1,this%sizev
         if (.not. i == 1) w(i-1) = 0.0_8
         w(i) = max( abs(q(i))*1.0e-8_8,  1.0e-12_8)
         rslt(:,i) = (this%SHAKELie_Phi(this%SHAKELie_qlpexphDqtilde(q,1.0_8, w)) - Phi0) / w(i)
         !w    = 0.0_8
         !w(i) = h
         !rslt(:,i) = (this%SHAKELie_g(q+this%SHAKELie_tilde(w),v,t) - g0)/h
      end do
   end function SHAKELie_num_B
! GUBED

   pure function SHAKELie_num_Kt_lambda(this, q, v, vd, l, t) result(rslt)
      ! input
      class(SHAKELie_problem),            intent(in)   :: this
      real(8), dimension(:),          intent(in)   :: q
      real(8), dimension(:),          intent(in)   :: v
      real(8), dimension(:),          intent(in)   :: vd
      real(8), dimension(:),          intent(in)   :: l
      real(8),                        intent(in)   :: t
      ! result
      real(8), dimension(this%sizev,this%sizev)    :: rslt
      ! internal
      integer                                      :: i
      integer                                      :: j
      integer                                      :: diags
      real(8), dimension(this%sizev)               :: w
      real(8), dimension(this%sizev)               :: g0
      !
      ERROR stop "Kt not needed"
#if 0
      ! Kt is the derivative of r wrt q. WE ASSUME THAT THE MASS MATRIX
      ! M DOES NOT DEPEND ON Q, OTHERWISE THIS APPROXIMATION WILL BE
      ! FALSE! TODO! TODO!
      ! Thus, it is sufficient to calculate the derivative of g.
      !
      ! Calculate g ! TODO: man könnte g0 durchschleifen oder so
      g0 = this%SHAKELie_g(q,v,t) + matmul(transpose(this%SHAKELie_B(q)),l)
#ifdef NOBANDEDDER
      if (.false.) then
#else
      ! if the iteration matrix is banded, Kt must be banded, as well
      if (this%opts%banded_iteration_matrix == 1) then
#endif
         ! most of Kt is zero
         rslt = 0.0_8
         ! the number of total off-diagonals
         diags = this%opts%nr_subdiag+this%opts%nr_superdiag
         ! loop over the first diags columns of Kt
         do i=1,min(diags+1,this%sizev)
            ! Construct vector w
            w = 0.0_8
            do j=i, this%sizev, 1+diags
               w(j) = max( abs(q(j))*1.0e-8_8,  1.0e-12_8)
            end do
            ! calculate finite differences
            rslt(:,i) = this%SHAKELie_g(this%SHAKELie_qlpexphDqtilde(q,1.0_8, w),v,t) &
               + matmul(transpose(this%SHAKELie_B(this%SHAKELie_qlpexphDqtilde(q,1.0_8,w))),l)- g0
            ! move the parts of the finite difference, that belog to
            ! an other column of Kt, also divide by h=w(k)
            do j=i, this%sizev, 1+diags
               rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), j) =     &
                     rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), i) &
                     / w(j)
               if (i .ne. j) then
                  rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), i) = 0.0_8
               end if
            end do
         end do
      else
         ! set w to zero
         w = 0.0_8
         ! loop over the columns of Ct
         do i=1,this%sizev
            if (.not. i == 1)  w(i-1) = 0.0_8
            w(i) = max( abs(q(i))*1.0e-8_8,  1.0e-12_8)
            !rslt(:,i) = (this%SHAKELie_g(q+this%SHAKELie_tilde(w),v,t) + matmul(transpose(this%SHAKELie_B(q+this%SHAKELie_tilde(w))),l) - g0)/h
            !w(i) = 1.0_8 ! TODO
            rslt(:,i) = (this%SHAKELie_g(this%SHAKELie_qlpexphDqtilde(q,1.0_8,w),v,t)  &
             + matmul(transpose(this%SHAKELie_B(this%SHAKELie_qlpexphDqtilde(q,1.0_8,w))),l) - g0)/w(i)
         end do
      end if
#endif
   end function SHAKELie_num_Kt_lambda

   subroutine SHAKELie_print_stats(this)
      ! input
      class(SHAKELie_problem), intent(in)  :: this
      !
      print *, 'time:          ', this%SHAKELie_stats%time
      print *, '#calls of g:   ', this%SHAKELie_stats%ngcalls
      print *, '#calls of B:   ', this%SHAKELie_stats%nBcalls
      print *, 'newt_steps_max:', this%SHAKELie_stats%newt_steps_max
      print *, 'newt_steps_avg:', this%SHAKELie_stats%newt_steps_avg
   end subroutine

   subroutine print_matrix(A,Aname)
      implicit none
      ! input
      real(8), dimension(:,:) :: A
      character(len=*)       :: Aname
      ! internal
      integer                 :: i
      !
      print *, '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print *, trim(Aname), ' = ['
      do i=1,ubound(A,1) - 1
         print *, A(i,:), ';'
      end do
      print *, A(ubound(A,1),:), '];'
   end subroutine print_matrix

   subroutine print_vector_int(A,Aname)
      implicit none
      ! input
      integer, dimension(:)   :: A
      character(len=*)        :: Aname
      !
      print *, '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print *, trim(Aname), ' = ['
      print *, A, '];'
   end subroutine print_vector_int

   subroutine print_vector(A,Aname)
      implicit none
      ! input
      real(8), dimension(:)   :: A
      character(len=*)        :: Aname
      !
      print *, '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print *, trim(Aname), ' = ['
      print *, A, '];'
   end subroutine print_vector

   subroutine SHAKELie_cleanup(this)
      implicit none
      ! input/output
      class(SHAKELie_problem), intent(inout)  :: this
      !
      this%opts%constrained = 0
      this%opts%const_mass_matrix = 0
      this%opts%diag_mass_matrix = 0
      this%opts%banded_iteration_matrix = 0
      this%opts%nr_subdiag = 0
      this%opts%nr_superdiag = 0
      this%opts%recalc_iteration_matrix = 0
      this%opts%use_num_Ct = 1
      this%opts%use_num_Kt = 1
      this%opts%atol = 1.0e-10_8
      this%opts%rtol = 1.0e-8_8
      this%opts%imax   = 5
      this%opts%t0 = 0.0_8
      this%opts%te = 1.0_8
      this%opts%nsteps = 100
      !
      this%SHAKELie_stats%newt_steps_curr = 0
      this%SHAKELie_stats%newt_steps_max = 0
      this%SHAKELie_stats%newt_steps_avg = 0.0_8
      this%SHAKELie_stats%ngcalls = 0
      this%SHAKELie_stats%nBcalls = 0
      this%SHAKELie_stats%time = 0.0_8
      !
      this%t = 0.0_8
      this%sizeq = 0
      this%sizev = 0
      this%sizel = 0
      !
      if (allocated(this%q))  deallocate(this%q)
      if (allocated(this%v))  deallocate(this%v)
      if (allocated(this%l))  deallocate(this%l)
      if (allocated(this%eta))  deallocate(this%eta)
      if (allocated(this%Dq))  deallocate(this%Dq)
      if (allocated(this%Dq_old))  deallocate(this%Dq_old)
      if (allocated(this%SHAKELie_const_M)) deallocate(this%SHAKELie_const_M)
      if (allocated(this%SHAKELie_const_diag_M)) deallocate(this%SHAKELie_const_diag_M)
      if (allocated(this%opts%jour)) deallocate(this%opts%jour)
   end subroutine SHAKELie_cleanup

#if 1
   function mycond(A) result(rslt)
      implicit none
      ! input
      real(8), intent(in)  :: A(:,:)
      ! result
      real(8)              :: rslt
      ! internal
      real(8)              :: AA(size(A,1),size(A,2))
      real(8)              :: S(min(size(A,1), size(A,2)))
      real(8), allocatable :: work(:)
      integer              :: lwork
      real(8)              :: U(size(A,1),size(A,1))
      real(8)              :: VT(size(A,2),size(A,2))
      integer              :: iwork(8*min(size(A,1), size(A,2)))
      integer              :: info
      !
      AA = A
      lwork = 3*min(size(A,1),size(A,2)) + max(max(size(A,1),size(A,2)),7*min(size(A,1),size(A,2)))
      allocate(work(lwork))
      call DGESDD('A',       &!JOBZ, calculate singular vectors TODO: Will fail if they are not calculated
                  size(A,1), &!M,
                  size(A,2), &!N,
                  AA,        &!A,
                  size(A,1), &!LDA,
                  S,         &!S,
                  U,         &!U,
                  size(A,1), &!LDU,
                  VT,        &!VT,
                  size(A,2), &!LDVT,
                  work,      &!WORK,
                  -1,        &!LWORK,
                  iwork,     &!IWORK
                  info)      !INFO )
      if (info /= 0) then
         errorstop "In cond: dgesdd did not succeed"
      end if
      AA = A
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
      call DGESDD('A',       &!JOBZ, calculate singular vectors TODO: Will fail if they are not calculated
                  size(A,1), &!M,
                  size(A,2), &!N,
                  AA,        &!A,
                  size(A,1), &!LDA,
                  S,         &!S,
                  U,         &!U,
                  size(A,1), &!LDU,
                  VT,        &!VT,
                  size(A,2), &!LDVT,
                  work,      &!WORK,
                  lwork,     &!LWORK,
                  iwork,     &!IWORK
                  info)      !INFO )
      if (info /= 0) then
         errorstop "In cond: dgesdd did not succeed"
      end if
      if (abs(S(min(size(A,1),size(A,2)))) < 1.0e-16_8) then
         errorstop "In cond: Matrix A is singular"
      end if
      rslt = S(1)/S(min(size(A,1),size(A,2)))
   end function mycond
#endif

end module SHAKELie
