! 2023-06-13: starting implementing 

module half_explicit
   implicit none

   ! definition of integrator options
   type  :: half_explicit_options
      ! system is constrained?
      integer :: constrained = 0
      ! use stabilized index two formulation
      integer :: stab2 = 0
      ! mass matrix is constant?
      integer :: const_mass_matrix = 0
      ! mass matrix is a diagonal matrix? If == 1 then half_explicit_diag_M
      ! is used in order to calculate the mass matrix rather than
      ! half_explicit_M
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
      ! pertube initial values?
      integer :: pertube = 0 ! DEBUG
      ! use numerical approximation for Ct and Kt
      integer :: use_num_Ct = 1
      integer :: use_num_Kt = 1
      ! omit Ct and Kt resp. in the iteration matrix
      integer :: no_Ct = 0
      integer :: no_Kt = 0
      ! variables for error control of newton method (absolute and relativ tolerance)
      real(8)  :: atol = 1.0e-10_8
      real(8)  :: rtol = 1.0e-8_8
      integer  :: imax = 5
      ! integration span
      real(8) :: t0 = 0.0_8
      real(8) :: te = 1.0_8
      logical :: step_size_control = .false.
      logical :: local_error_control = .false.
      integer :: nsteps = 100
      real(8) :: facmax = 5.0_8
      real(8) :: facmin = 0.3_8
      real(8) :: fac    = 0.8_8
   end type half_explicit_options

   ! definition of integrator statistics
   type  :: half_explicit_statistics
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
   end type half_explicit_statistics

   ! definition of abstract problem type
   type, abstract :: half_explicit_problem
      ! variables that define the state of the integrator
      real(8)                             :: t = 0.0_8
      real(8)                             :: err = 0.0_8
      real(8), dimension(:), allocatable  :: local_est_err
      integer                             :: sizeq = 0
      real(8), dimension(:), allocatable  :: q
      integer                             :: sizev = 0
      real(8), dimension(:), allocatable  :: v
      real(8), dimension(:), allocatable  :: vd
      integer                             :: sizel = 0 ! number of constraints, irrelevant for this%opts%constrained == 0
      real(8), dimension(:), allocatable  :: l         ! Lagrange multipliers, only needed in the constrained case
      real(8), dimension(:), allocatable  :: eta       ! Auxiliar variables eta, only needed in the stabilized index-2 case
      ! integrator options
      type(half_explicit_options)      :: opts  ! TODO: rename to half_explicit_opts
      ! solver statistics
      type(half_explicit_statistics)   :: half_explicit_stats
      ! constant mass matrix (if opts%const_mass_matrix == 1)
      real(8), dimension(:,:), allocatable   :: half_explicit_const_M
      ! constant diagonal mass matrix (if opts%diag_mass_matrix == 1, also)
      real(8), dimension(:), allocatable     :: half_explicit_const_diag_M
      ! internal variables
      real(8), dimension(:,:), allocatable   :: half_explicit_A
      real(8), dimension(:), allocatable     :: half_explicit_c
      real(8), dimension(:), allocatable     :: half_explicit_d
      real(8), dimension(:), allocatable     :: half_explicit_variable_step_coeff
      integer                                :: half_explicit_s
      integer                                :: half_explicit_s_bar
      ! integer                                :: half_explicit_sizeA
      integer                                :: half_explicit_order
      integer                                :: half_explicit_order_variable_step
   contains
         ! function $M$ in $M(q) \dot v = -g(q,v,t)$
      procedure(half_explicit_M),              deferred :: half_explicit_M
         ! function $M$ in $M(q) \dot v = -g(q,v,t)$, diagonal case
      procedure(half_explicit_diag_M),         deferred :: half_explicit_diag_M
         ! function $g$ in $M(q) \dot v = -g(q,v,t)$
      procedure(half_explicit_g),              deferred :: half_explicit_g
         ! operation in the lie space $q_n * \exp(\widetilde{\Theta_n})$
      procedure(half_explicit_qlpexphDqtilde), deferred :: half_explicit_qlpexphDqtilde
         ! operation in the lie algebra $inversetilde([\tilde{v},\tilde{w}])$, may be a dummy, if system is unconstrained
      procedure(half_explicit_itlbtvtw),       deferred :: half_explicit_itlbtvtw
         ! tilde operator in the lie group (may be a dummy function, if half_explicit_num_Kt is not used)  ! TODO: This is not used anymore
      procedure(half_explicit_tilde),          deferred :: half_explicit_tilde ! TODO
         ! tangent damping matrix $C_t$
      procedure(half_explicit_Ct),             deferred :: half_explicit_Ct
         ! tangent stiffness matrix $K_t$
      procedure(half_explicit_Kt),             deferred :: half_explicit_Kt
         ! tangent stiffness matrix $K_t$ in the constrained case
         ! (may depend on the Lagrange multipliers)
      procedure(half_explicit_Kt_lambda),      deferred :: half_explicit_Kt_lambda
         ! tangent operator of the exponential map
      procedure(half_explicit_Tg),             deferred :: half_explicit_Tg
         ! transpose of the inverse of the tangent operator of the exponential map
      procedure(half_explicit_Tg_inv_T),       deferred :: half_explicit_Tg_inv_T
         ! norm    ! TODO: this is actually not used anymore
      procedure(half_explicit_norm),           deferred :: half_explicit_norm   ! TODO
         ! subroutine for output, is called after every integration step
      procedure(half_explicit_outputFunction), deferred :: half_explicit_outputFunction
         ! subroutine to initialize the problem
      procedure(half_explicit_init),           deferred :: half_explicit_init
         ! function $\Phi(q)$
      procedure(half_explicit_phi),            deferred :: half_explicit_phi
         ! function $B(q)$, where $B(q)v = d/dq \Phi(q) * (DL_q(e) v)$
      procedure(half_explicit_b),              deferred :: half_explicit_b
         ! function $F(q,v) = d/(dq) (B(q)v) * (DL_q(e) v)$
      procedure(half_explicit_Z),              deferred :: half_explicit_Z
         ! function Z, where $Zw = Z(q(Dq)) (v(Dq + BT(q)eta), T w)$
      procedure(half_explicit_matZ),           deferred :: half_explicit_matZ
         ! subroutine to calculate correct initial values for vd and a
      ! procedure                  :: half_explicit_calcInitial => half_explicit_calcInitial ! TODO: private
         ! subroutine to calculate correct initial values for vd and a, also apply pertubing, if needed
      procedure                  :: half_explicit_calcInitialConstrained => half_explicit_calcInitialConstrained
         ! subroutine to integrate one time step
      procedure                  :: half_explicit_solveTimeStep => half_explicit_solveTimeStep ! TODO: private
         ! subroutine to integrate one time step in a constrained system
      procedure                  :: half_explicit_solveConstrainedTimeStep => half_explicit_solveConstrainedTimeStep
         ! subroutine for integration
      procedure                  :: half_explicit_integrate => half_explicit_integrate
      procedure                  :: half_explicit_print_stats => half_explicit_print_stats
         ! Clean up the half_explicit_problem object (only internal variables
         ! are reset, half_explicit_alpha_m, half_explicit_alpha_f, half_explicit_beta and
         ! half_explicit_gamma are NOT reset)
      procedure                  :: half_explicit_cleanup => half_explicit_cleanup
   end type half_explicit_problem

   ! abstract interface of the problem
   abstract interface
      pure function half_explicit_M(this, q) result(rslt)
         import                                    :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function half_explicit_M

      pure function half_explicit_diag_M(this, q) result(rslt)
         import                            :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in) :: this
         real(8), dimension(:), intent(in) :: q
         ! result
         real(8), dimension(this%sizev)    :: rslt
      end function half_explicit_diag_M

      pure function half_explicit_g(this, q, v, t) result(rslt)
         import                              :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8), dimension(:), intent(in)   :: v
         real(8),               intent(in)   :: t
         ! result
         real(8), dimension(this%sizev)      :: rslt
      end function half_explicit_g

      pure function half_explicit_qlpexphDqtilde(this, q, h, Dq) result(rslt)
         import                             :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in)  :: this
         real(8), dimension(:), intent(in)  :: q
         real(8),               intent(in)  :: h
         real(8), dimension(:), intent(in)  :: Dq
         ! result
         real(8), dimension(this%sizeq)     :: rslt
      end function half_explicit_qlpexphDqtilde

      pure function half_explicit_itlbtvtw(this, v, w) result(rslt)
         import                              :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         real(8), dimension(:), intent(in)   :: w
         ! result
         real(8), dimension(this%sizev)      :: rslt
      end function half_explicit_itlbtvtw

      pure function half_explicit_tilde(this, v) result(rslt)
         import                              :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         ! result
         real(8), dimension(this%sizeq)      :: rslt
      end function half_explicit_tilde

      pure function half_explicit_Ct(this, q, v, t) result(rslt)
         import                                    :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function half_explicit_Ct

      pure function half_explicit_Kt(this, q, v, vd, t) result(rslt)
         import                                    :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8), dimension(:), intent(in)         :: vd
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function half_explicit_Kt

      pure function half_explicit_Kt_lambda(this, q, v, vd, l, t) result(rslt)
         import                                    :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8), dimension(:), intent(in)         :: vd
         real(8), dimension(:), intent(in)         :: l
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function half_explicit_Kt_lambda

      pure function half_explicit_Tg(this, h, dq) result(rslt)
         import                                    :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in)         :: this
         real(8),               intent(in)         :: h
         real(8), dimension(:), intent(in)         :: dq
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function half_explicit_Tg

      pure function half_explicit_Tg_inv_T(this, dq) result(rslt)
      import                                    :: half_explicit_problem
      ! input
      class(half_explicit_problem),     intent(in)         :: this
      real(8), dimension(:), intent(in)         :: dq
      ! result
      real(8), dimension(this%sizev,this%sizev) :: rslt
      end function half_explicit_Tg_inv_T
      
      pure function half_explicit_norm(this, v) result(rslt)
         import                  :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         ! result
         real(8)                             :: rslt
      end function half_explicit_norm

      subroutine half_explicit_outputFunction(this,info)
         import                           :: half_explicit_problem
         ! input
         class(half_explicit_problem),   intent(in)  :: this
         integer,             intent(in)  :: info
      end subroutine half_explicit_outputFunction

      subroutine half_explicit_init(this)
         import                            :: half_explicit_problem
         ! input/output
         class(half_explicit_problem), intent(inout)  :: this
      end subroutine half_explicit_init

      pure function half_explicit_phi(this,q) result(rslt)
         import                            :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in) :: this
         real(8), dimension(:), intent(in) :: q
         ! result
         real(8), dimension(this%sizel)    :: rslt
      end function half_explicit_phi

      pure function half_explicit_B(this,q) result(rslt)
         import                                     :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in)          :: this
         real(8), dimension(:), intent(in)          :: q
         ! result
         real(8), dimension(this%sizel,this%sizev)  :: rslt
      end function half_explicit_B

      pure function half_explicit_Z(this,q,v) result(rslt)
         import                            :: half_explicit_problem
         ! input
         class(half_explicit_problem),     intent(in) :: this
         real(8), dimension(:), intent(in) :: q
         real(8), dimension(:), intent(in) :: v
         ! result
         real(8), dimension(this%sizel)    :: rslt
      end function half_explicit_Z

      pure function half_explicit_matZ(this,q,v,T) result(rslt)
         import                              :: half_explicit_problem
         ! input
         class(half_explicit_problem),       intent(in) :: this
         real(8), dimension(:),   intent(in) :: q
         real(8), dimension(:),   intent(in) :: v
         real(8), dimension(:,:), intent(in) :: T
         ! result
         real(8), dimension(this%sizel, this%sizev) :: rslt
      end function half_explicit_matZ

   end interface

   contains

   subroutine half_explicit_calcInitialConstrained(this)
      implicit none
      class(half_explicit_problem),   intent(inout)  :: this  ! problem object
      ! internal variables
      integer                                                              :: i
      real(8), dimension(this%sizev+this%sizel)                            :: vdl
      integer, dimension(this%sizev+this%sizel)                            :: ipiv0   ! pivot vector for dgesv for LU factorization of MBB0, this is actually needed
      integer                                                              :: info    ! info flag for dgesv
      real(8), dimension(this%sizev+this%sizel, this%sizev+this%sizel)     :: MBB0    ! Matrix
      real(8), dimension(this%sizel, this%sizev)                           :: B0      ! $B(q_0)$
      !
      ! calulate $B(q_0)$
      B0 = this%half_explicit_b(this%q)
      ! count calls
      this%half_explicit_stats%nBcalls = this%half_explicit_stats%nBcalls + 1

      ! calculate $MBB0$
      if (this%opts%diag_mass_matrix == 1) then
         ! Set the mass matrix part to zero beforehand
         MBB0(1:this%sizev, 1:this%sizev) = 0.0_8
         if (this%opts%const_mass_matrix == 1) then
            MBB0(1:this%sizev, 1) = this%half_explicit_const_diag_M
         else
            MBB0(1:this%sizev, 1) = this%half_explicit_diag_M(this%q)
         end if
         do concurrent (i=2:this%sizev)
            MBB0(i,i) = MBB0(i,1)
            MBB0(i,1) = 0.0_8
         end do
      else
         if (this%opts%const_mass_matrix == 1) then
            MBB0(1:this%sizev, 1:this%sizev) = this%half_explicit_const_M
         else
            MBB0(1:this%sizev, 1:this%sizev) = this%half_explicit_M(this%q)
         end if
      end if
      MBB0(this%sizev+1:this%sizev+this%sizel, 1:this%sizev)                       =           B0
      MBB0(1:this%sizev,                       this%sizev+1:this%sizev+this%sizel) = transpose(B0)
      MBB0(this%sizev+1:this%sizev+this%sizel, this%sizev+1:this%sizev+this%sizel) = 0.0_8

      ! we need to solve a linear equation, first calulate the rhs
      vdl(1:this%sizev)                       = -this%half_explicit_g(this%q, this%v, this%t)
      vdl(this%sizev+1:this%sizev+this%sizel) = -this%half_explicit_Z(this%q, this%v)
      ! count calls
      this%half_explicit_stats%ngcalls = this%half_explicit_stats%ngcalls + 1
      ! then solve the system
      call dgesv(          &! solve the System A*X=B and save the result in B
                  this%sizev+this%sizel,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                  1,       &! number of right hand sides (=size(B,2))
                  MBB0,    &! matrix A
                  this%sizev+this%sizel,   &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                  ipiv0,   &! integer pivot vector; it is not needed
                  vdl,     &! matrix B
                  this%sizev+this%sizel,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                  info)     ! integer information flag
      ! Now vdl actually contains $\dot v(t_0)$ and $\lambda(t_0)
      if (info .ne. 0)  print*, "calcInitialConstrained:  dgesv sagt info=", info ! TODO

      ! apply the calculated values
      this%vd = vdl(   1:   this%sizev)
      this%l  = vdl(this%sizev+1:this%sizev+this%sizel)
   end subroutine half_explicit_calcInitialConstrained

   ! subroutine for integrating one time step
   subroutine half_explicit_solveConstrainedTimeStep(this, t1)
      implicit none
      class(half_explicit_problem), intent(inout)  :: this ! problem object
      real(8)          , intent(in   )  :: t1   ! $t_{n+1}$

      ! internal integer variables
      integer                           :: i    ! for iteration
      integer                           :: j    ! for iteration
      integer, dimension(this%sizev)    :: ipiv ! pivot vector for dgesv
      integer                           :: info ! info flag for dgesv

      ! intenal real variables
      real(8)                                     :: h    ! step size
      real(8), dimension(this%sizev+this%sizel)   :: dVl
      ! real(8), dimension(this%sizev+this%sizel)   :: dVl_variable_step
      real(8), dimension(this%sizev)              :: Vcrr
      real(8), dimension(this%half_explicit_s_bar+1,   this%sizeq) :: Qn
      real(8), dimension(                            this%sizeq) :: Qn_local_error
      real(8), dimension(this%half_explicit_s_bar+1,   this%sizev) :: Vn
      real(8), dimension(                            this%sizev) :: Vn_local_error
      real(8), dimension(this%half_explicit_s_bar+1,   this%sizev) :: Thetan
      real(8), dimension(this%half_explicit_s_bar, this%sizev) :: dThetan
      real(8), dimension(                            this%sizev) :: Thetan_local_error
      real(8), dimension(this%half_explicit_s_bar, this%sizev) :: dVn
      ! real(8), dimension(                            this%sizev) :: dVn_variable_step
      real(8), dimension(this%half_explicit_s_bar, this%sizel) :: Lambdan
      ! real(8), dimension(                            this%sizel) :: Lambdan_variable_step
      real(8), dimension(this%sizev,  this%sizev) :: M
      real(8), dimension(this%sizel,  this%sizev) :: B0
      real(8), dimension(this%sizel,  this%sizev) :: B1
      ! real(8), dimension(this%sizel,  this%sizev) :: B1_variable_step
      real(8), dimension(this%sizev+this%sizel, this%sizev+this%sizel) :: MBB0
      ! real(8), dimension(this%sizev+this%sizel, this%sizev+this%sizel) :: MBB0_variable_step

      ! calculation of step size $h$
      h = t1 - this%t

      ! first stage
      ! print *, "fist stage"
      Thetan(1,1:this%sizev) = 0.0_8
      Qn(1,1:this%sizeq) = this%q
      ! Qn(1,1:this%sizeq) = this%half_explicit_qlpexphDqtilde(this%q, 1.0_8, Thetan(1,:))
      Vn(1,1:this%sizev) = this%v
      dThetan(1,:) = matmul(transpose(this%half_explicit_Tg_inv_T(Thetan(1,:))), Vn(1,:))
      Lambdan(1,:) = this%l
      ! calculate $\dotV_{m1}$ --> dVn(1,:)
      ! 1. MASS MATRIX
      if (this%opts%diag_mass_matrix == 1) then
         ! Set the mass matrix part to zero beforehand
         M = 0.0_8
         if (this%opts%const_mass_matrix == 1) then
            M(1:this%sizev, 1) = this%half_explicit_const_diag_M
         else
            M(1:this%sizev, 1) = this%half_explicit_diag_M(Qn(1,:))
         end if
         do concurrent (i=2:this%sizev)
            M(i,i) = M(i,1)
            M(i,1) = 0.0_8
         end do
      else
         if (this%opts%const_mass_matrix == 1) then
            M(:,:) = this%half_explicit_const_M
         else
            M(:,:) = this%half_explicit_M(Qn(1,:))
         end if
      end if
      ! we need to solve a linear equation, first calulate the rhs
      ! 2. RIGHT HAND SIDE
      dVn(1,:) = -(this%half_explicit_g(Qn(1,:), this%v, this%t + this%half_explicit_c(1)*h)+matmul(transpose(this%half_explicit_B(Qn(1,:))),this%l))
      ! count calls
      this%half_explicit_stats%ngcalls = this%half_explicit_stats%ngcalls + 1
      this%half_explicit_stats%nBcalls = this%half_explicit_stats%nBcalls + 1
      ! 3. SOLVE THE SYSTEM
      ! print *, "fist stage - linear system"
      call dgesv(         &! solve the System A*X=B and save the result in B
                  this%sizev,      &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                  1,       &! number of right hand sides (=size(B,2))
                  M(:,:),  &! matrix A
                  this%sizev,      &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                  ipiv,    &! integer pivot vector; it is not needed
                  dVn(1,:),&! matrix B
                  this%sizev,      &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                  info)     ! integer information flag
      ! Now dVn(1,:) actually contains $\dot v(t_0)$
      if (info .ne. 0)  print*, "TimeStep--fist step: dgesv sagt info=", info ! TODO

      ! second stage (only Theta_2 and Q_2)
      ! print *, "second stage - Theta and Q"
      Thetan(2,:) = h * this%half_explicit_A(2,1) * dThetan(1,:)
      Qn(2,:) = this%half_explicit_qlpexphDqtilde(this%q, 1.0_8, Thetan(2,:))
      ! Vn(2,:) = this%v + h * this%half_explicit_A(2,1) * dVn(1,:)

      ! following stages
      do i=2,this%half_explicit_s_bar + 1
         ! print *, "stage ", i
         ! evaluate V_i
         Vn(i,:) = this%v
         do j = 1,i-1
            Vn(i,:) = Vn(i,:) + h * this%half_explicit_A(i,j) * dVn(j,:)
         end do

         if ( i < this%half_explicit_s_bar + 1 ) then
            ! evaluate dTheta_i
            dThetan(i,:) = matmul(transpose(this%half_explicit_Tg_inv_T(Thetan(i,:))), Vn(i,:))

            ! evaluate Theta_{i+1} and Q_{i+1}
            Thetan(i+1,:) = 0.0_8
            do j=1,i
               Thetan(i+1,:)=Thetan(i+1,:) + h * this%half_explicit_A(i+1,j) * dThetan(j,:)
            end do
            Qn(i+1,:) = this%half_explicit_qlpexphDqtilde(this%q, 1.0_8, Thetan(i+1,:))

            ! calculate $\dotV_{m,i}$ and $\Lambda_{m,i}$ --> dVn(i,:), Lambdan(i,:)
            ! 1. MATRIX OF COEFFICIENTS
            MBB0 = 0.0_8
            if (this%opts%diag_mass_matrix == 1) then
               if (this%opts%const_mass_matrix == 1) then
                  MBB0(1:this%sizev,1) = this%half_explicit_const_diag_M
               else
                  MBB0(1:this%sizev,1) = this%half_explicit_diag_M(Qn(i,:))
               end if
               do concurrent (i=2:this%sizev)
                  MBB0(i,i) = MBB0(i,1)
                  MBB0(i,1) = 0.0_8
               end do
            else
               if (this%opts%const_mass_matrix == 1) then
                  MBB0(1:this%sizev,1:this%sizev) = this%half_explicit_const_M
               else
                  MBB0(1:this%sizev,1:this%sizev) = this%half_explicit_M(Qn(i,:))
               end if
            end if
            ! calulate $B(Qn(i))$ and $B(Qn(i+1))$
            B0 = this%half_explicit_B(Qn(i,:))
            B1 = this%half_explicit_B(Qn(i+1,:))
            MBB0(1:this%sizev, this%sizev+1:this%sizev+this%sizel) = transpose(B0)
            MBB0(this%sizev+1:this%sizev+this%sizel, 1:this%sizev) = h * this%half_explicit_A(i+1,i) * B1
            ! count calls
            this%half_explicit_stats%nBcalls = this%half_explicit_stats%nBcalls + 2
            ! we need to solve a linear equation, first calulate the rhs
            ! 2. RIGHT HAND SIDE
            ! V_{i+1} without the last value of dV_i, which we are going to evaluate NOW
            Vcrr = this%v
            do j = 1,i-1
               Vcrr = Vcrr + h * this%half_explicit_A(i+1,j) * dVn(j,:)
            end do
            dVl(1:this%sizev) = -this%half_explicit_g(Qn(i,:), Vn(i,:), this%t + h * this%half_explicit_c(i))
            dVl(this%sizev+1:this%sizev+this%sizel) = -matmul(B1, Vcrr)
            ! count calls
            this%half_explicit_stats%ngcalls = this%half_explicit_stats%ngcalls + 1
            ! then solve the system
            ! print *, "stage ", i, " - linear system"
            call dgesv(                       &! solve the System A*X=B and save the result in B
                       this%sizev+this%sizel, &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                       1,                     &! number of right hand sides (=size(B,2))
                       MBB0,                  &! matrix A
                       this%sizev+this%sizel, &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                       ipiv,                  &! integer pivot vector; it is not needed
                       dVl,                   &! matrix B
                       this%sizev+this%sizel, &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                       info)                   ! integer information flag
                     ! if (info .ne. 0)  print*, "TimeStep--following steps: dgesv sagt info=", info ! TODO
            dVn(i,:)     = dVl(1:this%sizev)
            Lambdan(i,:) = dVl(this%sizev+1:this%sizev+this%sizel)
         end if

         if ( (this%opts%local_error_control) .and. (i == this%half_explicit_s_bar + 1) ) then
            ! evaluate Theta_error_control and Q_error_control
            Thetan_local_error = 0.0_8
            do j=1,i-1
               Thetan_local_error = Thetan_local_error + h * this%half_explicit_variable_step_coeff(j) * dThetan(j,:)
            end do
            Qn_local_error = this%half_explicit_qlpexphDqtilde(this%q, 1.0_8, Thetan_local_error)
            Vn_local_error = this%v
            do j = 1,i-1
               Vn_local_error = Vn_local_error + h * this%half_explicit_variable_step_coeff(j) * dVn(j,:)
            end do

            ! if (allocated(this%local_est_err)) deallocate(this%local_est_err)
            ! allocate(this%local_est_err(this%sizeq + this%sizev))
            this%local_est_err(1:this%sizeq) = abs(Qn(this%half_explicit_s+1,:) - Qn_local_error)
            this%local_est_err(this%sizeq+1:this%sizeq+this%sizev) = abs(Vn(this%half_explicit_s+1,:) - Vn_local_error)

            this%err = 0.0_8
            ! evaluate local error
            do j = 1,this%sizeq+this%sizev
               if (j < this%sizeq + 1) then
                  this%err = this%err + ((Qn(this%half_explicit_s+1,j)-Qn_local_error(j))/(this%opts%atol + max(abs(Qn(this%half_explicit_s+1,j)),abs(Qn(1,j))) * this%opts%rtol))**2
               else
                  this%err = this%err + ((Vn(this%half_explicit_s+1,j-this%sizeq)-Vn_local_error(j-this%sizeq))/(this%opts%atol + max(abs(Vn(this%half_explicit_s+1,j-this%sizeq)),abs(Vn(1,j-this%sizeq))) * this%opts%rtol))**2
               end if
            end do
            this%err = sqrt(this%err / (this%sizeq + this%sizev))
         end if
      end do
      ! print *, "saving solutions"
      this%q = Qn(this%half_explicit_s+1,:)
      this%v = Vn(this%half_explicit_s+1,:)
      this%l = 0.0_8
      do i = 1, this%half_explicit_s_bar
         this%l = this%l + this%half_explicit_d(i) * Lambdan(i,:)
      end do
      this%t = t1
      ! print *, "exit"
      print *, "err=", this%err
   end subroutine half_explicit_solveConstrainedTimeStep

   ! subroutine for integrating one time step
   subroutine half_explicit_solveTimeStep(this, t1)
      implicit none
      class(half_explicit_problem), intent(inout)  :: this  ! problem object
      real(8)          , intent(in   )  :: t1    ! $t_{n+1}$

      ! internal integer variables
      integer                           :: i    ! for iteration
      integer                           :: j    ! for iteration
      integer, dimension(this%sizev)    :: ipiv ! pivot vector for dgesv
      integer                           :: info ! info flag for dgesv

      ! intenal real variables
      real(8)                           :: h       ! step size
      real(8), dimension(this%half_explicit_s_bar+1, this%sizeq) :: Qn
      real(8), dimension(this%half_explicit_s_bar+1, this%sizev) :: Vn
      real(8), dimension(                          this%sizeq) :: Qn_local_error
      real(8), dimension(                          this%sizev) :: Vn_local_error
      real(8), dimension(this%half_explicit_s_bar+1, this%sizev) :: Thetan
      real(8), dimension(this%half_explicit_s_bar+1, this%sizev) :: dThetan
      real(8), dimension(this%half_explicit_s_bar+1, this%sizev) :: dVn
      real(8), dimension(                          this%sizev) :: Thetan_local_error
      real(8), dimension(this%sizev,  this%sizev) :: M

      ! internal logical variables

      h = t1 - this%t

      ! first stage
      Thetan(1,1:this%sizev) = 0.0_8
      Qn(1,1:this%sizeq) = this%half_explicit_qlpexphDqtilde(this%q, 1.0_8, Thetan(1,:))
      Vn(1,1:this%sizev) = this%v
      dThetan(1,:) = matmul(transpose(this%half_explicit_Tg_inv_T(Thetan(1,:))), Vn(1,:))
      ! calculate $\dotV_{m1}$ --> dVn(1,:)
      if (this%opts%diag_mass_matrix == 1) then
         if (this%opts%const_mass_matrix == 1) then
            dVn(1,:) = -this%half_explicit_g(Qn(1,:), Vn(1,:), this%t + this%half_explicit_c(1)*h)/this%half_explicit_const_diag_M
         else
            dVn(1,:) = -this%half_explicit_g(Qn(1,:), Vn(1,:), this%t + this%half_explicit_c(1)*h)/this%half_explicit_diag_M(Qn(1,:))
         end if
         ! count calls
         this%half_explicit_stats%ngcalls = this%half_explicit_stats%ngcalls + 1
      else
         if (this%opts%const_mass_matrix == 1) then
            M(:,:) = this%half_explicit_const_M
         else
            M(:,:) = this%half_explicit_M(Qn(1,:))
         end if
         ! we need to solve a linear equation, first calulate the rhs
         dVn(1,:) = -this%half_explicit_g(Qn(1,:), this%v, this%t + this%half_explicit_c(1)*h)
         ! count calls
         this%half_explicit_stats%ngcalls = this%half_explicit_stats%ngcalls + 1
         ! then solve the system
         call dgesv(         &! solve the System A*X=B and save the result in B
                    this%sizev,      &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                    1,       &! number of right hand sides (=size(B,2))
                    M(:,:),  &! matrix A
                    this%sizev,      &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                    ipiv,    &! integer pivot vector; it is not needed
                    dVn(1,:),&! matrix B
                    this%sizev,      &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                    info)     ! integer information flag
         ! Now dVn(1,:) actually contains $\dot v(t_0)$
         if (info .ne. 0)  print*, "TimeStep--fist step: dgesv sagt info=", info ! TODO
      end if

      ! following stages
      do i=2,this%half_explicit_s_bar+1
         Thetan(i,:) = 0.0_8
         Vn(i,:) = this%v
         if ((this%opts%step_size_control) .and. (i == this%half_explicit_s_bar+1)) then
            Thetan_local_error = 0.0_8
            Vn_local_error = this%v
         end if
         do j=1,i-1
            Thetan(i,:)=Thetan(i,:) + h * this%half_explicit_A(i,j) * dThetan(j,:)
            Vn(i,:) = Vn(i,:) + h * this%half_explicit_A(i,j) * dVn(j,:)
            if ((this%opts%step_size_control) .and. (i == this%half_explicit_s_bar+1)) then
               Thetan_local_error=Thetan_local_error + h * this%half_explicit_variable_step_coeff(j) * dThetan(j,:)
               Vn_local_error = Vn_local_error + h * this%half_explicit_variable_step_coeff(j) * dVn(j,:)
            end if
         end do
         Qn(i,:) = this%half_explicit_qlpexphDqtilde(this%q, 1.0_8, Thetan(i,:))
         if ((this%opts%step_size_control) .and. (i == this%half_explicit_s_bar+1)) then
            Qn_local_error = this%half_explicit_qlpexphDqtilde(this%q, 1.0_8, Thetan_local_error)

            this%err = 0.0_8
            ! evaluate local error
            do j = 1,this%sizeq+this%sizev
               if (j < this%sizeq + 1) then
                  this%err = this%err + ((Qn(i,j)-Qn_local_error(j))/(this%opts%atol + max(abs(Qn(i,j)),abs(Qn(1,j))) * this%opts%rtol))**2
               else
                  this%err = this%err + ((Vn(i,j-this%sizeq)-Vn_local_error(j-this%sizeq))/(this%opts%atol + max(abs(Vn(i,j-this%sizeq)),abs(Vn(1,j-this%sizeq))) * this%opts%rtol))**2
               end if
            end do
            this%err = sqrt(this%err / (this%sizeq + this%sizev))
            ! this%err = sqrt(this%err / (this%sizev))
         end if
         dThetan(i,:) = matmul(transpose(this%half_explicit_Tg_inv_T(Thetan(i,:))), Vn(i,:))
         ! calculate $\dotV_{m,i}$ --> dVn(i,:)
         if (this%opts%diag_mass_matrix == 1) then
            if (this%opts%const_mass_matrix == 1) then
               dVn(i,:) = -this%half_explicit_g(Qn(i,:), Vn(i,:), this%t + this%half_explicit_c(i)*h)/this%half_explicit_const_diag_M
            else
               dVn(i,:) = -this%half_explicit_g(Qn(i,:), Vn(i,:), this%t + this%half_explicit_c(i)*h)/this%half_explicit_diag_M(Qn(i,:))
            end if
            ! count calls
            this%half_explicit_stats%ngcalls = this%half_explicit_stats%ngcalls + 1
         else
            if (this%opts%const_mass_matrix == 1) then
               M(:,:) = this%half_explicit_const_M
            else
               M(:,:) = this%half_explicit_M(Qn(i,:))
            end if
            ! we need to solve a linear equation, first calulate the rhs
            dVn(i,:) = -this%half_explicit_g(Qn(i,:), Vn(i,:), this%t + h * this%half_explicit_c(i))
            ! count calls
            this%half_explicit_stats%ngcalls = this%half_explicit_stats%ngcalls + 1
            ! then solve the system
            call dgesv(         &! solve the System A*X=B and save the result in B
                       this%sizev,      &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                       1,       &! number of right hand sides (=size(B,2))
                       M(:,:),  &! matrix A
                       this%sizev,      &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                       ipiv,    &! integer pivot vector; it is not needed
                       dVn(i,:),&! matrix B
                       this%sizev,      &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                       info)     ! integer information flag
            ! Now dVn(1,:) actually contains $\dot v(t_0)$
            if (info .ne. 0)  print*, "TimeStep--fist step: dgesv sagt info=", info ! TODO
         end if
         
      end do

      this%q = Qn(this%half_explicit_s+1,:)
      this%v = Vn(this%half_explicit_s+1,:)
      this%t = t1

   end subroutine half_explicit_solveTimeStep

   ! subroutine to integrate
   subroutine half_explicit_integrate(this)
      implicit none
      class(half_explicit_problem), intent(inout) :: this
      integer                          :: n = 1 ! needed for iteration
      real(8)                          :: h     ! step size $h$
      real(8)                          :: mulfac
      real(8)                          :: h_new ! step size $h$ of step size control
      real(8)                          :: h_old
      real(8), dimension(this%sizeq)   :: q_old
      real(8), dimension(this%sizev)   :: v_old
      real(8), dimension(this%sizel)   :: l_old
      real(8)                          :: t1 = 0.0_8 ! next time $t_{n+1}$
      real(8)                          :: t_old
      real(4)                          :: times(2), time ! for dtime

      integer :: i = 0 ! counter

      ! problem initialization
      call this%half_explicit_init()

      ! initialize output function
      call this%half_explicit_outputFunction(0)

      ! Calculate step size $h$
      h = (this%opts%te - this%opts%t0)/this%opts%nsteps

      ! Set stats of solver to zero
      this%half_explicit_stats%newt_steps_curr = 0
      this%half_explicit_stats%newt_steps_sum  = 0
      this%half_explicit_stats%newt_steps_max  = 0
      this%half_explicit_stats%newt_steps_avg  = 0
      this%half_explicit_stats%ngcalls = 0
      this%half_explicit_stats%nBcalls = 0

      ! if mass matrix is constant, calculate it
      if (this%opts%const_mass_matrix == 1) then
         if (this%opts%diag_mass_matrix == 1) then
            if (.not. allocated(this%half_explicit_const_diag_M)) then
               allocate(this%half_explicit_const_diag_M(this%sizev))
            end if
            this%half_explicit_const_diag_M = this%half_explicit_diag_M(this%q)
         else
            if (.not. allocated(this%half_explicit_const_M)) then
               allocate(this%half_explicit_const_M(this%sizev,this%sizev))
            end if
            this%half_explicit_const_M = this%half_explicit_M(this%q)
         end if
      end if

      ! start stopwatch
      time = dtime(times)

      if (this%opts%constrained == 0) then

         print *, 'This is an explicit Runge-Kutta method of order ', this%half_explicit_order
         
         ! no need for vd, but needs initialization
         this%vd = 0.0_8
         ! output for the first time
         call this%half_explicit_outputFunction(1)
         
         ! integration loop
         do while (t1 <= this%opts%te .or. n <= this%opts%nsteps)

            ! Calculate the next time $t_{n+1}$
            if (this%opts%step_size_control) then
               if (n == 1) then
                  t1 = this%opts%t0 + h
                  h_new = h
                  n = this%opts%nsteps + 1
                  ! save in auxiliary variable current solution
                  q_old = this%q
                  v_old = this%v
                  l_old = this%l
                  t_old = this%t
                  h_old = h_new
               else
                  ! checking local error
                  ! always evaluate new step size h
                  if (this%err > 1) then
                     ! if err > 1 --> not accepted, this%t = t_old
                     this%q = q_old
                     this%v = v_old
                     this%l = l_old
                     this%t = t_old
                     mulfac = min(this%opts%facmax, max(this%opts%facmin, (this%opts%fac * (1/this%err)**(1/(min(this%half_explicit_order,this%half_explicit_order_variable_step)+1.0_8)))))
                     h_new = h_old * mulfac
                     t1 = this%t + h_new
                     h_old = h_new
                  else
                     ! if err <= 1 --> accepted, t1 = t1
                     mulfac = min(this%opts%facmax, max(this%opts%facmin, (this%opts%fac * (1/this%err)**(1/(min(this%half_explicit_order,this%half_explicit_order_variable_step)+1.0_8)))))
                     h_new = h_old * mulfac
                     t1 = this%t + h_new

                     ! save in auxiliary variable current solution
                     q_old = this%q
                     v_old = this%v
                     l_old = this%l
                     t_old = this%t
                     h_old = h_new

                     ! output normally
                     call this%half_explicit_outputFunction(1)
                  end if
               end if
            else
               t1 = this%opts%t0 + n*h
               n = n+1
               if ( n > 2 ) then
                  ! output normally
                  call this%half_explicit_outputFunction(1)
               endif
            end if

            ! ! Check the error!!
            ! t1 = this%opts%t0 + n*h
            ! n = n+1
            ! if ( n > 2 ) then
            !    ! output normally
            !    call this%half_explicit_outputFunction(1)
            ! endif
            ! print *, 'err=', this%err

            ! ! DEBUG
            ! print *, 't1=', t1
            ! print *, 'hn=', h_new
            ! if (n>1) then
            !    print *, 'er=', this%err
            ! endif
            ! solve time step
            call this%half_explicit_solveTimeStep(t1)
         end do
      else
         ! print *, "initial integration"
         ! calculate correct initial values
         call this%half_explicit_calcInitialConstrained()
         ! output for the first time
         call this%half_explicit_outputFunction(1)

         ! integration loop
         i = 0
         do while (t1 <= this%opts%te .or. n <= this%opts%nsteps)

            ! Calculate the next time $t_{n+1}$
            if (this%opts%step_size_control) then
               if (n == 1) then
                  t1 = this%opts%t0 + h
                  h_new = h
                  n = this%opts%nsteps + 1
                  ! save in auxiliary variable current solution
                  q_old = this%q
                  v_old = this%v
                  l_old = this%l
                  t_old = this%t
                  h_old = h_new
               else
                  ! checking local error
                  ! always evaluate new step size h
                  if (this%err > 1) then
                     ! if err > 1 --> not accepted, this%t = t_old
                     this%q = q_old
                     this%v = v_old
                     this%l = l_old
                     this%t = t_old
                     mulfac = min(this%opts%facmax, max(this%opts%facmin, (this%opts%fac * (1/this%err)**(1/(min(this%half_explicit_order,this%half_explicit_order_variable_step)+1.0_8)))))
                     h_new = h_old * mulfac
                     t1 = this%t + h_new
                     h_old = h_new
                  else
                     ! if err <= 1 --> accepted, t1 = t1
                     mulfac = min(this%opts%facmax, max(this%opts%facmin, (this%opts%fac * (1/this%err)**(1/(min(this%half_explicit_order,this%half_explicit_order_variable_step)+1.0_8)))))
                     h_new = h_old * mulfac
                     t1 = this%t + h_new

                     h_old = h_new
                     q_old = this%q
                     v_old = this%v
                     l_old = this%l
                     t_old = this%t

                     ! output normally
                     call this%half_explicit_outputFunction(1)
                  end if
               end if
            else
               t1 = this%opts%t0 + n*h
               n = n+1
               if ( n > 2 ) then
                  ! output normally
                  call this%half_explicit_outputFunction(1)
               endif
            end if

            ! ! Check the error!!
            ! t1 = this%opts%t0 + n*h
            ! n = n+1
            ! if ( n > 2 ) then
            !    ! output normally
            !    call this%half_explicit_outputFunction(1)
            ! endif
            ! print *, 'err=', this%err

            ! DEBUG
            ! print *, 't1=', t1
            ! print *, 'hn=', h_new
            ! if (n>1) then
            !    print *, 'er=', this%err
            ! endif
            ! print *, "starting time step ", i+1
            call this%half_explicit_solveConstrainedTimeStep(t1)
            i = i + 1
            ! print *, "time step ", i, " complete!"
         end do
      end if
      ! stop stopwatch
      this%half_explicit_stats%time = dtime(times) - time
      ! output to terminate
      call this%half_explicit_outputFunction(99)

   end subroutine half_explicit_integrate

   subroutine half_explicit_print_stats(this)
      ! input
      class(half_explicit_problem), intent(in)  :: this
      !
      print *, 'time:          ', this%half_explicit_stats%time
      print *, '#calls of g:   ', this%half_explicit_stats%ngcalls
      print *, '#calls of B:   ', this%half_explicit_stats%nBcalls
      print *, 'newt_steps_max:', this%half_explicit_stats%newt_steps_max
      print *, 'newt_steps_avg:', this%half_explicit_stats%newt_steps_avg
   end subroutine half_explicit_print_stats

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

   subroutine half_explicit_cleanup(this)
      implicit none
      ! input/output
      class(half_explicit_problem), intent(inout)  :: this
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
      this%opts%step_size_control = .false.
      this%opts%nsteps = 100
      !
      this%half_explicit_stats%newt_steps_curr = 0
      this%half_explicit_stats%newt_steps_max = 0
      this%half_explicit_stats%newt_steps_avg = 0.0_8
      this%half_explicit_stats%ngcalls = 0
      this%half_explicit_stats%nBcalls = 0
      this%half_explicit_stats%time = 0.0_8
      !
      this%t = 0.0_8
      this%sizeq = 0
      this%sizev = 0
      this%sizel = 0
      !
      if (allocated(this%q))  deallocate(this%q)
      if (allocated(this%v))  deallocate(this%v)
      if (allocated(this%vd)) deallocate(this%vd)
      if (allocated(this%l))  deallocate(this%l)
      if (allocated(this%eta)) deallocate(this%eta)
      if (allocated(this%half_explicit_const_M)) deallocate(this%half_explicit_const_M)
      if (allocated(this%half_explicit_const_diag_M)) deallocate(this%half_explicit_const_diag_M)
      if (allocated(this%opts%jour)) deallocate(this%opts%jour)
   end subroutine half_explicit_cleanup

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

end module half_explicit
