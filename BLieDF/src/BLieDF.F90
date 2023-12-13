
module BLieDF
   implicit none

   ! coefficients for allowed k_bdf
   real(8), dimension(3,2), parameter :: alphas = reshape(&
      [1.0_8, -1.0_8, 0.0_8,& ! k=1 ! DEBUG
       1.5_8, -2.0_8, 0.5_8]& ! k=2
      , [3,2])


   ! definition of integrator options
   type  :: BLieDF_options
      ! system is constrained?
      integer :: constrained = 0
      ! use stabilized index two formulation
      integer :: stab2 = 0
      ! mass matrix is constant?
      integer :: const_mass_matrix = 0
      ! mass matrix is a diagonal matrix? If == 1 then BLieDF_diag_M
      ! is used in order to calculate the mass matrix rather than
      ! BLieDF_M
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
   end type BLieDF_options

   ! definition of integrator statistics
   type  :: BLieDF_statistics
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
   end type BLieDF_statistics

   ! definition of abstract problem type
   type, abstract :: BLieDF_problem
      ! variables that define the state of the integrator
      real(8)                             :: t = 0.0_8
      integer                             :: sizeq = 0
      real(8), dimension(:), allocatable  :: q
      integer                             :: sizev = 0
      real(8), dimension(:), allocatable  :: v
      real(8), dimension(:), allocatable  :: vd
      real(8), dimension(:), allocatable  :: Dq ! this is just for output
      real(8), dimension(:,:), allocatable  :: v_old ! these are just for input
      real(8), dimension(:,:), allocatable  :: Dq_old !  -- "" --
      integer                             :: sizel = 0 ! number of constraints, irrelevant for this%opts%constrained == 0
      real(8), dimension(:), allocatable  :: l ! Lagrange multipliers, only needed in the constrained case
      real(8), dimension(:), allocatable  :: eta ! Auxiliar variables eta, only needed in the stabilized index-2 case
      ! integrator options
      type(BLieDF_options)                :: opts
      ! solver statistics
      type(BLieDF_statistics)             :: BLieDF_stats
      ! constant mass matrix (if opts%const_mass_matrix == 1)
      real(8), dimension(:,:), allocatable   :: BLieDF_const_M
      ! constant diagonal mass matrix (if opts%diag_mass_matrix == 1, also)
      real(8), dimension(:), allocatable     :: BLieDF_const_diag_M
      ! internal variables:
      ! Number of steps, coefficients
      integer                             :: k_bdf = -1
      real(8), dimension(:), allocatable  :: alpha
      real(8), dimension(:), allocatable  :: gamma
   ! definition of deferred procedures
   contains
         ! function $M$ in $M(q) \dot v = -g(q,v,t)$
      procedure(BLieDF_M),              deferred :: BLieDF_M
         ! function $M$ in $M(q) \dot v = -g(q,v,t)$, diagonal case
      procedure(BLieDF_diag_M),         deferred :: BLieDF_diag_M
         ! function $g$ in $M(q) \dot v = -g(q,v,t)$
      procedure(BLieDF_g),              deferred :: BLieDF_g
         ! operation in the lie space $q_n * \exp(h\cdot\widetilde{\Delta q_n})$
      procedure(BLieDF_qlpexphDqtilde), deferred :: BLieDF_qlpexphDqtilde
         ! operation in the lie algebra $inversetilde([\tilde{v},\tilde{w}])$, may be a dummy, if system is unconstrained
      procedure(BLieDF_itlbtvtw),         deferred :: BLieDF_itlbtvtw
         ! tangent damping matrix $C_t$
      procedure(BLieDF_Ct),             deferred :: BLieDF_Ct
         ! tangent stiffness matrix $K_t$
      procedure(BLieDF_Kt),             deferred :: BLieDF_Kt
         ! tangent stiffness matrix $K_t$ in the constrained case
         ! (may depend on the Lagrange multipliers)
      procedure(BLieDF_Kt_lambda),      deferred :: BLieDF_Kt_lambda
         ! tangent operator of the exponential map
      procedure(BLieDF_Tg),             deferred :: BLieDF_Tg
         ! subroutine for output, is called after every integration step
      procedure(BLieDF_outputFunction), deferred :: BLieDF_outputFunction
         ! subroutine to initialize the problem
      procedure(BLieDF_init),           deferred :: BLieDF_init
         ! function $\Phi(q)$
      procedure(BLieDF_phi),            deferred :: BLieDF_phi
         ! function $B(q)$, where $B(q)v = d/dq \Phi(q) * (DL_q(e) v)$
      procedure(BLieDF_b),              deferred :: BLieDF_b
         ! function $F(q,v) = d/(dq) (B(q)v) * (DL_q(e) v)$
      procedure(BLieDF_Z),              deferred :: BLieDF_Z
         ! function Z, where $Zw = Z(q(Dq)) (v(Dq + BT(q)eta), T w)$
      procedure(BLieDF_matZ),           deferred :: BLieDF_matZ
         ! subroutine to calculate correct initial values for vd and a
      procedure                  :: BLieDF_calcInitial => BLieDF_calcInitial ! TODO: private
         ! subroutine to calculate correct initial values for vd and a
      procedure                  :: BLieDF_calcInitialConstrained => BLieDF_calcInitialConstrained
         ! subroutine to integrate one time step
      procedure                  :: BLieDF_solveTimeStep => BLieDF_solveTimeStep ! TODO: private
         ! subroutine to integrate one time step in a constrained system
      procedure                  :: BLieDF_solveConstrainedTimeStep => BLieDF_solveConstrainedTimeStep
         ! subroutine to integrate one time step in a constrained system
      procedure                  :: BLieDF_solveConstrainedTimeStep_stab2 => BLieDF_solveConstrainedTimeStep_stab2
         ! subroutine for integration
      procedure                  :: BLieDF_integrate => BLieDF_integrate
         ! subroutine for numerical approximation of Ct
      procedure                  :: BLieDF_num_Ct => BLieDF_num_Ct
         ! subroutine for numerical approximation of Kt
      procedure                  :: BLieDF_num_Kt => BLieDF_num_Kt
         ! subroutine for numerical approximation of Kt in the constrained case
      procedure                  :: BLieDF_num_Kt_lambda => BLieDF_num_Kt_lambda
         ! In order to use the numerical approximations of Ct and Kt
         ! respectively, the procedures BLieDF_Ct and BLieDF_Kt as well as
         ! their constrained pendants gen_CT_lambda and BLieDF_Kt_lambde
         ! have to be referenced to BLieDF_num_Ct and BLieDF_num_Kt resp.
      procedure                  :: BLieDF_print_stats => BLieDF_print_stats
         ! Clean up the BLieDF_problem object (only internal variables
         ! are reset, BLieDF_alpha_m, BLieDF_alpha_f, BLieDF_beta and
         ! BLieDF_gamma are NOT reset)
      procedure                  :: BLieDF_cleanup => BLieDF_cleanup
   end type BLieDF_problem

   ! abstract interface of the problem
   abstract interface
      pure function BLieDF_M(this, q) result(rslt)
         import                                          :: BLieDF_problem
         ! input
         class(BLieDF_problem),   intent(in)               :: this
         real(8), dimension(:), intent(in)               :: q
         ! result
         real(8), dimension(this%sizev,this%sizev)       :: rslt
      end function BLieDF_M

      pure function BLieDF_diag_M(this, q) result(rslt)
         import                                          :: BLieDF_problem
         ! input
         class(BLieDF_problem),   intent(in)               :: this
         real(8), dimension(:), intent(in)               :: q
         ! result
         real(8), dimension(this%sizev)                  :: rslt
      end function BLieDF_diag_M

      pure function BLieDF_g(this, q, v, t) result(rslt)
         import                              :: BLieDF_problem
         ! input
         class(BLieDF_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8), dimension(:), intent(in)   :: v
         real(8),               intent(in)   :: t
         ! result
         real(8), dimension(this%sizev)      :: rslt
      end function BLieDF_g

      pure function BLieDF_qlpexphDqtilde(this, q, h, dq) result(rslt)
         import                              :: BLieDF_problem
         ! input
         class(BLieDF_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8),               intent(in)   :: h
         real(8), dimension(:), intent(in)   :: dq
         ! result
         real(8), dimension(this%sizeq)      :: rslt
      end function BLieDF_qlpexphDqtilde

      pure function BLieDF_itlbtvtw(this, v, w) result(rslt)
         import                              :: BLieDF_problem
         ! input
         class(BLieDF_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         real(8), dimension(:), intent(in)   :: w
         ! result
         real(8), dimension(this%sizev)      :: rslt
      end function BLieDF_itlbtvtw

      pure function BLieDF_Ct(this, q, v, t) result(rslt)
         import                                    :: BLieDF_problem
         ! input
         class(BLieDF_problem),   intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function BLieDF_Ct

      pure function BLieDF_Kt(this, q, v, vd, t) result(rslt)
         import                                    :: BLieDF_problem
         ! input
         class(BLieDF_problem),   intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8), dimension(:), intent(in)         :: vd
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function BLieDF_Kt

      pure function BLieDF_Kt_lambda(this, q, v, vd, l, t) result(rslt)
         import                                    :: BLieDF_problem
         ! input
         class(BLieDF_problem),   intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8), dimension(:), intent(in)         :: vd
         real(8), dimension(:), intent(in)         :: l
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function BLieDF_Kt_lambda

      pure function BLieDF_Tg(this, h, dq) result(rslt)
         import                                    :: BLieDF_problem
         ! input
         class(BLieDF_problem),   intent(in)         :: this
         real(8),               intent(in)         :: h
         real(8), dimension(:), intent(in)         :: dq
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function BLieDF_Tg

      pure function BLieDF_norm(this, v) result(rslt)
         import                  :: BLieDF_problem
         ! input
         class(BLieDF_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         ! result
         real(8)                             :: rslt
      end function BLieDF_norm

      subroutine BLieDF_outputFunction(this,info)
         import                           :: BLieDF_problem
         ! input
         class(BLieDF_problem), intent(in)  :: this
         integer,             intent(in)  :: info
      end subroutine BLieDF_outputFunction

      subroutine BLieDF_init(this)
         import                              :: BLieDF_problem
         ! input/output
         class(BLieDF_problem), intent(inout)  :: this
      end subroutine BLieDF_init

      pure function BLieDF_phi(this,q) result(rslt)
         import                                       :: BLieDF_problem
         ! input
         class(BLieDF_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         ! result
         real(8), dimension(this%sizel)               :: rslt
      end function BLieDF_phi

      pure function BLieDF_B(this,q) result(rslt)
         import                                       :: BLieDF_problem
         ! input
         class(BLieDF_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         ! result
         real(8), dimension(this%sizel,this%sizev)    :: rslt
      end function BLieDF_B

      pure function BLieDF_Z(this,q,v) result(rslt)
         import                                       :: BLieDF_problem
         ! input
         class(BLieDF_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         real(8), dimension(:),            intent(in) :: v
         ! result
         real(8), dimension(this%sizel)               :: rslt
      end function BLieDF_Z

      pure function BLieDF_matZ(this,q,v,T) result(rslt)
         import                                       :: BLieDF_problem
         ! input
         class(BLieDF_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         real(8), dimension(:),            intent(in) :: v
         real(8), dimension(:,:),          intent(in) :: T
         ! result
         real(8), dimension(this%sizel, this%sizev)   :: rslt
      end function BLieDF_matZ

   end interface

   contains

   ! subroutine for pertubing initial values
   subroutine BLieDF_calcInitial(this, h)
      implicit none
      class(BLieDF_problem),  intent(inout)  :: this  ! problem object
      real(8),                intent(in   )  :: h     ! step size
      ! internal variables
      integer                                :: info   ! info flag for dgesv
      integer                                :: ipiv(this%sizev) ! pivot vector for output from dgesv
      real(8)                                :: M(this%sizev,this%sizev)
      !

      associate ( t  => this%t,   &
                  q  => this%q,   &
                  v  => this%v,   &
                  vd => this%vd,  &
                  sv => this%sizev)
         ! calculate $\dot v(t_0)$
         if (this%opts%diag_mass_matrix == 1) then
            if (this%opts%const_mass_matrix == 1) then
               vd = -this%BLieDF_g(q, v, t)/this%BLieDF_const_diag_M
            else
               vd = -this%BLieDF_g(q, v, t)/this%BLieDF_diag_M(q)
            end if
            ! count calls
            this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + 1
         else
            if (this%opts%const_mass_matrix == 1) then
               M = this%BLieDF_const_M
            else
               M = this%BLieDF_M(q)
            end if
            ! we need to solve a linear equation, first calulate the rhs
            vd = -this%BLieDF_g(q, v, t)
            ! count calls
            this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + 1
            ! then solve the system
            call dgesv(         &! solve the System A*X=B and save the result in B
                       sv,      &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                       1,       &! number of right hand sides (=size(B,2))
                       M,       &! matrix A
                       sv,      &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                       ipiv,    &! integer pivot vector; it is not needed
                       vd,      &! matrix B
                       sv,      &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                       info)     ! integer information flag
            ! Now vd actually contains $\dot v(t_0)$
            if (info .ne. 0)  print*, "dgesv sagt info=", info ! TODO
         end if

      end associate
   end subroutine BLieDF_calcInitial

   subroutine BLieDF_calcInitialConstrained(this, h) !TODO TODO
      implicit none
      class(BLieDF_problem),  intent(inout)  :: this  ! problem object
      real(8),                intent(in   )  :: h     ! step size
      ! internal variables
      integer                                                              :: i
      real(8), dimension(this%sizev+this%sizel)                            :: vdl
      integer, dimension(this%sizev+this%sizel)                            :: ipiv0   ! pivot vector for dgesv for LU factorization of MBB0
      integer                                                              :: info    ! info flag for dgesv
      real(8), dimension(this%sizev+this%sizel, this%sizev+this%sizel)     :: MBB0    ! Matrix
      real(8), dimension(this%sizel, this%sizev)                           :: B0      ! $B(q_0)$
      !

      associate ( t  => this%t,            &
                  q  => this%q,            &
                  v  => this%v,            &
                  vd => this%vd,           &
                  l  => this%l,            &
                  sv => this%sizev,        &
                  sl => this%sizel         )
         ! calulate $B(q_0)$
         B0 = this%BLieDF_b(q)
         ! count calls
         this%BLieDF_stats%nBcalls = this%BLieDF_stats%nBcalls + 1

         ! calculate $MBB0$
         if (this%opts%diag_mass_matrix == 1) then
            ! Set the mass matrix part to zero beforehand
            MBB0(1:sv, 1:sv) = 0.0_8
            if (this%opts%const_mass_matrix == 1) then
               MBB0(1:sv, 1) = this%BLieDF_const_diag_M
            else
               MBB0(1:sv, 1) = this%BLieDF_diag_M(q)
            end if
            do concurrent (i=2:sv)
               MBB0(i,i) = MBB0(i,1)
               MBB0(i,1) = 0.0_8
            end do
         else
            if (this%opts%const_mass_matrix == 1) then
               MBB0(1:sv, 1:sv) = this%BLieDF_const_M
            else
               MBB0(1:sv, 1:sv) = this%BLieDF_M(q)
            end if
         end if
         MBB0(sv+1:sv+sl, 1:sv) =           B0
         MBB0(1:sv, sv+1:sv+sl) = transpose(B0)
         MBB0(sv+1:sv+sl, sv+1:sv+sl) = 0.0_8

         ! we need to solve a linear equation, first calulate the rhs
         vdl(1:sv)       = -this%BLieDF_g(q, v, t)
         vdl(sv+1:sv+sl) = -this%BLieDF_Z(q, v)
         ! count calls
         this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + 1
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
   end subroutine BLieDF_calcInitialConstrained

   ! subroutine for integrating one time step
   subroutine BLieDF_solveConstrainedTimeStep(this,t1)
      implicit none
      class(BLieDF_problem), intent(inout)  :: this  ! problem object
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


      ! associate construct for better readability
      associate (alpha  => this%alpha,    &
                 gamma  => this%gamma,    &
                 k      => this%k_bdf,    &
                 sv     => this%sizev,    &
                 sl     => this%sizel,    &
                 v      => this%v,        &
                 vd     => this%vd,       &
                 q      => this%q,        &
                 t      => this%t,        &
                 l      => this%l,        &
                 Dq     => this%Dq,       &
                 Dq_old => this%Dq_old,   &
                 v_old  => this%v_old     )
         ! calculation of step size $h$
         h = t1 - t


         ! initialization of the new values $\dot v_{n+1}$, $a_{n+1}$ and $v_{n+1}$ and the value $\Delta q_n$ and $\lambda_{n+1}$
         if (k==1) then
            Dq = v_old(:,k) ! DEBUG
            v1 = gamma(1)*Dq
         else
            Dq  = Dq_old(:,k-1)
            v1  = matmul(Dq_old, gamma(k:2:-1)) + gamma(1)*Dq ! DEBUG
         end if
         ! Calculate new $v_{n+1}$
         vd1 = (matmul(v_old, alpha(k+1:2:-1)) + alpha(1)*v1)/h ! DEBUG

         ! loop of the newton method
         converged = .false.
         do i=1,this%opts%imax
            ! calculate the new value $q_{n+1}$
            q1  = this%BLieDF_qlpexphDqtilde(q,h,Dq)

            ! caluclate the residue $res$
            res(1:sv) = this%BLieDF_g(q1,v1,t1) + matmul(transpose(this%BLieDF_B(q1)),l)
            ! count calls
            this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + 1
            this%BLieDF_stats%nBcalls = this%BLieDF_stats%nBcalls + 1
            !
            if (this%opts%const_mass_matrix == 1) then
               if (this%opts%diag_mass_matrix == 1) then
                  res(1:sv) = res(1:sv) + this%BLieDF_const_diag_M*vd1
               else
                  res(1:sv) = res(1:sv) + matmul(this%BLieDF_const_M,vd1)
               end if
            else
               if (this%opts%diag_mass_matrix == 1) then
                  res(1:sv) = res(1:sv) + this%BLieDF_diag_M(q1)*vd1
               else
                  res(1:sv) = res(1:sv) + matmul(this%BLieDF_M(q1),vd1)
               end if
            end if
            ! scale
            res(1:sv) = h*res(1:sv)
            res(sv+1:sv+sl) = this%BLieDF_phi(q1)/h

            ! solve the linear System $S_t \cdot \Delta xl = -res$
            Dxl = -res ! the result will be saved in the right hand side's spot

            if (i == 1 .or. this%opts%recalc_iteration_matrix == 1) then
               ! calculate iteration matrix
               if (this%opts%const_mass_matrix == 1) then
                  if (this%opts%diag_mass_matrix == 1) then
                     St(1:sv,1:sv) = 0.0_8
                     forall (i=1:sv)
                        St(    i   ,     i   ) = this%BLieDF_const_diag_M(i) * alpha(1)*gamma(1)
                     end forall
                  else
                     St(   1:sv   ,   1:sv   ) = this%BLieDF_const_M * alpha(1)*gamma(1)
                  end if
               else
                  if (this%opts%diag_mass_matrix == 1) then
                     St(1:sv,1:sv) = 0.0_8
                     St(   1:sv   ,   1   ) = this%BLieDF_diag_M(q1) * alpha(1)*gamma(1)
                     forall (i=2:sv)
                        St(i,i) = St(i,1)
                        St(i,1) = 0.0_8
                     end forall
                  else
                     St(   1:sv   ,   1:sv   ) = this%BLieDF_M(q1) * alpha(1)*gamma(1)
                  end if
               end if
               St(   1:sv   ,sv+1:sv+sl) = transpose(this%BLieDF_B(q1)) ! TODO: dont calculate twice
               St(sv+1:sv+sl,   1:sv   ) = matmul(this%BLieDF_B(q1), this%BLieDF_Tg(h,Dq))
               St(sv+1:sv+sl,sv+1:sv+sl) = 0.0_8
               ! count calls
               this%BLieDF_stats%nBcalls = this%BLieDF_stats%nBcalls + 2

               if (this%opts%no_Ct == 0) then
                  if (this%opts%use_num_Ct == 1) then
                     St(1:sv,1:sv) = St(1:sv,1:sv) + this%BLieDF_num_Ct(q1,v1,t1) * h*gamma(1)
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + this%sizev + 1
                     end if
                  else
                     St(1:sv,1:sv) = St(1:sv,1:sv) + this%BLieDF_Ct(q1,v1,t1) * h*gamma(1)
                  end if
               end if

               if (this%opts%no_Kt == 0) then
                  if (this%opts%use_num_Kt == 1) then
                     St(1:sv,1:sv) = St(1:sv,1:sv) + matmul(this%BLieDF_num_Kt_lambda(q1,v1,vd1,l,t1),this%BLieDF_Tg(h,Dq)) * h*h ! TODO: dgemm oder sowas
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                        this%BLieDF_stats%nBcalls = this%BLieDF_stats%nBcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + this%sizev + 1
                        this%BLieDF_stats%nBcalls = this%BLieDF_stats%nBcalls + this%sizev + 1
                     end if
                  else
                     St(1:sv,1:sv) = St(1:sv,1:sv) + matmul(this%BLieDF_Kt_lambda(q1,v1,vd1,l,t1),this%BLieDF_Tg(h,Dq)) * h*h ! TODO: dgemm oder sowas
                  end if
               end if

               ! Actually solve the linear System $S_t \cdot \Delta xl = -res$
               if (this%opts%banded_iteration_matrix == 1) then
#ifdef DEBUG_PRINT_ITERATION_MATRIX_AT
                  if (t >= DEBUG_PRINT_ITERATION_MATRIX_AT) then
                     call print_matrix(St,'St')
                     call print_vector_int(this%opts%jour,'jour')
                     stop "Successfully printed iteration matrix"
                  end if
#endif
                  ! TODO
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
            Dq  = Dq  + Dxl(1:sv)                       REL
            v1  = v1  + gamma(1) * Dxl(1:sv)            REL
            vd1 = vd1 + alpha(1)*gamma(1)/h * Dxl(1:sv) REL
            l   = l   + Dxl(sv+1:sv+sl)/h               REL
#undef REL

            ! update solver stats
            this%BLieDF_stats%newt_steps_curr = i

            ! exit loop, if converged
            if ( converged ) exit

         end do ! end newton method

         ! DEBUG
         if (this%BLieDF_stats%newt_steps_curr >= this%opts%imax) then
            print *, ' ' ! DEBUG
            print *, "Newton did not converge at t=", t
            !! DEBUG !
            call this%BLieDF_print_stats
            print *, "Iteration matrix:"
            call print_matrix(St,'St')
            errorstop
            !! GUBED !
         end if


         ! The integration step is done; now the new values have to be
         ! applied to the problem object, and the new value for a can
         ! be calculated
         t  = t1
         q  = q1
         v  = v1
         vd = vd1

      end associate
   end subroutine BLieDF_solveConstrainedTimeStep

   ! subroutine for integrating one time step
   subroutine BLieDF_solveConstrainedTimeStep_stab2(this,t1)
      implicit none
      class(BLieDF_problem), intent(inout)  :: this  ! problem object
      real(8)            , intent(in   )  :: t1    ! $t_{n+1}$

      ! internal integer variables
      integer                                      :: i     ! for iteration
      integer, dimension(this%sizev+2*this%sizel)  :: ipiv  ! needed for dgesv
      integer                                      :: info  ! needed for dgesv

      ! intenal real variables
      real(8)                                                              :: h       ! step size
      real(8), dimension(this%sizev)                                       :: vd1     ! $\dot v_{n+1}$
      real(8), dimension(this%sizev)                                       :: v1      ! $v_{n+1}$
      real(8), dimension(this%sizev)                                       :: Dq      ! $\Delta q_n$
      real(8), dimension(this%sizeq)                                       :: q1      ! $q_{n+1}}
      real(8), dimension(this%sizev+2*this%sizel)                          :: res     ! $res(\dots)$
      real(8), dimension(this%sizev+2*this%sizel,this%sizev+2*this%sizel)  :: St      ! $S_t$
      real(8), dimension(this%sizev,this%sizev)                            :: Mstar   ! $M^\ast$
      real(8), dimension(this%sizel,this%sizev)                            :: Bqn     ! $B(q_n)$
      real(8), dimension(this%sizel,this%sizev)                            :: B       ! $B$
      real(8), dimension(:,:), allocatable                                 :: Stb     ! permuted $S_t$ in banded format
      real(8), dimension(this%sizev+2*this%sizel)                          :: Dxleta  ! $(\Delta x, \Delta \lambda, \Delta \eta)^T$

      ! internal logical
      logical                                                              :: converged

      ! associate construct for better readability
      associate (alpha  => this%alpha,    &
                 gamma  => this%gamma,    &
                 k      => this%k_bdf,    &
                 sv     => this%sizev,    &
                 sl     => this%sizel,    &
                 v      => this%v,        &
                 vd     => this%vd,       &
                 q      => this%q,        &
                 t      => this%t,        &
                 l      => this%l,        &
                 Dq     => this%Dq,       &
                 Dq_old => this%Dq_old,   &
                 v_old  => this%v_old,    &
                 eta    => this%eta       )
         ! calculation of step size $h$
         h = t1 - t

         ! calculate B(q_n)
         Bqn = this%BLieDF_B(q)
         ! count calls
         this%BLieDF_stats%nBcalls = this%BLieDF_stats%nBcalls + 1

         ! initialization
         if (k==1) then
            Dq = v_old(:,k) ! DEBUG
            v1 = gamma(1)*Dq
         else
            Dq  = Dq_old(:,k-1) - matmul(transpose(Bqn),eta)
            v1 = matmul(Dq_old, gamma(k:2:-1)) + gamma(1)*Dq ! DEBUG
         end if
         vd1 = (matmul(v_old, alpha(k+1:2:-1)) + alpha(1)*v1)/h

         ! initialization of and $\lambda_{n+1}$, $\eta_{n+1}
         !l   = l    ! this is not necessary to do, obviously
         !eta = eta  ! dito

         ! loop of the newton method
         converged = .false.
         do i=1,this%opts%imax
            ! calculate the new value $q_{n+1}$
            q1  = this%BLieDF_qlpexphDqtilde(q,h,Dq)

            ! calculate new B(q1)
            B = this%BLieDF_B(q1)
            ! count calls
            this%BLieDF_stats%nBcalls = this%BLieDF_stats%nBcalls + 1

            ! caluclate the residue $res_h$
            res(1:sv) = this%BLieDF_g(q1,v1,t1) + matmul(transpose(B),l)
            ! count calls
            this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + 1
            !
            if (this%opts%const_mass_matrix == 1) then
               if (this%opts%diag_mass_matrix == 1) then
                  res(1:sv) = res(1:sv) + this%BLieDF_const_diag_M*vd1
               else
                  res(1:sv) = res(1:sv) + matmul(this%BLieDF_const_M,vd1)
               end if
            else
               if (this%opts%diag_mass_matrix == 1) then
                  res(1:sv) = res(1:sv) + this%BLieDF_diag_M(q1)*vd1
               else
                  res(1:sv) = res(1:sv) + matmul(this%BLieDF_M(q1),vd1)
               end if
            end if
            res(1:sv) = res(1:sv) * h
            res(sv+1:sv+sl) = this%BLieDF_phi(q1) / h
            res(sv+sl+1:sv+sl+sl) = matmul(B, v1)

            ! solve the linear System $S_t \cdot \Delta xl = -res$
            Dxleta = -res ! the result will be saved in the right hand side's spot

            if (i == 1 .or. this%opts%recalc_iteration_matrix == 1) then
               ! calculate iteration matrix
               if (this%opts%const_mass_matrix == 1) then
                  if (this%opts%diag_mass_matrix == 1) then
                     Mstar = 0.0_8
                     forall (i=1:sv)
                        Mstar(i, i) = this%BLieDF_const_diag_M(i)
                     end forall
                  else
                     Mstar = this%BLieDF_const_M
                  end if
               else
                  if (this%opts%diag_mass_matrix == 1) then
                     Mstar = 0.0_8
                     Mstar(1:sv, 1) = this%BLieDF_diag_M(q1)
                     forall (i=2:sv)
                        Mstar(i,i) = Mstar(i,1)
                        Mstar(i,1) = 0.0_8
                     end forall
                  else
                     Mstar = this%BLieDF_M(q1)
                  end if
               end if
               Mstar = Mstar*alpha(1)*gamma(1)

               if (this%opts%no_Ct == 0) then
                  ! Add damping matrix to Mstar
                  if (this%opts%use_num_Ct == 1) then
                     Mstar = Mstar + h*gamma(1) * this%BLieDF_num_Ct(q1,v1,t1)
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + this%sizev + 1
                     end if
                  else
                     Mstar = Mstar + h*gamma(1) * this%BLieDF_Ct(q1,v1,t1)
                  end if
               end if

               St(      1:sv,             1:sv      ) = Mstar
               St(      1:sv      ,    sv+1:sv+sl   ) = transpose(B)
               St(sv   +1:sv+sl   ,       1:sv      ) = matmul(B, this%BLieDF_Tg(h,Dq))
               St(sv   +1:sv+sl   ,    sv+1:sv+sl   ) = 0.0_8
#ifdef NOZ
               St(sv+sl+1:sv+sl+sl,       1:sv      ) = gamma(1)*B
#else
               St(sv+sl+1:sv+sl+sl,       1:sv      ) = gamma(1)*B + h*this%BLieDF_matZ(q1,v1,this%BLieDF_Tg(h,Dq))
#endif
               St(sv+sl+1:sv+sl+sl, sv   +1:sv+sl   ) = 0.0_8
               St(sv+sl+1:sv+sl+sl, sv+sl+1:sv+sl+sl) = gamma(1)*matmul(B,transpose(Bqn))
               St(sv   +1:sv+sl   , sv+sl+1:sv+sl+sl) = 0.0_8
               St(      1:sv      , sv+sl+1:sv+sl+sl) = matmul(Mstar,transpose(Bqn))

               if (this%opts%no_Kt == 0) then
                  if (this%opts%use_num_Kt == 1) then
                     St(1:sv,1:sv) = St(1:sv,1:sv) + h*h*matmul(this%BLieDF_num_Kt_lambda(q1,v1,vd1,l,t1),this%BLieDF_Tg(h,Dq)) ! TODO: dgemm oder sowas
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                        this%BLieDF_stats%nBcalls = this%BLieDF_stats%nBcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + this%sizev + 1
                        this%BLieDF_stats%nBcalls = this%BLieDF_stats%nBcalls + this%sizev + 1
                     end if
                  else
                     St(1:sv,1:sv) = St(1:sv,1:sv) + h*h*matmul(this%BLieDF_Kt_lambda(q1,v1,vd1,l,t1),this%BLieDF_Tg(h,Dq)) ! TODO: dgemm oder sowas
                  end if
               end if

               ! Actually solve the linear System $S_t \cdot \Delta xleta = -res$
               if (this%opts%banded_iteration_matrix == 1) then
#ifdef DEBUG_PRINT_ITERATION_MATRIX_AT
                  if (t >= DEBUG_PRINT_ITERATION_MATRIX_AT) then
                     call print_matrix(St,'St')
                     call print_vector_int(this%opts%jour,'jour')
                     stop "Successfully printed iteration matrix"
                  end if
#endif
                  ! TODO
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag)
                     if (.not. allocated(Stb)) then
                        allocate(Stb(lo*2+hi+1, sv+2*sl))
                     end if
                     !Stb = 0.0_8 ! DEBUG
                     forall (i=1:sv+2*sl)
                        Stb((max(1, i-hi)-i+lo+1+hi):(min(sv+sl,i+lo)-i+lo+1+hi),i) &
                            = St(this%opts%jour(max(1, i-hi):min(sv+2*sl, i + lo)), this%opts%jour(i))
                     end forall
                     Dxleta = Dxleta(this%opts%jour)
                     call dgbsv(           &! solve the system A*X=B and save the result in B
                                 sv+2*sl,  &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,       &!
                                 hi,       &!
                                 1,        &! number of right hand sides (=size(B,2))
                                 Stb,      &! matrix A
                                 2*lo+hi+1,&! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv,     &! integer pivot vector; it is not needed
                                 Dxleta,   &! matrix B
                                 sv+2*sl,  &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)      ! integer information flag
                     Dxleta(this%opts%jour) = Dxleta
                  end associate
               else
                  call dgesv(          &! solve the System A*X=B and save the result in B
                              sv+sl+sl,&! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                              1,       &! number of right hand sides (=size(B,2))
                              St,      &! matrix A
                              sv+sl+sl,&! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv,    &! integer pivot vector; it is not needed
                              Dxleta,  &! matrix B
                              sv+sl+sl,&! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                              info)     ! integer information flag
               end if

               if (info.ne.0) print *, "dgesv sagt info=", info

               ! update $\dot v_{n+1}$, $v_{n+1}$ and $\Delta q_n$
            else ! (St was not (re)calculated)
               ! Actually solve the linear System $S_t \cdot \Delta xleta = -res$
               if (this%opts%banded_iteration_matrix == 1) then
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag, &
                             uStb => ubound(Stb,1))
                     Dxleta = Dxleta(this%opts%jour) !TODO
                     call dgbtrs(               &
                                 'No transpose',&! solve the System A*X=B and save the result in B
                                 sv+2*sl,       &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,            &!
                                 hi,            &!
                                 1,             &! number of right hand sides (=size(B,2))
                                 Stb,           &! matrix A
                                 2*lo+hi+1,     &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv,          &! integer pivot vector; it is not needed
                                 Dxleta,        &! matrix B
                                 sv+2*sl,       &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)           ! integer information flag
                     Dxleta(this%opts%jour) = Dxleta
                  end associate
               else
                  call dgetrs(                &
                              'No transpose', & ! solve the System A*X=B and save the result in B
                              sv + sl + sl,   & ! number of linear equations (=size(A,1))  ! Vorsicht: double precision muss real(8) sein
                              1,              & ! number of right hand sides (=size(B,2))
                              St,             & ! matrix A
                              sv + sl + sl,   & ! leading dimesnion of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv,           & ! integer pivot vector; holds the pivoting that was calculated in the dgesv call
                              Dxleta,         & ! matrix B
                              sv + sl + sl,   & ! leading dimension of B, in thie case is equal to the number of linear equations (=size(A,1))
                              info)             ! integer information flag
               end if
               ! DEBUG
               if (info.ne.0) print*, 'dgetrs sagt info=', info
               ! GUBED
            end if

            ! Check for convergence
            if ( ( sum((Dxleta(      1:      sv) / (this%opts%atol + this%opts%rtol*abs(Dq )))**2 )  &
                  +sum((Dxleta(sv+   1:   sv+sl) / (this%opts%atol + this%opts%rtol*abs(h*l)))**2 )  &
                  +sum((Dxleta(sv+sl+1:sv+sl+sl) / (this%opts%atol + this%opts%rtol*abs(eta)))**2 ) )&
                   / (sv+2*sl)   <=  1.0_8 ) then
               converged = .true.
            end if

            ! update $\dot v_{n+1}$, $v_{n+1}$, $\Delta q_n$ and $\lambda_{n+1}$, $\eta_{n+1}$
! DEBUGGETY
#define REL *1
            Dq  = Dq  + Dxleta(1:sv)                                                                             REL
            v1  = v1  + gamma(1)*(Dxleta(1:sv) + matmul(transpose(Bqn),Dxleta(sv+sl+1:sv+sl+sl)))                REL
            vd1 = vd1 + alpha(1)*gamma(1)/h * (Dxleta(1:sv) + matmul(transpose(Bqn),Dxleta(sv+sl+1:sv+sl+sl)))   REL
            l   = l   + Dxleta(sv   +1:sv+sl  )/h                                                                REL
            eta = eta + Dxleta(sv+sl+1:sv+sl+sl)                                                                 REL
#undef REL

            ! update solver stats
            this%BLieDF_stats%newt_steps_curr = i

            ! exit loop, if converged
            if ( converged ) exit

         end do ! end newton method

         ! DEBUG
         if (this%BLieDF_stats%newt_steps_curr >= this%opts%imax) then
            print *, ' ' ! DEBUG
            print *, "Newton did not converge at t=", t
            !! DEBUG !
            call this%BLieDF_print_stats
            print *, "Iteration matrix:"
            !call print_matrix(St,'St')
            errorstop
            !! GUBED !
         end if

         ! The integration step is done; now the new values have to be
         ! applied to the problem object, and the new value for a can
         ! be calculated
         t  = t1
         q  = q1
         v  = v1
         vd = vd1

      end associate
   end subroutine BLieDF_solveConstrainedTimeStep_stab2

   ! subroutine for integrating one time step
   subroutine BLieDF_solveTimeStep(this,t1)
      implicit none
      class(BLieDF_problem), intent(inout)      :: this  ! problem object
      real(8),               intent(in   )      :: t1    ! $t_{n+1}$

      ! internal integer variables
      integer                                   :: i     ! for iteration
      integer, dimension(this%sizev)            :: ipiv  ! needed for dgesv
      integer                                   :: info  ! needed for dgesv

      ! internal real variables
      real(8)                                   :: h       ! step size
      real(8), dimension(this%sizev)            :: vd1     ! $\dot v_{n+1}$
      real(8), dimension(this%sizev)            :: v1      ! $v_{n+1}$
      real(8), dimension(this%sizeq)            :: q1      ! $q_{n+1}}
      real(8), dimension(this%sizev)            :: res     ! $res$
      real(8), dimension(this%sizev,this%sizev) :: St      ! $S_t$
      real(8), dimension(:,:), allocatable      :: Stb     ! $S_t$ in banded format
      real(8), dimension(this%sizev)            :: Dx      ! $\Delta x$

      ! internal logical
      logical                                   :: converged


      ! associate construct for better readability
      associate (alpha  => this%alpha,           &
                 gamma  => this%gamma,           &
                 k      => this%k_bdf,           &
                 sv     => this%sizev,           &
                 v      => this%v,               &
                 vd     => this%vd,              &
                 q      => this%q,               &
                 t      => this%t,               &
                 Dq     => this%Dq,              &
                 Dq_old => this%Dq_old,          &
                 v_old  => this%v_old            )
         ! calculation of step size $h$
         h = t1 - t

         ! initialization TODO! DEBUG!
         if (k==1) then
            Dq = v_old(:,k) ! DEBUG
            v1 = gamma(1)*Dq
         else
            Dq  = Dq_old(:,k-1)
            v1  = matmul(Dq_old, gamma(k:2:-1)) + gamma(1)*Dq ! DEBUG
         end if
         ! calculate new value $v_{n+1}$
         ! calculate new value $\dot v_{n+1}$
         vd1 = (matmul(v_old, alpha(k+1:2:-1)) + alpha(1)*v1)/h ! DEBUG

         ! loop of the newton method
         converged = .false.
         do i=1,this%opts%imax
            ! calculate the new value $q_{n+1}$
            q1  = this%BLieDF_qlpexphDqtilde(q,h,Dq)
            ! caluclate the residue $res$
            if (this%opts%const_mass_matrix == 1) then
               if (this%opts%diag_mass_matrix == 1) then
                  res = this%BLieDF_const_diag_M*vd1 + this%BLieDF_g(q1,v1,t1)
               else
                  res = matmul(this%BLieDF_const_M,vd1) + this%BLieDF_g(q1,v1,t1)
               end if
            else
               if (this%opts%diag_mass_matrix == 1) then
                  res = this%BLieDF_diag_M(q1)*vd1 + this%BLieDF_g(q1,v1,t1)
               else
                  res = matmul(this%BLieDF_M(q1),vd1) + this%BLieDF_g(q1,v1,t1)
               end if
            end if
            ! scale
            res = h*res
            ! count calls
            this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + 1


            Dx = -res ! the result will be saved in the right hand side's spot

            if (i == 1  .or. this%opts%recalc_iteration_matrix == 1) then
               ! calculate iteration matrix
               if (this%opts%const_mass_matrix == 1) then
                  if (this%opts%diag_mass_matrix == 1) then
                     St = 0.0_8
                     forall (i=1:sv)
                        St(i,i) = this%BLieDF_const_diag_M(i) * alpha(1)*gamma(1)
                     end forall
                  else
                     St = this%BLieDF_const_M * alpha(1)*gamma(1)
                  end if
               else
                  if (this%opts%diag_mass_matrix == 1) then
                     St = 0.0_8
                     St(:,1) = this%BLieDF_diag_M(q1) * alpha(1)*gamma(1)
                     forall (i=2:sv)
                        St(i,i) = St(i,1)
                        St(i,1) = 0.0_8
                     end forall
                  else
                     St = this%BLieDF_M(q1) * alpha(1)*gamma(1)
                  end if
               end if

               if (this%opts%no_Ct == 0) then
                  if (this%opts%use_num_Ct == 1) then
                     St = St + this%BLieDF_num_Ct(q1,v1,t1) * h*gamma(1)
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + this%sizev + 1
                     end if
                  else
                     St = St + this%BLieDF_Ct(q1,v1,t1) * h*gamma(1)
                  end if
               end if

               if (this%opts%no_Kt == 0) then
                  if (this%opts%use_num_Kt == 1) then
                     St = St + matmul(this%BLieDF_num_Kt(q1,v1,vd1,t1),this%BLieDF_Tg(h,Dq)) * h*h! TODO: dgemm oder sowas: Idee: Man kann per compilerflags einschalten, dass für matmul direkt eine blas-routine benutzt wird
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%BLieDF_stats%ngcalls = this%BLieDF_stats%ngcalls + this%sizev + 1
                     end if
                  else
                     St = St + matmul(this%BLieDF_Kt(q1,v1,vd1,t1),this%BLieDF_Tg(h,Dq)) * h*h ! TODO: dgemm oder sowas
                  end if
               end if

               !! DEBUG !
               !print *, 'i=', i
               !print *, 't=', t
               !call print_vector(q,'q')
               !call print_vector(v,'v')
               !call print_vector(vd,'vd')
               !call print_vector(a,'a')
               !!call print_matrix(St,'Stn')
               !call print_vector(res, 'res')
               !!call print_vector(this%BLieDF_const_diag_M,'Mn')
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
            v1  = v1  + gamma(1) * Dx            REL
            vd1 = vd1 + alpha(1)*gamma(1)/h * Dx REL
#undef REL
! GUBED!!


            ! update solver stats
            this%BLieDF_stats%newt_steps_curr = i

            ! exit loop, if converged
            if ( converged ) exit

         end do ! end newton method

         !! DEBUG
         if (this%BLieDF_stats%newt_steps_curr >= this%opts%imax) then
            print *, ' ' ! DEBUG
            print *, "Newton did not converge at t=", t
            !!! DEBUG !
            !call this%BLieDF_print_stats
            !print *, "Iteration matrix:"
            !call print_matrix(St,'St')
            errorstop
            !!! GUBED !
         end if


         ! The integration step is done; now the new values have to be
         ! applied to the problem object, and the new value for a can
         ! be calculated
         t  = t1
         q  = q1
         v  = v1
         vd = vd1


      end associate

   end subroutine BLieDF_solveTimeStep


   pure function cum_sum_omit_last(vec) result(rslt)
      ! input
      real(8), intent(in)  :: vec(:)
      ! result
      real(8)              :: rslt(size(vec)-1)
      ! internal
      integer              :: i
      !
      rslt(1) = vec(1)
      do i=2,size(vec)-1
         rslt(i) = rslt(i-1) + vec(i)
      end do
   end function cum_sum_omit_last


   ! subroutine to integrate
   subroutine BLieDF_integrate(this)
      implicit none
      ! input/output
      class(BLieDF_problem),  intent(inout)  :: this
      ! internal
      integer  :: n  ! needed for iteration
      real(8)  :: h  ! step size $h$
      real(8)  :: t1 ! next time $t_{n+1}$
      real(4)  :: times(2), time ! for dtime
      integer  :: the_k_bdf
      real(8), allocatable :: the_Dq_old(:,:)
      real(8), allocatable :: the_v_old(:,:)

      ! problem initialization
      call this%BLieDF_init()

      ! TODO: More error checking
      if (this%opts%constrained == 1             .and. &
          this%opts%banded_iteration_matrix == 1 .and. &
          .not. allocated(this%opts%jour)               ) then
         print *, '[ERROR] opts%jour not allocated, aborting'
         errorstop
      end if
      ! TODO: t0 und t vergleichen (müssen gleich sein)

      ! initialize output function
      call this%BLieDF_outputFunction(0)

      ! Calculate step size $h$
      h = (this%opts%te - this%opts%t0)/this%opts%nsteps

      ! Set stats of solver to zero
      this%BLieDF_stats%newt_steps_curr = 0
      this%BLieDF_stats%newt_steps_sum = 0
      this%BLieDF_stats%newt_steps_max = 0
      this%BLieDF_stats%newt_steps_avg = 0
      this%BLieDF_stats%ngcalls = 0
      this%BLieDF_stats%nBcalls = 0

      ! if mass matrix is constant, calculate it
      if (this%opts%const_mass_matrix == 1) then
         if (this%opts%diag_mass_matrix == 1) then
            if (.not. allocated(this%BLieDF_const_diag_M)) then
               allocate(this%BLieDF_const_diag_M(this%sizev))
            end if
            this%BLieDF_const_diag_M = this%BLieDF_diag_M(this%q)
         else
            if (.not. allocated(this%BLieDF_const_M)) then
               allocate(this%BLieDF_const_M(this%sizev,this%sizev))
            end if
            this%BLieDF_const_M = this%BLieDF_M(this%q)
         end if
      end if

      ! start stopwatch
      time = dtime(times)

      ! calculate initial values
      if (this%opts%constrained == 1) then
         call this%BLieDF_calcInitialConstrained(h)
      else
         call this%BLieDF_calcInitial(h)
      end if

      ! allocate space for Dq, the_Dq_old and the_v_old
      allocate(this%Dq(this%sizev))
      allocate(the_Dq_old(this%sizev,this%k_bdf-1))
      allocate(the_v_old(this%sizev,this%k_bdf))

      ! output for the first time
      call this%BLieDF_outputFunction(1)

      ! Start the integrator by doing steps with the lower-order methods
      ! increasin the k in each step until enough data is aquired
      the_k_bdf = this%k_bdf
      the_v_old(:,1) = this%v
      do n=1,the_k_bdf-1
         ! Calculate the next time $t_{n+1}$
         t1 = this%opts%t0 + n*h

         ! At this point, we only have information about $n$ previous steps,
         ! so we use the $n$-BLieDF method
         this%k_bdf = n

         ! Allocate and set the method parameters
         allocate(this%alpha(n+1))
         this%alpha = alphas(1:n+1,n)
         !
         allocate(this%gamma(n))
         this%gamma = cum_sum_omit_last(this%alpha)

         ! Allocate and fill $\delta q$ informations from the previous steps
         allocate(this%Dq_old(this%sizev,n-1))
         this%Dq_old(:,1:n-1) = the_Dq_old(:,1:n-1)

         ! Allocate and fill $v$ informations from the previous steps
         allocate(this%v_old(this%sizev,n))
         this%v_old(:,1:n) = the_v_old(:,1:n)

         ! solve time step
         if (this%opts%constrained == 1) then
            if (this%opts%stab2 == 1) then
               call this%BLieDF_solveConstrainedTimeStep_stab2(t1)
            else
               call this%BLieDF_solveConstrainedTimeStep(t1)
            end if
         else
            call this%BLieDF_solveTimeStep(t1)
         end if
         ! update solver stats
         this%BLieDF_stats%newt_steps_sum = this%BLieDF_stats%newt_steps_sum &
            + this%BLieDF_stats%newt_steps_curr
         this%BLieDF_stats%newt_steps_avg = real(this%BLieDF_stats%newt_steps_sum,8)/n
         this%BLieDF_stats%newt_steps_max = max( &
            this%BLieDF_stats%newt_steps_max,    &
            this%BLieDF_stats%newt_steps_curr  )
         ! output normally
         call this%BLieDF_outputFunction(1)

         ! Fill in new information
         the_Dq_old(:,n) = this%Dq
         the_v_old(:,n+1) = this%v

         ! Deallocate alpha, gamma, Dq_old and v_old, because they need to
         ! be bigger in the next step
         deallocate(this%alpha)
         deallocate(this%gamma)
         deallocate(this%Dq_old)
         deallocate(this%v_old)
      end do

      ! Set up final integrator
      this%k_bdf = the_k_bdf

      ! Allocate and set the method parameters
      allocate(this%alpha(this%k_bdf+1))
      this%alpha = alphas(1:this%k_bdf+1, this%k_bdf)
      !
      allocate(this%gamma(this%k_bdf))
      this%gamma = cum_sum_omit_last(this%alpha)

      ! Allocate and fill $\delta q$ informations from the previous steps
      allocate(this%Dq_old(this%sizev,this%k_bdf-1))
      this%Dq_old = the_Dq_old

      ! Allocate and fill $v$ informations from the previous steps
      allocate(this%v_old(this%sizev,this%k_bdf))
      this%v_old = the_v_old

      ! Deallocate the_Dq_old and the_v_old
      deallocate(the_Dq_old)
      deallocate(the_v_old)

      ! integration loop
      do n=the_k_bdf,this%opts%nsteps
         ! Calculate the next time $t_{n+1}$
         t1 = this%opts%t0 + n*h
         ! solve time step
         if (this%opts%constrained == 1) then
            if (this%opts%stab2 == 1) then
               call this%BLieDF_solveConstrainedTimeStep_stab2(t1)
            else
               call this%BLieDF_solveConstrainedTimeStep(t1)
            end if
         else
            call this%BLieDF_solveTimeStep(t1)
         end if
         ! update solver stats
         this%BLieDF_stats%newt_steps_sum = this%BLieDF_stats%newt_steps_sum &
            + this%BLieDF_stats%newt_steps_curr
         this%BLieDF_stats%newt_steps_avg = real(this%BLieDF_stats%newt_steps_sum,8)/n
         this%BLieDF_stats%newt_steps_max = max( &
            this%BLieDF_stats%newt_steps_max,    &
            this%BLieDF_stats%newt_steps_curr  )
         ! output normally
         call this%BLieDF_outputFunction(1)
         ! update Dq_old and v_old
         if (.not. this%k_bdf == 1) then
            this%Dq_old(:,1:this%k_bdf-2) = this%Dq_old(:,2:this%k_bdf-1)
            this%Dq_old(:,this%k_bdf-1)   = this%Dq
         end if
         this%v_old(:,1:this%k_bdf-1)  = this%v_old(:,2:this%k_bdf)
         this%v_old(:,this%k_bdf)      = this%v


      end do

      ! stop stopwatch
      this%BLieDF_stats%time = dtime(times) - time

      ! output to terminate
      call this%BLieDF_outputFunction(99)

   end subroutine BLieDF_integrate

   pure function BLieDF_num_Ct(this, q, v, t) result(rslt)
      ! input
      class(BLieDF_problem),   intent(in)         :: this
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
      ! Ct is the derivative of r wrt v. But v only appears in g, so
      ! it is sufficient to calculate the derivative of g.
      !
      ! Calculate g ! TODO: man könnte g0 durchschleifen oder so
      g0 = this%BLieDF_g(q,v,t)
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
            rslt(:,i) = this%BLieDF_g(q,v+w,t) - g0
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
            rslt(:,i) = (this%BLieDF_g(q,v+w,t) - g0)/w(i)
         end do
      end if
   end function BLieDF_num_Ct

   pure function BLieDF_num_Kt(this, q, v, vd, t) result(rslt)
      ! input
      class(BLieDF_problem),   intent(in)         :: this
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
      ! Kt is the derivative of r wrt q. WE ASSUME THAT THE MASS MATRIX
      ! M DOES NOT DEPEND ON Q, OTHERWISE THIS APPROXIMATION WILL BE
      ! FALSE! TODO! TODO!
      ! Thus, it is sufficient to calculate the derivative of g.
      !
      ! Calculate g ! TODO: man könnte g0 durchschleifen oder so
      g0 = this%BLieDF_g(q,v,t)
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
            rslt(:,i) = this%BLieDF_g(this%BLieDF_qlpexphDqtilde(q,1.0_8, w),v,t) - g0
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
            rslt(:,i) = (this%BLieDF_g(this%BLieDF_qlpexphDqtilde(q,1.0_8, w),v,t) - g0) / w(i)
            !w    = 0.0_8
            !w(i) = h
            !rslt(:,i) = (this%BLieDF_g(q+this%BLieDF_tilde(w),v,t) - g0)/h
         end do
      end if
   end function BLieDF_num_Kt

! DEBUG
   pure function BLieDF_num_B(this, q) result(rslt)
      ! input
      class(BLieDF_problem),   intent(in)         :: this
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
      Phi0 = this%BLieDF_Phi(q)

      ! set w to zero
      w = 0.0_8
      ! loop over the columns of Ct
      do i=1,this%sizev
         if (.not. i == 1) w(i-1) = 0.0_8
         w(i) = max( abs(q(i))*1.0e-8_8,  1.0e-12_8)
         rslt(:,i) = (this%BLieDF_Phi(this%BLieDF_qlpexphDqtilde(q,1.0_8, w)) - Phi0) / w(i)
         !w    = 0.0_8
         !w(i) = h
         !rslt(:,i) = (this%BLieDF_g(q+this%BLieDF_tilde(w),v,t) - g0)/h
      end do
   end function BLieDF_num_B
! GUBED

   pure function BLieDF_num_Kt_lambda(this, q, v, vd, l, t) result(rslt)
      ! input
      class(BLieDF_problem),            intent(in)   :: this
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
      ! Kt is the derivative of r wrt q. WE ASSUME THAT THE MASS MATRIX
      ! M DOES NOT DEPEND ON Q, OTHERWISE THIS APPROXIMATION WILL BE
      ! FALSE! TODO! TODO!
      ! Thus, it is sufficient to calculate the derivative of g.
      !
      ! Calculate g ! TODO: man könnte g0 durchschleifen oder so
      g0 = this%BLieDF_g(q,v,t) + matmul(transpose(this%BLieDF_B(q)),l)
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
            rslt(:,i) = this%BLieDF_g(this%BLieDF_qlpexphDqtilde(q,1.0_8, w),v,t) &
               + matmul(transpose(this%BLieDF_B(this%BLieDF_qlpexphDqtilde(q,1.0_8,w))),l)- g0
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
            !rslt(:,i) = (this%BLieDF_g(q+this%BLieDF_tilde(w),v,t) + matmul(transpose(this%BLieDF_B(q+this%BLieDF_tilde(w))),l) - g0)/h
            !w(i) = 1.0_8 ! TODO
            rslt(:,i) = (this%BLieDF_g(this%BLieDF_qlpexphDqtilde(q,1.0_8,w),v,t)  &
             + matmul(transpose(this%BLieDF_B(this%BLieDF_qlpexphDqtilde(q,1.0_8,w))),l) - g0)/w(i)
         end do
      end if
   end function BLieDF_num_Kt_lambda

   subroutine BLieDF_print_stats(this)
      ! input
      class(BLieDF_problem), intent(in)  :: this
      !
      print *, 'time:          ', this%BLieDF_stats%time
      print *, '#calls of g:   ', this%BLieDF_stats%ngcalls
      print *, '#calls of B:   ', this%BLieDF_stats%nBcalls
      print *, 'newt_steps_max:', this%BLieDF_stats%newt_steps_max
      print *, 'newt_steps_avg:', this%BLieDF_stats%newt_steps_avg
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

   subroutine BLieDF_cleanup(this)
      implicit none
      ! input/output
      class(BLieDF_problem), intent(inout)  :: this
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
      this%BLieDF_stats%newt_steps_curr = 0
      this%BLieDF_stats%newt_steps_max = 0
      this%BLieDF_stats%newt_steps_avg = 0.0_8
      this%BLieDF_stats%ngcalls = 0
      this%BLieDF_stats%nBcalls = 0
      this%BLieDF_stats%time = 0.0_8
      !
      this%t = 0.0_8
      this%sizeq = 0
      this%sizev = 0
      this%sizel = 0
      !
      if (allocated(this%alpha))  deallocate(this%alpha)
      if (allocated(this%gamma))  deallocate(this%gamma)
      if (allocated(this%Dq_old))   deallocate(this%Dq_old)
      if (allocated(this%v_old))    deallocate(this%v_old)
      if (allocated(this%q))  deallocate(this%q)
      if (allocated(this%v))  deallocate(this%v)
      if (allocated(this%vd)) deallocate(this%vd)
      if (allocated(this%Dq)) deallocate(this%Dq)
      if (allocated(this%l))  deallocate(this%l)
      if (allocated(this%eta)) deallocate(this%eta)
      if (allocated(this%BLieDF_const_M)) deallocate(this%BLieDF_const_M)
      if (allocated(this%BLieDF_const_diag_M)) deallocate(this%BLieDF_const_diag_M)
      if (allocated(this%opts%jour)) deallocate(this%opts%jour)
   end subroutine BLieDF_cleanup

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

end module BLieDF
