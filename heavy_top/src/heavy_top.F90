#define GL(x) x

! describe the heavy top in different Lie groups

module heavy_top
   use GL(INTEGRATOR)
   use lie_group_functions
   use s3sdr3_functions
   implicit none



   ! extension of the type GL(INTEGRATOR)_problem
   type, extends(GL(INTEGRATOR)_problem)   :: heavy_top_t
      !
      ! Output: file name and luns
      character(len=256)      :: out_fname
      integer                 :: out_bin_lun
      integer                 :: out_misc_lun
      !
      ! Which Lie group formulation to use: (See config.lua)
      integer                 :: liegroup
      !
      ! Mass of the top
      real(8)                 :: mass
      !
      ! Diagonal elements of the inertial tensor (wrt. point of origin or
      !    center of mass of the top; depending on the formulation)
      real(8)                 :: inerJ(3)
      !
      ! Gravity including its direction
      real(8)                 :: gravity(3)
      !
      ! Initial configuration
      real(8)                 :: x0(3)
      real(8)                 :: p0(4)
      real(8)                 :: Om0(3)
      !
      ! Reference position (i.e. position where the top is at when it is not rotated)
      real(8)                 :: refX(3)
      !
      ! How to output the calculated data (see config.lua)
      integer                 :: output_type
      ! Output in t direction
      integer                 :: output_t_at
      real(8)                 :: t_output_at_multiples_of
   contains

      ! referencing the former deferred procedures
      procedure   :: GL(INTEGRATOR)_M              => M
      procedure   :: GL(INTEGRATOR)_diag_M         => diagM
#if !defined(INT_RATTLie) && !defined(INT_SHAKELie) && !defined(INT_varint4lie)
      procedure   :: GL(INTEGRATOR)_g              => g
#endif
      procedure   :: GL(INTEGRATOR)_qlpexphDqtilde => qlpexphDqtilde
      procedure   :: GL(INTEGRATOR)_itlbtvtw       => itlbtvtw
      procedure   :: GL(INTEGRATOR)_tilde          => tilde ! dummy
      procedure   :: GL(INTEGRATOR)_Ct             => Ct
      procedure   :: GL(INTEGRATOR)_Kt             => Kt
      procedure   :: GL(INTEGRATOR)_Kt_lambda      => Kt_lambda
      procedure   :: GL(INTEGRATOR)_Tg             => Tg
      procedure   :: GL(INTEGRATOR)_norm           => norm ! dummy
      procedure   :: GL(INTEGRATOR)_outputFunction => outputFunction
      procedure   :: GL(INTEGRATOR)_init           => init
      procedure   :: GL(INTEGRATOR)_Phi            => Phi
      !! DEBUG
      procedure   :: GL(INTEGRATOR)_B              => B
      !procedure   :: GL(INTEGRATOR)_B              => GL(INTEGRATOR)_num_B
      !! GUBED
      procedure   :: GL(INTEGRATOR)_Z              => Z
      procedure   :: GL(INTEGRATOR)_matZ           => matZ ! dummy
      !
#if defined(INT_RATTLie) || defined(INT_SHAKELie) || defined(INT_varint4lie)
      procedure   :: GL(INTEGRATOR)_inertial       => inertial
      procedure   :: GL(INTEGRATOR)_f              => f
      procedure   :: GL(INTEGRATOR)_Tg_inv_T       => Tg_inv_T
      procedure   :: GL(INTEGRATOR)_d_Tg_inv_T     => d_Tg_inv_T
#endif
   end type heavy_top_t

   ! Functions and subroutines that the module contains:
   ! The procedures we need in order to solve the problem are
   ! actually implemented here.
   contains

      ! In linear spaces, tilde = id
      pure function tilde(this, v) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)  :: this
         real(8), dimension(:),  intent(in)  :: v
         ! result
         real(8), dimension(this%sizeq)      :: rslt
         !
         ! dummy function
         ERROR STOP "tilde is a dummy function and should not be reached"
      end function tilde

      ! q*exp(tilde(h Dq))
      pure function qlpexphDqtilde(this, q, h, Dq) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)  :: this
         real(8), dimension(:),  intent(in)  :: q
         real(8),                intent(in)  :: h
         real(8), dimension(:),  intent(in)  :: Dq
         ! result
         real(8), dimension(size(q))         :: rslt
         !
         select case (this%liegroup)
            case (1) ! SO(3)
               rslt = matmatv(q, so3exptv(h*Dq))
            case (2) ! S^3
               rslt = qp(q, expt_s3(h*Dq))
            case (3) ! direct product with SO(3)
               associate (Rv =>  q( 1: 9),&
                           x =>  q(10:12),&
                          DR => Dq( 1: 3),&
                          Dx => Dq( 4: 6))
                  rslt( 1: 9) = matmatv(Rv, so3exptv(h*DR))
                  rslt(10:12) = x + h*Dx
               end associate
            case (4) ! direct product with S^3
               associate ( p => q(1:4), &
                           x => q(5:7), &
                          Dp => Dq(1:3),&
                          Dx => Dq(4:6))
                  rslt(1:4) = qp(p, expt_s3(h*Dp))
                  rslt(5:7) = x + h*Dx
               end associate
            case (5) ! semidirect product with SO(3) (SE(3))
               associate (Dp => Dq(1:3),&
                          Dx => Dq(4:6))
                  rslt(1:9) = so3exptv(h*Dp)
                  !DEBUG print *, matmatv([rslt(1),rslt(4),rslt(7),rslt(2),rslt(5),rslt(8),rslt(3),rslt(6),rslt(9)],rslt(1:9))
                  rslt(10:12) = se3A(h*Dp, h*Dx)
                  rslt = se3pv(q, rslt)
               end associate

            case (6) ! semidirect product with S^3
               rslt = expt_s3sdr3(h*Dq)
               rslt = lp_s3sdr3(q, rslt)
            case (7) ! unit dual quaternions
               associate (Dp => Dq(1:3),&
                          Dx => Dq(4:6))
                  rslt(1:4) = expt_s3(h*Dp)
                  ! intermediate result
                  rslt(5) = 0.0_8
                  rslt(6:8) = se3A(h*Dp, h*Dx)/2
                  ! final
                  rslt(5:8) = qp(rslt(5:8),rslt(1:4))
                  rslt = dqp(q, rslt)
               end associate
         end select
      end function qlpexphDqtilde

      ! Essentially the Lie bracket of the vector space in which the velocities are measured
      pure function itlbtvtw(this, v, w) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)  :: this
         real(8), dimension(:),  intent(in)  :: v
         real(8), dimension(:),  intent(in)  :: w
         ! result
         real(8), dimension(size(v))         :: rslt
         select case (this%liegroup)
         case (1, 2)
               rslt = cross(v,w)
            case (3)
               associate (Om1 => v(1:3),&
                          Om2 => w(1:3))
                  rslt(1:3) = cross(Om1, Om2)
                  rslt(4:6) = 0.0_8
               end associate
            case (4)
               rslt(1:3) = cross(v(1:3), w(1:3))
               rslt(4:6) = 0.0_8
            case (5, 6, 7)
               associate (Om1 => v(1:3),&
                          Om2 => w(1:3),&
                           U1 => v(4:6),&
                           U2 => w(4:6))
                  rslt(1:3) = cross(Om1, Om2)
                  rslt(4:6) = cross(Om1, U2) - cross(Om2, U1)
               end associate
         end select
      end function itlbtvtw

      ! Tangent operator
      pure function Tg(this, h, dq) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)     :: this
         real(8),                intent(in)     :: h
         real(8), dimension(:),  intent(in)     :: dq
         ! result
         real(8), dimension(size(dq),size(dq))  :: rslt
         !
         select case (this%liegroup)
            case (1, 2)
               rslt = tan_op_s3(h*dq)
            case (3, 4)
               rslt(1:3,1:3) = tan_op_s3(h*dq(1:3))
               rslt(1:3,4:6) = 0.0_8
               rslt(4:6,1:3) = 0.0_8
               rslt(4:6,4) = [1.0_8,0.0_8,0.0_8]
               rslt(4:6,5) = [0.0_8,1.0_8,0.0_8]
               rslt(4:6,6) = [0.0_8,0.0_8,1.0_8]
            case (5, 6, 7)
               rslt = tan_op_s3sdr3(h*dq)
         end select
      end function Tg

      ! Transpose of inverse tangent operator
      pure function Tg_inv_T(this, dq) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)     :: this
         real(8), dimension(:),  intent(in)     :: dq
         ! result
         real(8), dimension(size(dq),size(dq))  :: rslt
         !
         select case (this%liegroup)
            case (1, 2)
               rslt = tan_tr_inv_s3(dq)
            case (3, 4)
               rslt(1:3,1:3) = tan_tr_inv_s3(dq)
               rslt(1:3,4:6) = 0.0_8
               rslt(4:6,1:3) = 0.0_8
               rslt(4:6,4) = [1.0_8,0.0_8,0.0_8]
               rslt(4:6,5) = [0.0_8,1.0_8,0.0_8]
               rslt(4:6,6) = [0.0_8,0.0_8,1.0_8]
            case (5, 6, 7)
               rslt = tan_tr_inv_s3sdr3(dq)
         end select
      end function Tg_inv_T

      ! derivative of (Transpose of inverse tangent operator multiplied by a vector)
      pure function d_Tg_inv_T(this, v, w) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)     :: this
         real(8), dimension(:),  intent(in)     :: v, w
         ! result
         real(8), dimension(this%sizev,this%sizev)  :: rslt
         !
         select case (this%liegroup)
            case (1, 2)
               rslt = d_tan_tr_inv_s3(v, w)
            case (3, 4)
               rslt = 0.0_8
               rslt(1:3,1:3) = d_tan_tr_inv_s3(v(1:3), w(1:3))
            case (5, 6, 7)
               rslt = d_tan_tr_inv_s3sdr3(v,w)
         end select
      end function d_Tg_inv_T

      pure function M(this, q) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)              :: this
         real(8), dimension(:),  intent(in)              :: q
         ! result
         real(8), dimension(this%sizev,this%sizev)       :: rslt
         !
         select case (this%liegroup)
            case (1, 2)
               rslt(1:3,1) = [this%inerJ(1),0.0_8,0.0_8]
               rslt(1:3,2) = [0.0_8,this%inerJ(2),0.0_8]
               rslt(1:3,3) = [0.0_8,0.0_8,this%inerJ(3)]
               rslt = rslt - this%mass * matmul(skw(this%refX),skw(this%refX))
            case (3:7)
               rslt(1:3,1) = [this%inerJ(1),0.0_8,0.0_8]
               rslt(1:3,2) = [0.0_8,this%inerJ(2),0.0_8]
               rslt(1:3,3) = [0.0_8,0.0_8,this%inerJ(3)]
               rslt(1:3,4:6) = 0.0_8
               rslt(4:6,1:3) = 0.0_8
               rslt(4:6,4) = [this%mass,0.0_8,0.0_8]
               rslt(4:6,5) = [0.0_8,this%mass,0.0_8]
               rslt(4:6,6) = [0.0_8,0.0_8,this%mass]
         end select
      end function M

      pure function diagM(this, q) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)  :: this
         real(8), dimension(:),  intent(in)  :: q
         ! result
         real(8), dimension(this%sizev)      :: rslt
         !
         ERROR STOP "function diagM should not be reached"
      end function diagM

      ! inertial forces (repeats code from mass matrix M)
      pure function inertial(this, v) result(rslt)
         ! input
         class(heavy_top_t), intent(in)  :: this
         real(8)         , intent(in)  :: v(:)
         ! result
         real(8)                       :: rslt(size(v))
         ! internal
         real(8)                       :: massm(3,3)
         !
         select case (this%liegroup)
            case (1, 2)
               massm(1:3,1) = [this%inerJ(1),0.0_8,0.0_8]
               massm(1:3,2) = [0.0_8,this%inerJ(2),0.0_8]
               massm(1:3,3) = [0.0_8,0.0_8,this%inerJ(3)]
               massm = massm - this%mass * matmul(skw(this%refX),skw(this%refX))
               rslt = cross(v, matmul(massm, v))
            case (3, 4)
               associate (Om => v(1:3))
                  rslt(1:3) = cross(Om, this%inerJ*Om)
                  rslt(4:6) = 0.0_8
               end associate
            case (5, 6, 7)
               associate ( Om => v( 1: 3),&
                            U => v( 4: 6))
                  rslt(1:3) = cross(Om, this%inerJ*Om)
                  rslt(4:6) = this%mass*cross(Om, U)
               end associate
         end select
         rslt = -rslt
      end function inertial

      pure function f(this, q, v, t) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)  :: this
         real(8), dimension(:),  intent(in)  :: q
         real(8), dimension(:),  intent(in)  :: v
         real(8),                intent(in)  :: t
         ! result
         real(8), dimension(size(v))         :: rslt
         !
         select case (this%liegroup)
            case (1)
               rslt =  this%mass*cross(this%refX, matTvecv(q, this%gravity))
            case (2)
               rslt =  this%mass*cross(this%refX, apply_conj_quat(q,this%gravity))
            case (3, 4)
                  rslt(1:3) = 0.0_8
                  rslt(4:6) = this%mass * this%gravity
            case (5)
               associate ( Rv => q( 1: 9),&
                            x => q(10:12),&
                           Om => v( 1: 3),&
                            U => v( 4: 6))
                  rslt(1:3) = 0.0_8
                  rslt(4:6) = matTvecv(Rv, this%mass*this%gravity)
               end associate
            case (6)
               associate ( p => q(1:4),&
                           x => q(5:7),&
                          Om => v(1:3),&
                           U => v(4:6))
                  rslt(1:3) = 0.0_8
                  rslt(4:6) = this%mass*apply_conj_quat(p,this%gravity)
               end associate
            case (7)
               associate (pR => q(1:4),&
                          pD => q(5:8),&
                          Om => v(1:3),&
                           U => v(4:6))
                  rslt(1:3) = 0.0_8
                  rslt(4:6) = this%mass*apply_conj_quat(pR,this%gravity)
               end associate
         end select
         ! for legacy reasons, we have to give back the negative
         rslt = -rslt
      end function f

      pure function g(this, q, v, t) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)  :: this
         real(8), dimension(:),  intent(in)  :: q
         real(8), dimension(:),  intent(in)  :: v
         real(8),                intent(in)  :: t
         ! result
         real(8), dimension(size(v))         :: rslt
         !
         rslt = -inertial(this, v) + f(this, q, v, t)
      end function g


      pure function Phi(this,q) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)  :: this
         real(8), dimension(:),  intent(in)  :: q
         ! result
         real(8), dimension(this%sizel)      :: rslt
         !
         select case (this%liegroup)
            case (1, 2)
               ERROR STOP "Function Phi should not be reached for unconstrained formulations"
            case (4, 6)
               associate (p => q(1:4),&
                          x => q(5:7))
                  rslt = apply_conj_quat(p, x) - this%refX
               end associate
            case (3, 5)
               associate (Rv => q( 1: 9),&
                           x => q(10:12))
                 rslt = matTvecv(Rv,x) - this%refX
               end associate
            case (7)
               associate (pR => q(1:4),&
                          pD => q(5:8))
                  rslt = 2*qikpq(pR, pD) - this%refX
               end associate
         end select
      end function Phi

! DEBUG
   pure function num_B(this, q) result(rslt)
      ! input
      class(heavy_top_t),     intent(in)        :: this
      real(8), dimension(:),  intent(in)        :: q
      ! result
      real(8), dimension(this%sizel,this%sizev) :: rslt
      ! internal
      integer                                   :: i
      real(8), dimension(this%sizev)            :: w
      real(8), dimension(this%sizev)            :: Phi0
      !
      !
      ! Calculate Phi
      Phi0 = this%GL(INTEGRATOR)_Phi(q)

      ! set w to zero
      w = 0.0_8
      ! loop over the columns of Ct
      do i=1,this%sizev
         if (.not. i == 1) w(i-1) = 0.0_8
         w(i) = max( abs(q(i))*1.0e-8_8,  1.0e-12_8)*1.0e8
         rslt(:,i) = (this%GL(INTEGRATOR)_Phi(this%GL(INTEGRATOR)_qlpexphDqtilde(q,1.0_8, w)) - Phi0) / w(i)
         !w    = 0.0_8
         !w(i) = h
         !rslt(:,i) = (this%GL(INTEGRATOR)_g(q+this%GL(INTEGRATOR)_tilde(w),v,t) - g0)/h
      end do
   end function num_B
! GUBED

      pure function B(this,q) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)           :: this
         real(8), dimension(:),  intent(in)           :: q
         ! result
         real(8), dimension(this%sizel,this%sizev)    :: rslt
         !
         select case (this%liegroup)
            case (1, 2)
               ERROR STOP "function B should not be reached in unconstrained formulations"
            case (3)
               associate (Rv => q(1:9))
                  rslt(:,1:3) = skw(this%refX)
                  rslt(1,4:6) = Rv(1:3)
                  rslt(2,4:6) = Rv(4:6)
                  rslt(3,4:6) = Rv(7:9)
               end associate
            case (4)
               associate (p => q(1:4))
                  rslt(:,1:3) = skw(this%refX)
                  rslt(:,4:6) = transpose(qR(p))
               end associate
            case (5, 6, 7)
               rslt(:,1:3) = skw(this%refX)
               rslt(:,4) = [1.0_8, 0.0_8, 0.0_8]
               rslt(:,5) = [0.0_8, 1.0_8, 0.0_8]
               rslt(:,6) = [0.0_8, 0.0_8, 1.0_8]
         end select
      end function B


      pure function Z(this,q,v) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)  :: this
         real(8), dimension(:),  intent(in)  :: q
         real(8), dimension(:),  intent(in)  :: v
         ! result
         real(8), dimension(this%sizel)      :: rslt
         !
         select case (this%liegroup)
            case (1, 2)
               ERROR STOP "function Z should not be reached in unconstrained formulations"
            case (3)
               associate (Rv => q( 1: 9),&
                           x => q(10:12),&
                          Om => v( 1: 3),&
                           u => v( 4: 6))
                  rslt = cross(matTvecv(Rv, u), Om)
               end associate
            case (4)
               associate ( p => q(1:4),&
                           x => q(5:7),&
                          Om => v(1:3),&
                           u => v(4:6))
                  rslt = cross(apply_conj_quat(p, u), Om)
               end associate
            case (5, 6, 7)
               rslt = 0.0_8
         end select
      end function Z

      pure function matZ(this,q,v,T) result(rslt)
         ! input
         class(heavy_top_t),        intent(in)        :: this
         real(8), dimension(:),     intent(in)        :: q
         real(8), dimension(:),     intent(in)        :: v
         real(8), dimension(:,:),   intent(in)        :: T
         ! result
         real(8), dimension(this%sizel, this%sizev)   :: rslt
         !
         ERROR STOP "function matZ should not be reached"
      end function matZ

      pure function norm(this, v) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)  :: this
         real(8), dimension(:),  intent(in)  :: v
         ! result
         real(8)                             :: rslt
         !
         ERROR STOP "function norm should not be reached"
      end function norm

      ! The derivative  of g wrt. v
      pure function Ct(this, q, v, t) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)  :: this
         real(8), dimension(:),  intent(in)  :: q
         real(8), dimension(:),  intent(in)  :: v
         real(8),                intent(in)  :: t
         ! result
         real(8), dimension(size(v),size(v)) :: rslt
         ! internal
         integer                             :: i
         !
         select case (this%liegroup)
            case (1, 2)
               rslt = matmul(skw(v), M(this,q)) - skw(matmul(M(this,q), v))
            case (3, 4)
               associate (Om => v(1:3))
                  rslt(1:3,1:3) = skw(Om)
                  forall (i=1:3)
                     rslt(1:3,i) = rslt(1:3,i)*this%inerJ(i)
                  end forall
                  rslt(1:3,1:3) = rslt(1:3,1:3) - skw(this%inerJ * Om)
                  rslt(4:6,1:3) = 0.0_8
                  rslt(:,4:6) = 0.0_8
               end associate
            case (5, 6, 7)
               associate (Om => v(1:3),&
                           U => v(4:6))
                  rslt(1:3,1:3) = skw(Om)
                  forall (i=1:3)
                     rslt(1:3,i) = rslt(1:3,i)*this%inerJ(i)
                  end forall
                  rslt(1:3,1:3) = rslt(1:3,1:3) - skw(this%inerJ * Om)
                  !
                  rslt(4:6,1:3) = - skw(this%mass*U)
                  rslt(1:3,4:6) = 0.0_8
                  rslt(4:6,4:6) = skw(this%mass*Om)
               end associate
         end select
      end function Ct

      ! The derivative of g wrt q
      pure function Kt(this, q, v, vd, t) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)  :: this
         real(8), dimension(:),  intent(in)  :: q
         real(8), dimension(:),  intent(in)  :: v
         real(8), dimension(:),  intent(in)  :: vd
         real(8),                intent(in)  :: t
         ! result
         real(8), dimension(size(v),size(v)) :: rslt
         !
         select case (this%liegroup)
            case (1)
               rslt = - this%mass * &
                      matmul(skw(this%refX), skw(matTvecv(q, this%gravity)))
            case (2)
               rslt = - this%mass * &
                      matmul(skw(this%refX), skw(apply_conj_quat(q,this%gravity)))
            case (3:7)
               ERROR STOP "function Kt should not be reached for constrained formulations"
         end select
      end function Kt

      ! The derivative of g wrt q
      pure function Kt_lambda(this, q, v, vd, l, t) result(rslt)
         ! input
         class(heavy_top_t),     intent(in)  :: this
         real(8), dimension(:),  intent(in)  :: q
         real(8), dimension(:),  intent(in)  :: v
         real(8), dimension(:),  intent(in)  :: vd
         real(8), dimension(:),  intent(in)  :: l
         real(8),                intent(in)  :: t
         ! result
         real(8), dimension(size(v),size(v)) :: rslt
         ! internal
         integer                             :: m
         !
         select case (this%liegroup)
            case (1:2)
               ERROR STOP "function Kt_lambda should not be reached for unconstrained formulations"
            case (3)
               associate (Rv => q(1:9))
                  rslt(1:3,1:3) = 0.0_8
                  rslt(4:6,1:3) = matvmat2mat(Rv, skw(l))
                  rslt( : ,4:6) = 0.0_8
               end associate
            case (4)
               associate (p => q(1:4))
                  rslt(1:3,1:3) = 0.0_8
                  rslt(4:6,1:3) = matmul(qR(p), skw(l))
                  rslt( : ,4:6) = 0.0_8
               end associate
            case (5)
               associate (Rv => q(1:9))
                  rslt(1:3,1:3) = 0.0_8
                  rslt(4:6,1:3) = - skw(this%mass * matTvecv(Rv,this%gravity))
                  rslt( : ,4:6) = 0.0_8
               end associate
            case (6)
               associate (p => q(1:4))
                  rslt(1:3,1:3) = 0.0_8
                  rslt(4:6,1:3) = - skw(this%mass * apply_conj_quat(p, this%gravity))
                  rslt( : ,4:6) = 0.0_8
               end associate
            case (7)
               associate (pR => q(1:4))
                  rslt(1:3,1:3) = 0.0_8
                  rslt(4:6,1:3) = - skw(this%mass * apply_conj_quat(pR, this%gravity))
                  rslt( : ,4:6) = 0.0_8
               end associate
         end select
      end function Kt_lambda

      pure function modmod(r1,r2) result(rslt)
         ! input
         real(8), intent(in)  :: r1
         real(8), intent(in)  :: r2
         ! result
         real(8)              :: rslt
         !
         rslt = r1 - nint(r1/r2)*r2
      end function modmod

      subroutine outputFunction(this,info)
         ! input
         class(heavy_top_t),     intent(in)  :: this
         integer,                intent(in)  :: info

         select case (info)
            case (0)
               ! initialization:

               ! Binary output file and misc output file were opened in
               ! main, because 'this' is intent(in)

            case (1)
               if (this%output_t_at == 0 .or. &
                   abs(modmod(this%t,this%t_output_at_multiples_of))  < 1.0e-9_8) then ! TODO: Besserer Wert? DEBUG

                  if (this%output_type == 1) then
                     ! normal output:
                     ! write time to file
                     write (this%out_bin_lun) this%t
                     !
                     if (this%opts%constrained == 1) then
                        if (this%opts%stab2 == 1) then
                           write (this%out_bin_lun) this%q, this%v, &
#if defined(INT_RATTLie) || defined(INT_varint4lie)
                              this%l, this%lm, &
#elif defined(INT_SHAKELie)
                              this%l, this%eta, &
#else
                              this%vd, this%l, this%eta, &
#endif
                              this%GL(INTEGRATOR)_phi(this%q), matmul(this%GL(INTEGRATOR)_B(this%q), this%v)
                        else
#if defined(INT_RATTLie)
                              ERROR STOP "RATTLie does not support index-3"
#elif defined(INT_varint4lie)
                              ERROR STOP "varint4lie does not support index-3"
#elif defined(INT_SHAKELie)
                           write (this%out_bin_lun) this%q, this%v, &
                              this%l, this%GL(INTEGRATOR)_phi(this%q), matmul(this%GL(INTEGRATOR)_B(this%q), this%v)
#else
                           write (this%out_bin_lun) this%q, this%v, &
                              this%vd, this%l, this%GL(INTEGRATOR)_phi(this%q), matmul(this%GL(INTEGRATOR)_B(this%q), this%v)
#endif
                        end if
                     else
#if defined(INT_RATTLie) || defined(INT_SHAKELie) || defined(INT_varint4lie)
                        write (this%out_bin_lun) this%q, this%v
#else
                        write (this%out_bin_lun) this%q, this%v, this%vd
#endif
                     end if
                     !! DEBUG
                     !flush (this%out_bin_lun)
                     !! GUBED
                  else if (this%output_type == 2) then
                     ERROR STOP "output_type = 2 not implemented"
                  end if

                  !! write misc output (DEBUG)
                  !write (this%out_misc_lun,*), this%t, this%v, !this%GL(INTEGRATOR)_stats%newt_steps_curr, this%GL(INTEGRATOR)_stats%newt_steps_max, !this%GL(INTEGRATOR)_stats%newt_steps_avg
                  !flush (this%out_misc_lun)
               end if
            case (99)
               ! termination:

               ! Closing binary output file is done in main

            case default
               print *, 'FATAL Error: Got unsupported info flag ', info, ' in outputFunction'
               ERROR STOP
         end select
      end subroutine outputFunction

      subroutine init(this)
         ! input/output
         class(heavy_top_t),     intent(inout)  :: this

         ! Set lengths
         select case (this%liegroup)
            case (1) ! SO(3)
               this%sizeq = 9
               this%sizev = 3
            case (2) ! S^3
               this%sizeq = 4
               this%sizev = 3
            case (3) ! SO(3) + R^3
               this%sizeq = 9+3
               this%sizev = 3+3
               this%sizel = 3
            case (4) ! S^3 x R^3
               this%sizeq = 4+3
               this%sizev = 3+3
               this%sizel = 3
            case (5) ! SE(3)
               this%sizeq = 9+3
               this%sizev = 3+3
               this%sizel = 3
            case (6) ! S^3 |x R^3
               this%sizeq = 4+3
               this%sizev = 3+3
               this%sizel = 3
            case (7) ! UDQ
               this%sizeq = 4+4
               this%sizev = 3+3
               this%sizel = 3
         end select

         ! Allocate position, velocity, acceleration, acceleration-like-varible
         allocate(this%q(this%sizeq))
         allocate(this%v(this%sizev))
#if defined(INT_RATTLie) || defined(INT_varint4lie)
         allocate(this%p(this%sizev))
#endif
#if !defined(INT_RATTLie) && !defined(INT_SHAKELie) && !defined(INT_varint4lie)
         allocate(this%vd(this%sizev))
#endif
#ifdef INT_gena
         allocate(this%a(this%sizev))
#endif
         ! Allocate Lagrange multipliers and auxiliar variable eta if applicable
         if (this%opts%constrained == 1) then
            ! allocate this%l
            if (allocated(this%l)) then
               deallocate(this%l)
            end if
            allocate(this%l(this%sizel))

#if defined(INT_RATTLie) || defined(INT_varint4lie)
            ! allocate this%l
            if (allocated(this%lm)) then
               deallocate(this%lm)
            end if
            allocate(this%lm(this%sizel))
            ! allocate this%l
            if (allocated(this%lp)) then
               deallocate(this%lp)
            end if
            allocate(this%lp(this%sizel))
#endif
#if !defined(INT_RATTLie) && !defined(INT_varint4lie)
            ! In the stabilized index-2 case
            if ( this%opts%stab2 == 1 ) then
               ! allocate eta
               if (allocated(this%eta)) then
                  deallocate(this%eta)
               end if
               allocate(this%eta(this%sizel))
               ! eta vanishes analytically, so \eta(t_0) = 0
               this%eta = 0.0_8
            end if
#endif
         end if

         ! Calculate the top's reference position
         this%refX = apply_conj_quat(this%p0, this%x0)

         ! Set initial values
         this%t = this%opts%t0
         select case (this%liegroup)
            case (1) ! SO(3)
               this%q = qRv(this%p0)
               this%v = this%Om0
            case (2) ! S^3
               this%q = this%p0
               this%v = this%Om0
            case (3) ! SO(3) + R^3
               this%q( 1: 9) = qRv(this%p0)
               this%q(10:12) = this%x0
               this%v( 1: 3) = this%Om0
               this%v( 4: 6) = matTvecv(this%q(1:9), cross(this%Om0, this%refX))
            case (4) ! S^3 x R^3
               this%q(1:4) = this%p0
               this%q(5:7) = this%x0
               this%v(1:3) = this%Om0
               this%v(4:6) = apply_quat(this%p0, cross(this%Om0, this%refX))
            case (5) ! SE(3)
               this%q( 1: 9) = qRv(this%p0)
               this%q(10:12) = this%x0
               this%v( 1: 3) = this%Om0
               this%v(4:6) = cross(this%Om0, this%refX)
            case (6) ! S^3 |x R^3
               this%q(1:4) = this%p0
               this%q(5:7) = this%x0
               this%v(1:3) = this%Om0
               this%v(4:6) = cross(this%Om0, this%refX)
            case (7) ! UDQ
               this%q(1:4) = this%p0
               ! intermediate result
               this%q(5) = 0.0_8
               this%q(6:8) = this%x0/2
               ! final
               this%q(5:8) = qp(this%q(5:8), this%q(1:4))
               this%v(1:3) = this%Om0
               this%v(4:6) = cross(this%Om0, this%refX)
         end select
      end subroutine init

end module heavy_top
