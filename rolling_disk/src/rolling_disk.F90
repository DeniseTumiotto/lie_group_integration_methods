#define GL(x) x

module testprobrolling_disk
   use GL(INTEGRATOR)
   use s3sdr3_functions
   implicit none

   ! extension of the type GL(INTEGRATOR)_problem
   type, extends(GL(INTEGRATOR)_problem)   :: myproblem
      ! additional parameters
      !
      ! name of the problem
      character(len=256)      :: probname
      !
      character(len=256)      :: out_fname
      !
      integer                 :: out_bin_lun
      integer                 :: out_misc_lun
      !
      real(8)                 :: mass
      real(8), dimension(3)   :: inerJ ! diagonal elements of inertial tensor (which is a diag matrix)
      real(8), dimension(3)   :: gravity ! direction and magnitude of gravity
      !
      real(8), dimension(3)   :: x0 ! initial position
      real(8), dimension(4)   :: p0 ! initial orientation
      real(8), dimension(3)   :: Om0 ! initial angular velocity (initial velocity is calculated from this)
      !
      integer                 :: output_t_at = 0
      real(8)                 :: t_output_at_multiples_of = 0.0_8
   contains

      ! referencing the former deferred procedures
      procedure   :: GL(INTEGRATOR)_M              => myM
      procedure   :: GL(INTEGRATOR)_diag_M         => mydiagM
      procedure   :: GL(INTEGRATOR)_f              => rolling_disk_f
      procedure   :: GL(INTEGRATOR)_inertial       => rolling_disk_inertial
      procedure   :: GL(INTEGRATOR)_qlpexphDqtilde => myqlpexphDqtilde
      procedure   :: GL(INTEGRATOR)_itlbtvtw       => myitlbtvtw
      procedure   :: GL(INTEGRATOR)_tilde          => myembeddingop
      procedure   :: GL(INTEGRATOR)_Ct             => myCt
      procedure   :: GL(INTEGRATOR)_Kt             => myKt
      procedure   :: GL(INTEGRATOR)_Kt_lambda      => myKtl
      procedure   :: GL(INTEGRATOR)_Tg_inv_T       => my_Tg_inv_T
      procedure   :: GL(INTEGRATOR)_d_Tg_inv_T     => my_d_Tg_inv_T
      procedure   :: GL(INTEGRATOR)_Tg             => myTg
      procedure   :: GL(INTEGRATOR)_norm           => mynorm
      procedure   :: GL(INTEGRATOR)_outputFunction => myoutputFunction
      procedure   :: GL(INTEGRATOR)_init           => myinit
      procedure   :: GL(INTEGRATOR)_phi            => myphi
      procedure   :: GL(INTEGRATOR)_b              => myb
      procedure   :: GL(INTEGRATOR)_Bnh            => myBnh
      procedure   :: GL(INTEGRATOR)_Z              => myZ
      procedure   :: GL(INTEGRATOR)_Znh            => myZnh
      procedure   :: GL(INTEGRATOR)_matZ           => mymatZ
   end type myproblem

   ! functions and subroutines that the module contains
   ! here the procedures we need in order to solve the problem are
   ! actually implemented
   contains

      pure function myembeddingop(this, v) result(rslt)
         ! input
         class(myproblem),      intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         ! result
         real(8), dimension(this%sizeq)      :: rslt
         !
         ERROR STOP "tilde is a dummy function and should not be reached"
      end function myembeddingop

      pure function myqlpexphDqtilde(this, q, h, Dq) result(rslt)
         ! input
         class(myproblem),      intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8),               intent(in)   :: h
         real(8), dimension(:), intent(in)   :: Dq
         ! result
         real(8), dimension(size(q))         :: rslt
         !
         rslt = lp_s3sdr3(q, expt_s3sdr3(h*Dq))
      end function myqlpexphDqtilde

      pure function myitlbtvtw(this, v, w) result(rslt)
         ! input
         class(myproblem),      intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         real(8), dimension(:), intent(in)   :: w
         ! result
         real(8), dimension(size(v))         :: rslt
         !
         rslt = lie_bracket_s3sdr3(v, w)
      end function myitlbtvtw

      pure function myTg(this, h, dq) result(rslt)
         ! input
         class(myproblem),      intent(in)      :: this
         real(8),               intent(in)      :: h
         real(8), dimension(:), intent(in)      :: dq
         ! result
         real(8), dimension(size(dq),size(dq))  :: rslt
         !
         rslt = tan_op_s3sdr3(h*dq)
      end function myTg

      pure function my_Tg_inv_T(this, dq) result(rslt)
         ! input
         class(myproblem),      intent(in)         :: this
         real(8), dimension(:), intent(in)         :: dq
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
         !
         rslt = tan_tr_inv_s3sdr3(dq)
      end function my_Tg_inv_T

      pure function my_d_Tg_inv_T(this, v, w) result(rslt)
         ! input
         class(myproblem),      intent(in)         :: this
         real(8), dimension(:), intent(in)         :: v, w
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
         !
         rslt = d_tan_tr_inv_s3sdr3(v, w)
      end function my_d_Tg_inv_T

      pure function myM(this, q) result(rslt)
         ! input
         class(myproblem),      intent(in)               :: this
         real(8), dimension(:), intent(in)               :: q
         ! result
         real(8), dimension(this%sizev,this%sizev)       :: rslt
         !
         rslt(1:3,1) = [this%inerJ(1),0.0_8,0.0_8]
         rslt(1:3,2) = [0.0_8,this%inerJ(2),0.0_8]
         rslt(1:3,3) = [0.0_8,0.0_8,this%inerJ(3)]
         rslt(1:3,4:6) = 0.0_8
         rslt(4:6,1:3) = 0.0_8
         rslt(4:6,4) = [this%mass,0.0_8,0.0_8]
         rslt(4:6,5) = [0.0_8,this%mass,0.0_8]
         rslt(4:6,6) = [0.0_8,0.0_8,this%mass]
      end function myM

      pure function mydiagM(this, q) result(rslt)
         ! input
         class(myproblem),      intent(in)               :: this
         real(8), dimension(:), intent(in)               :: q
         ! result
         real(8), dimension(this%sizev)                  :: rslt
         !
         rslt(1:3) = this%inerJ
         rslt(4:6) = this%mass
      end function mydiagM

      pure function rolling_disk_f(this, q, v, t) result(rslt)
         ! input
         class(myproblem),      intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8), dimension(:), intent(in)   :: v
         real(8),               intent(in)   :: t
         ! result
         real(8), dimension(size(v))         :: rslt
         !
         rslt(1:3) = 0.0_8
         rslt(4:6) = this%mass * apply_conj_quat(q(1:4),this%gravity)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! Don't forget the minus
         rslt = -rslt

      end function rolling_disk_f

      pure function rolling_disk_inertial(this, v) result(rslt)
         ! input
         class(myproblem), intent(in)  :: this
         real(8)         , intent(in)  :: v(:)
         ! result
         real(8)                       :: rslt(size(v))
         !
         associate (Om => v(1:3), U => v(4:6))
            rslt(1:3) = cross(Om, this%inerJ*Om)
            rslt(4:6) = this%mass*cross(Om,U)
         end associate
      end function rolling_disk_inertial

      pure function myphi(this,q) result(rslt)
         ! input
         class(myproblem),                 intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         ! result
         real(8), dimension(this%sizel)               :: rslt
         ! internal
         integer                                      :: i
         !
#ifdef ONLY_NONHOL
         real(8), dimension(0)   :: my
         !
         rslt = my
#else
         associate (p0 => q(1), p1 => q(2), p2 => q(3), p3 => q(4), &
                    y1 => q(5), y2 => q(6), y3 => q(7))
            rslt = -2*(p1**2 + p2**2)*(p0**2 + p3**2) + y3**2/2
            !rslt = y3 - 1.0_8
         end associate
#endif
      end function myphi

      pure function myB(this,q) result(rslt)
         ! input
         class(myproblem), intent(in)  :: this
         real(8),          intent(in)  :: q(:)
         ! result
         real(8)                       :: rslt(this%sizel,this%sizev)
         !
#ifdef ONLY_NONHOL
         real(8)  :: my(0,6)
         !
         rslt = my
#else
         associate (p0 => q(1), p1 => q(2), p2 => q(3), p3 => q(4), &
                    y1 => q(5), y2 => q(6), y3 => q(7))
            rslt(1,:) = [2*(p0*p1 + p2*p3)*(-p0**2 + p1**2 + p2**2 - p3**2),-2*(p0*p2 - p1*p3)*(p0**2 - p1**2 - p2**2 + p3**2),0.0_8,(-2*p0*p2 + 2*p1*p3)*y3,2*(p0*p1 + p2*p3)*y3,(p0**2 - p1**2 - p2**2 + p3**2)*y3]
            !rslt(1,1:3) = 0.0_8
            !rslt(1,4:6) = [-2*p0*p2 + 2*p1*p3,2*(p0*p1 + p2*p3),p0**2 - p1**2 - p2**2 + p3**2]
         end associate
#endif
      end function myb

      pure function myZ(this,q,v) result(rslt)
         ! input
         class(myproblem), intent(in)  :: this
         real(8),          intent(in)  :: q(:)
         real(8),          intent(in)  :: v(:)
         ! result
         real(8)                       :: rslt(this%sizel)
         !
#ifdef ONLY_NONHOL
         real(8) :: my(0)
         !
         rslt = my
#else
         associate (p0 => q(1), p1 => q(2), p2 => q(3), p3 => q(4), &
                    y1 => q(5), y2 => q(6), y3 => q(7),             &
                    Om1 => v(1), Om2 => v(2), Om3 => v(3),          &
                    V1 => v(4), V2 => v(5), V3 => v(6)               )
            rslt(1) = 2*Om2*Om3*(p0*p1 + p2*p3)*(p0**2 - p1**2 - p2**2 + p3**2) - Om2**2*(p0**2 - p1**2 + 2*p0*p2 - p2**2 - 2*p1*p3 + p3**2)*(p0**2 - p1**2 - 2*p0*p2 - p2**2 + 2*p1*p3 + p3**2) - Om1**2*(p0**2 - 2*p0*p1 - p1**2 - p2**2 - 2*p2*p3 + p3**2)*(p0**2 + 2*p0*p1 - p1**2 - p2**2 + 2*p2*p3 + p3**2) - 2*Om1*(p0*p2 - p1*p3)*(-4*Om2*(p0*p1 + p2*p3) + Om3*(p0**2 - p1**2 - p2**2 + p3**2)) + 4*p0**2*p2**2*V1**2 - 8*p0*p1*p2*p3*V1**2 + 4*p1**2*p3**2*V1**2 - 8*p0**2*p1*p2*V1*V2 + 8*p0*p1**2*p3*V1*V2 - 8*p0*p2**2*p3*V1*V2 + 8*p1*p2*p3**2*V1*V2 + 4*p0**2*p1**2*V2**2 + 8*p0*p1*p2*p3*V2**2 + 4*p2**2*p3**2*V2**2 - 4*p0**3*p2*V1*V3 + 4*p0*p1**2*p2*V1*V3 + 4*p0*p2**3*V1*V3 + 4*p0**2*p1*p3*V1*V3 - 4*p1**3*p3*V1*V3 - 4*p1*p2**2*p3*V1*V3 - 4*p0*p2*p3**2*V1*V3 + 4*p1*p3**3*V1*V3 + 4*p0**3*p1*V2*V3 - 4*p0*p1**3*V2*V3 - 4*p0*p1*p2**2*V2*V3 + 4*p0**2*p2*p3*V2*V3 - 4*p1**2*p2*p3*V2*V3 - 4*p2**3*p3*V2*V3 + 4*p0*p1*p3**2*V2*V3 + 4*p2*p3**3*V2*V3 + p0**4*V3**2 - 2*p0**2*p1**2*V3**2 + p1**4*V3**2 - 2*p0**2*p2**2*V3**2 + 2*p1**2*p2**2*V3**2 + p2**4*V3**2 + 2*p0**2*p3**2*V3**2 - 2*p1**2*p3**2*V3**2 - 2*p2**2*p3**2*V3**2 + p3**4*V3**2 + 2*Om3*(p0*p1*V1 + p2*p3*V1 + p0*p2*V2 - p1*p3*V2)*y3 + Om2*(-(p0**2*V1) + (p1**2 + p2**2 - p3**2)*V1 - 2*p0*p2*V3 + 2*p1*p3*V3)*y3 + Om1*((p0**2 - p1**2 - p2**2 + p3**2)*V2 - 2*(p0*p1 + p2*p3)*V3)*y3
            !rslt(1) = -2*Om2*(p0*p2 + p1*p3)*V1 + Om3*(2*p1*p2*V1 - 2*p0*p3*V1 - p0**2*V2 - p1**2*V2 + (p2**2 + p3**2)*V2) + Om2*(p0**2 + p1**2 - p2**2 - p3**2)*V3 + 2*Om1*(p0*p2*V2 + p1*p3*V2 - p1*p2*V3 + p0*p3*V3)
         end associate
#endif
      end function myZ

      pure function myBnh(this,q) result(rslt)
         ! input
         class(myproblem), intent(in)  :: this
         real(8),          intent(in)  :: q(:)
         ! result
         real(8)                       :: rslt(this%sizelnh,this%sizev)
         !
         associate (p0 => q(1), p1 => q(2), p2 => q(3), p3 => q(4), &
                    y1 => q(5), y2 => q(6), y3 => q(7))
#ifdef ONLY_NONHOL
            rslt(1,:) = [0.0_8,0.0_8,2*(p0*p1 + p2*p3),y3,0.0_8,0.0_8]
            rslt(2,:) = [0.0_8,0.0_8,2*p0*p2 - 2*p1*p3,0.0_8,y3,0.0_8]
            rslt(3,:) = [-2*(p0*p1 + p2*p3),-2*p0*p2 + 2*p1*p3,0.0_8,0.0_8,0.0_8,y3]
#else
            rslt(1,:) = [-4*(p0*p2 + p1*p3)*(p0*p1 + p2*p3),-4*p0**2*p2**2 + 4*p1**2*p3**2,2*(p0*p1 - p2*p3)*(p0**2 + p1**2 + p2**2 + p3**2),(p0**2 + p1**2 - p2**2 - p3**2)*y3,2*(p1*p2 - p0*p3)*y3,2*(p0*p2 + p1*p3)*y3]
            rslt(2,:) = [4*(p0*p1 - p2*p3)*(p0*p1 + p2*p3),4*(p0*p2 - p1*p3)*(p0*p1 - p2*p3),2*(p0*p2 + p1*p3)*(p0**2 + p1**2 + p2**2 + p3**2),2*(p1*p2 + p0*p3)*y3,(p0**2 - p1**2 + p2**2 - p3**2)*y3,(-2*p0*p1 + 2*p2*p3)*y3]

            !rslt(1,:) = [1.0_8, 0.0_8, 0.0_8, y3, 0.0_8, 0.0_8]
            !rslt(2,:) = [0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, y3]
#endif

         end associate
      end function mybnh

      pure function myZnh(this,q,v) result(rslt)
         ! input
         class(myproblem), intent(in)  :: this
         real(8),          intent(in)  :: q(:)
         real(8),          intent(in)  :: v(:)
         ! result
         real(8)                       :: rslt(this%sizelnh)
         !
         associate (p0 => q(1), p1 => q(2), p2 => q(3), p3 => q(4), &
                    y1 => q(5), y2 => q(6), y3 => q(7),             &
                    Om1 => v(1), Om2 => v(2), Om3 => v(3),          &
                    V1 => v(4), V2 => v(5), V3 => v(6)               )
#ifdef ONLY_NONHOL
            rslt = [Om1*Om3*(p0**2 - p1**2 - p2**2 + p3**2) + 2*(p0*p2 - p1*p3)*(Om3 - V1)*(Om3 + V1) + 2*(p0*p1 + p2*p3)*V1*V2 + (p0**2 - p1**2 - p2**2 + p3**2)*V1*V3,-2*Om3**2*(p0*p1 + p2*p3) + Om2*Om3*(p0**2 - p1**2 - p2**2 + p3**2) + 2*V2*(-(p0*p2*V1) + p1*p3*V1 + p0*p1*V2 + p2*p3*V2) + (p0**2 - p1**2 - p2**2 + p3**2)*V2*V3,2*Om1*Om3*(-(p0*p2) + p1*p3) + Om1**2*(-p0**2 + p1**2 + p2**2 - p3**2) + Om2*(2*Om3*(p0*p1 + p2*p3) + Om2*(-p0**2 + p1**2 + p2**2 - p3**2)) + 2*(-(p0*p2*V1) + p1*p3*V1 + p0*p1*V2 + p2*p3*V2)*V3 + (p0**2 - p1**2 - p2**2 + p3**2)*V3**2]
#else
            rslt(1) = -2*Om2*Om3*(p1*p2 - p0*p3)*(-p0**2 + p1**2 + p2**2 - p3**2) + 4*Om2**2*(-(p0**3*p2) + p0*p2**3 + p1*(p1 - p3)*p3*(p1 + p3)) + 2*Om1**2*(-(p0**3*p2) - 3*p0**2*p1*p3 + p0*p2*(3*p1**2 + p2**2 - 3*p3**2) + p1*p3*(p1**2 + 3*p2**2 - p3**2)) - 2*p0**3*p2*V1**2 - 2*p0*p1**2*p2*V1**2 + 2*p0*p2**3*V1**2 + 2*p0**2*p1*p3*V1**2 + 2*p1**3*p3*V1**2 - 2*p1*p2**2*p3*V1**2 + 2*p0*p2*p3**2*V1**2 - 2*p1*p3**3*V1**2 + 2*p0**3*p1*V1*V2 + 2*p0*p1**3*V1*V2 - 6*p0*p1*p2**2*V1*V2 + 6*p0**2*p2*p3*V1*V2 + 6*p1**2*p2*p3*V1*V2 - 2*p2**3*p3*V1*V2 - 6*p0*p1*p3**2*V1*V2 - 2*p2*p3**3*V1*V2 + 4*p0*p1**2*p2*V2**2 - 4*p0**2*p1*p3*V2**2 + 4*p1*p2**2*p3*V2**2 - 4*p0*p2*p3**2*V2**2 + p0**4*V1*V3 - p1**4*V1*V3 - 6*p0**2*p2**2*V1*V3 + p2**4*V1*V3 + 6*p1**2*p3**2*V1*V3 - p3**4*V1*V3 + 6*p0**2*p1*p2*V2*V3 - 2*p1**3*p2*V2*V3 - 2*p1*p2**3*V2*V3 - 2*p0**3*p3*V2*V3 + 6*p0*p1**2*p3*V2*V3 + 6*p0*p2**2*p3*V2*V3 + 6*p1*p2*p3**2*V2*V3 - 2*p0*p3**3*V2*V3 + 2*p0**3*p2*V3**2 - 2*p0*p1**2*p2*V3**2 - 2*p0*p2**3*V3**2 + 2*p0**2*p1*p3*V3**2 - 2*p1**3*p3*V3**2 - 2*p1*p2**2*p3*V3**2 + 2*p0*p2*p3**2*V3**2 + 2*p1*p3**3*V3**2 + Om3*(2*p1*p2*V1 - 2*p0*p3*V1 - p0**2*V2 - p1**2*V2 + (p2**2 + p3**2)*V2)*y3 + Om2*(-2*(p0*p2 + p1*p3)*V1 + (p0**2 + p1**2 - p2**2 - p3**2)*V3)*y3 + Om1*(2*Om2*(-(p0*p1*(p0**2 + p1**2 - 3*p2**2)) + p2*(-3*(p0**2 + p1**2) + p2**2)*p3 + 3*p0*p1*p3**2 + p2*p3**3) + Om3*((p0**2 - p2**2)**2 - (p1**2 - p3**2)**2) + 2*(p0*p2*V2 + p1*p3*V2 - p1*p2*V3 + p0*p3*V3)*y3)
            rslt(2) = 2*Om1*Om3*(p1*p2 + p0*p3)*(p0**2 - p1**2 - p2**2 + p3**2) + 4*Om1**2*(p0**3*p1 - p0*p1**3 + p2*(p2 - p3)*p3*(p2 + p3)) + 2*Om2**2*(p0**3*p1 - 3*p0**2*p2*p3 + p2*p3*(3*p1**2 + p2**2 - p3**2) - p0*p1*(p1**2 + 3*(p2 - p3)*(p2 + p3))) - 4*p0*p1*p2**2*V1**2 - 4*p0**2*p2*p3*V1**2 + 4*p1**2*p2*p3*V1**2 + 4*p0*p1*p3**2*V1**2 - 2*p0**3*p2*V1*V2 + 6*p0*p1**2*p2*V1*V2 - 2*p0*p2**3*V1*V2 + 6*p0**2*p1*p3*V1*V2 - 2*p1**3*p3*V1*V2 + 6*p1*p2**2*p3*V1*V2 + 6*p0*p2*p3**2*V1*V2 - 2*p1*p3**3*V1*V2 + 2*p0**3*p1*V2**2 - 2*p0*p1**3*V2**2 + 2*p0*p1*p2**2*V2**2 + 2*p0**2*p2*p3*V2**2 - 2*p1**2*p2*p3*V2**2 + 2*p2**3*p3*V2**2 - 2*p0*p1*p3**2*V2**2 - 2*p2*p3**3*V2**2 + 6*p0**2*p1*p2*V1*V3 - 2*p1**3*p2*V1*V3 - 2*p1*p2**3*V1*V3 + 2*p0**3*p3*V1*V3 - 6*p0*p1**2*p3*V1*V3 - 6*p0*p2**2*p3*V1*V3 + 6*p1*p2*p3**2*V1*V3 + 2*p0*p3**3*V1*V3 + p0**4*V2*V3 - 6*p0**2*p1**2*V2*V3 + p1**4*V2*V3 - p2**4*V2*V3 + 6*p2**2*p3**2*V2*V3 - p3**4*V2*V3 - 2*p0**3*p1*V3**2 + 2*p0*p1**3*V3**2 + 2*p0*p1*p2**2*V3**2 + 2*p0**2*p2*p3*V3**2 - 2*p1**2*p2*p3*V3**2 - 2*p2**3*p3*V3**2 - 2*p0*p1*p3**2*V3**2 + 2*p2*p3**3*V3**2 + Om3*((p0**2 - p1**2 + p2**2 - p3**2)*V1 - 2*(p1*p2 + p0*p3)*V2)*y3 + Om1*(-2*p0*p1*V2 + 2*p2*p3*V2 - p0**2*V3 - p2**2*V3 + (p1**2 + p3**2)*V3)*y3 + Om2*(2*Om1*(p0*p2*(p0**2 - 3*p1**2 + p2**2) + p1*(-3*p0**2 + p1**2 - 3*p2**2)*p3 - 3*p0*p2*p3**2 + p1*p3**3) + Om3*((p0**2 - p1**2)**2 - (p2**2 - p3**2)**2) + 2*(p0*p1*V1 - p2*p3*V1 + p1*p2*V3 + p0*p3*V3)*y3)

            !rslt(1) = 0.0_8 !Om1*Om3*(p0**2 - p1**2 - p2**2 + p3**2) + 2*(p0*p2 - p1*p3)*(Om3 - V1)*(Om3 + V1) + 2*(p0*p1 + p2*p3)*V1*V2 + (p0**2 - p1**2 - p2**2 + p3**2)*V1*V3
            !rslt(2) = 0.0_8 !2*Om1*Om3*(-(p0*p2) + p1*p3) + Om1**2*(-p0**2 + p1**2 + p2**2 - p3**2) + Om2*(2*Om3*(p0*p1 + p2*p3) + Om2*(-p0**2 + p1**2 + p2**2 - p3**2)) + 2*(-(p0*p2*V1) + p1*p3*V1 + p0*p1*V2 + p2*p3*V2)*V3 + (p0**2 - p1**2 - p2**2 + p3**2)*V3**2
#endif
         end associate
      end function myZnh

      pure function mymatZ(this,q,v,T) result(rslt)
         ! input
         class(myproblem),                 intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         real(8), dimension(:),            intent(in) :: v
         real(8), dimension(:,:),          intent(in) :: T
         ! result
         real(8), dimension(this%sizel, this%sizev)   :: rslt
         !
         ERROR STOP "function mymatZ is not implemented and should not be reached"
         !rslt = 37.77777777777777777777777777777777777777777_8
      end function mymatZ

      pure function mynorm(this, v) result(rslt)
         ! input
         class(myproblem),      intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         ! result
         real(8)                             :: rslt
         !
         ERROR STOP "function mynorm is not implemented and should not be reached"
      end function mynorm

      pure function myCt(this, q, v, t) result(rslt)
         ! input
         class(myproblem),      intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8), dimension(:), intent(in)   :: v
         real(8),               intent(in)   :: t
         ! result
         real(8), dimension(size(v),size(v)) :: rslt
         !
         ERROR STOP "function myCt is not implemented and should not be reached"
      end function myCt

      pure function myKt(this, q, v, vd, t) result(rslt)
         ! input
         class(myproblem),      intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8), dimension(:), intent(in)   :: v
         real(8), dimension(:), intent(in)   :: vd
         real(8),               intent(in)   :: t
         ! result
         real(8), dimension(size(v),size(v)) :: rslt
         !
         ERROR STOP "function myKt is not implemented and should not be reached"
      end function myKt

      pure function myKtl(this, q, v, vd, l, t) result(rslt)
         ! input
         class(myproblem),      intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8), dimension(:), intent(in)   :: v
         real(8), dimension(:), intent(in)   :: vd
         real(8), dimension(:), intent(in)   :: l
         real(8),               intent(in)   :: t
         ! result
         real(8), dimension(size(v),size(v)) :: rslt
         !
         ERROR STOP "function myKtl is not implemented and should not be reached"
      end function myKtl

      pure function modmod(r1,r2) result(rslt)
         ! input
         real(8), intent(in)  :: r1
         real(8), intent(in)  :: r2
         ! result
         real(8)              :: rslt
         !
         rslt = r1 - nint(r1/r2)*r2
      end function modmod

      subroutine myoutputFunction(this,info)
         ! input
         class(myproblem), intent(in)  :: this
         integer,          intent(in)  :: info

         select case (info)
            case (0)
               ! initialization:

               ! Binary output file was opened in main, because 'this'
               ! is intent(in)

               !! opening a file for misc output (DEBUG)
               !open (newunit=126, file="out/testprobrolling_disk_misc")
            case (1)
               if (this%output_t_at == 0 .or. &
                   abs(modmod(this%t,this%t_output_at_multiples_of))  < 1.0e-9_8) then ! TODO: Besserer Wert? DEBUG

                  ! normal output:
                  ! write time to file
                  write (this%out_bin_lun) real(this%t,8)
                  ! write desired values to file
                  if (this%opts%constrained == 1) then
                     if (this%opts%stab2 == 1) then
                        write (this%out_bin_lun) real(this%q,8), real(this%v,8), &
                           real(this%l,8), real(this%lm,8), real(this%lnh,8), &
                           real(this%GL(INTEGRATOR)_phi(this%q),8), real(matmul(this%GL(INTEGRATOR)_B(this%q), this%v),8), &
                           real(matmul(this%GL(INTEGRATOR)_Bnh(this%q), this%v),8)
                     else
                           ERROR STOP "RATTLie does not support index-3"
                     end if
                  else
                     ERROR STOP "example is constrained, but constrained == 0"
                  end if

                  !! write misc output (DEBUG)
                  !write (this%out_misc_lun,*), this%t, this%GL(INTEGRATOR)_stats%newt_steps_curr, this%GL(INTEGRATOR)_stats%newt_steps_max, this%GL(INTEGRATOR)_stats%newt_steps_avg
                  !flush (this%out_misc_lun)
               end if
            case (99)
               ! termination:

               ! Closing binary output file is done in main

            case default
               print *, '[ERROR] Got unsupported info flag ', info, ' in myoutputFunction'
               ERROR STOP 'Got unsupported info flat in myoutputFunction'
         end select
      end subroutine myoutputFunction

      subroutine calculate_jour(this)
         ! input/output
         class(myproblem), intent(inout)     :: this
         !
         ERROR STOP "Calculate_jour not supported"
      end subroutine calculate_jour

      subroutine myinit(this)
         ! input/output
         class(myproblem), intent(inout)        :: this
         ! internal
         real(8)  :: aa(3,3)
         real(8)  :: bb(3,6)
         real(8)  :: rhs(3)
         integer  :: ipiv(3)
         integer  :: info

         ! Set current time point (t0)
         this%t  = this%opts%t0

         ! Set the sizes of q and v
         this%sizeq = 7
         this%sizev = 6

         ! Allocate memory for q, v, vd and a
         if (allocated(this%q)) then
            deallocate(this%q)
         end if
         allocate(this%q(this%sizeq))
         !
         if (allocated(this%v)) then
            deallocate(this%v)
         end if
         allocate(this%v(this%sizev))
         if (allocated(this%p)) then
            deallocate(this%p)
         end if
         allocate(this%p(this%sizev))

         ! Set size of Lagrange multiplier
#ifdef ONLY_NONHOL
         this%sizel = 0
#else
         this%sizel = 1
#endif

         ! allocate Lagrange multipliers
         if (allocated(this%lm)) then
            deallocate(this%lm)
         end if
         allocate(this%lm(this%sizel))
         ! allocate this%l
         if (allocated(this%lp)) then
            deallocate(this%lp)
         end if
         allocate(this%lp(this%sizel))
         ! allocate this%l
         if (allocated(this%l)) then
            deallocate(this%l)
         end if
         allocate(this%l(this%sizel))
         ! Initial values for lambda are calculated by integrator

         ! Set size of nonhol Lagrange multiplier
#ifdef ONLY_NONHOL
         this%sizelnh = 3
#else
         this%sizelnh = 2
#endif

         ! allocate nonhol Lagrange multipliers
         if (allocated(this%lnh)) then
            deallocate(this%lnh)
         end if
         allocate(this%lnh(this%sizelnh))
         ! Initial values for lambda are calculated by integrator


         ! Fill in the initial values
         this%q(1:4) = this%p0
         this%q(5:6) = this%x0(1:2)
         this%v(1:3) = this%Om0

         associate (p0 => this%p0(1), p1 => this%p0(2), p2 => this%p0(3), p3 => this%p0(4))
            this%q(7) = 2*sqrt((p0**2 + p3**2)*(p1**2 + p2**2))
         !   this%v(4:6) = cross([-2*p0*p2+2*p1*p3, 2*p0*p1 + 2*p2*p3,0.0_8],this%Om0)/this%x0(3)
         end associate

         !! Calculate initial velocity
         print *, this%sizel
         bb(1:this%sizel,:) = myB(this, this%q)
         bb(this%sizel+1:this%sizel+this%sizelnh,:) = myBnh(this, this%q)
         aa = bb(:,4:6)
                  call print_matrix(bb, 'bb')
         rhs = -matmul(bb(:,1:3),this%Om0)
         call dgesv(      &! solve the System A*X=B and save the result in B
                     3,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                     1,   &! number of right hand sides (=size(B,2))
                     aa,  &! matrix A
                     3,   &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                     ipiv,&! integer pivot vector; it is not needed
                     rhs, &! matrix B
                     3,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                     info)  ! integer information flag
         this%v(4:6) = rhs

         if (info.ne.0) then
            ERROR STOP "Could not initialize"
         end if

      end subroutine myinit

      subroutine print_m_list(vec)
         ! input
         real(8), intent(in)  :: vec(:)
         ! internal
         integer              :: i,j
         character(len=26)    :: tmpstr
         logical              :: hasE
         !
         write(*,'(A)',advance='no') '{'
         do i=1,size(vec)
            write(tmpstr,'(E26.16E3)') vec(i)
            hasE = .false.
            do j=1,26
               if (tmpstr(j:j)=='E') then
                  write(*,'(A)',advance='no') '*10^('
                  hasE = .true.
               else
                  if (.not. tmpstr(j:j)==' ') write(*,'(A)',advance='no') tmpstr(j:j)
               end if
            end do
            if (hasE) write(*,'(A)',advance='no') ')'
            if (.not. i==size(vec)) write(*,'(A)',advance='no') ','
         end do
         write(*,'(A)',advance='no') '}'
      end subroutine print_m_list

      subroutine print_m_vector(vec, nam)
         ! input
         real(8), intent(in)  :: vec(:)
         character(len=*)     :: nam
         !
         write(*,'(2A)',advance='no') nam, '='
         call print_m_list(vec)
         write(*,'(A)') ';'
      end subroutine print_m_vector

      subroutine print_m_matrix(mat, nam)
         ! input
         real(8), intent(in)  :: mat(:,:)
         character(len=*)     :: nam
         ! internal
         integer              :: i
         !
         write(*,'(2A)',advance='no') nam, '={'
         do i=1,size(mat,1)
            call print_m_list(mat(i,:))
            if (.not. i==size(mat,1)) write(*,'(A)',advance='no') ','
         end do
         write(*,'(A)') '};'
      end subroutine print_m_matrix

end module testprobrolling_disk
