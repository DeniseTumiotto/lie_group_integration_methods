#define GL(x) x

module testprobcrm
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
      integer                 :: n
      !
      integer                 :: n1
      !
      integer                 :: fixed_x0
      real(8)                 :: fixed_x0_position(3)
      !
      integer                 :: fixed_xn
      real(8)                 :: fixed_xN_position(3)
      !
      integer                 :: fixed_p0
      real(8)                 :: fixed_p0_orientation(4)
      !
      integer                 :: fixed_pn
      real(8)                 :: fixed_pn_orientation(4)
      !
      real(8)                 :: ds
      !
      real(8), dimension(3)   :: mI
      real(8)                 :: m
      !
      real(8), dimension(3)   :: CGam
      real(8), dimension(3)   :: CGamd
      !
      real(8), dimension(3)   :: CK
      real(8), dimension(3)   :: CKd
      ! Initial configuration
      real(8), allocatable    :: x0(:,:)
      real(8), allocatable    :: V0(:,:)
      real(8), allocatable    :: p0(:,:)
      real(8), allocatable    :: Om0(:,:)
      !
      integer                 :: output_s_at = 0
      real(8), allocatable    :: output_s(:)
      !! TODO: PRIVATE!
      integer, allocatable    :: output_s_id(:)
      real(8), allocatable    :: output_s_del(:)
      integer, allocatable    :: output_s_id_l(:)
      real(8), allocatable    :: output_s_del_l(:)
      !
      integer                 :: output_t_at = 0
      real(8)                 :: t_output_at_multiples_of = 0.0_8
      !!
      integer                 :: kirchhoff = 0
      !integer                 :: inextensible = 0
      !
      character(len=256)      :: externalfm_name
      integer                 :: externalfm_type
                                 ! type = 1: name = 'gravity'
      real(8)                 :: externalfm_gravity_g ! usually -9.81
                                 ! type = 2: name = 'flying_spaghetti'
      real(8)                 :: externalfm_flying_spaghetti_increasing_time   ! usually 2.5
      real(8)                 :: externalfm_flying_spaghetti_decreasing_time   ! usually 2.5
      real(8)                 :: externalfm_flying_spaghetti_maximum_height    ! usually 200
      real(8)                 :: externalfm_flying_spaghetti_force_factors(3)  ! usually [1/10, 0, 0]
      real(8)                 :: externalfm_flying_spaghetti_moment_factors(3) ! usually [0, -1, -1/2]
                                 ! type = 3: name = 'roll-up'
      real(8)                 :: externalfm_roll_up_factor
   contains

      ! referencing the former deferred procedures
      procedure   :: GL(INTEGRATOR)_M              => myM
      procedure   :: GL(INTEGRATOR)_diag_M         => mydiagM
#if defined(INT_RATTLie) || defined(INT_SHAKELie) || defined(INT_varint4lie)
      procedure   :: GL(INTEGRATOR)_f              => crm_f
#else
      procedure   :: GL(INTEGRATOR)_g              => crm_g
#endif
#if defined(INT_RATTLie)
      procedure   :: GL(INTEGRATOR)_inertial       => crm_inertial
#endif
      procedure   :: GL(INTEGRATOR)_qlpexphDqtilde => myqlpexphDqtilde
      procedure   :: GL(INTEGRATOR)_itlbtvtw       => myitlbtvtw
      procedure   :: GL(INTEGRATOR)_tilde          => myembeddingop
      procedure   :: GL(INTEGRATOR)_Ct             => myCt
      procedure   :: GL(INTEGRATOR)_Kt             => myKt
      procedure   :: GL(INTEGRATOR)_Kt_lambda      => myKtl
#if defined(INT_RATTLie) || defined(INT_SHAKELie) || defined(INT_varint4lie)
      procedure   :: GL(INTEGRATOR)_Tg_inv_T       => my_Tg_inv_T
      procedure   :: GL(INTEGRATOR)_d_Tg_inv_T     => my_d_Tg_inv_T
#endif
      procedure   :: GL(INTEGRATOR)_Tg             => myTg
      procedure   :: GL(INTEGRATOR)_norm           => mynorm
      procedure   :: GL(INTEGRATOR)_outputFunction => myoutputFunction
      procedure   :: GL(INTEGRATOR)_init           => myinit
      procedure   :: GL(INTEGRATOR)_phi            => myphi
      !! DEBUG
      procedure   :: GL(INTEGRATOR)_b              => myb
      !procedure   :: GL(INTEGRATOR)_B              => GL(INTEGRATOR)_num_B
      !! GUBED
      procedure   :: GL(INTEGRATOR)_Z              => myZ
      procedure   :: GL(INTEGRATOR)_matZ           => mymatZ

      ! Own procedures
      procedure   :: calculate_jour      => calculate_jour
      !procedure   :: external_force      => external_force
      !procedure   :: external_moment     => external_moment
   end type myproblem

   ! functions and subroutines that the module contains
   ! here the procedures we need in order to solve the problem are
   ! actually implemented
   contains
      pure subroutine px_from_q(this,p,x,q)
         ! input
         class(myproblem),                intent(in   )  :: this
         real(8), dimension(:),           intent(in   )  :: q
         ! output
         real(8), dimension(4,this%n1),   intent(  out)  :: p
         real(8), dimension(3,this%n1),   intent(  out)  :: x
         ! internal
         integer                                         :: point
         integer                                         :: i
         !
         point = 1
         !
         if (this%fixed_p0 == 1) then
            p(:,1) = this%fixed_p0_orientation
            point = point
         else
            p(:,1) = q(1:4)
            point = point + 4
         end if
         !
         if (this%fixed_x0 == 1) then
            x(:,1) = this%fixed_x0_position
            point = point
         else
            x(:,1) = q(5:7)
            point = point + 3
         end if

         do i = 2,this%n
            p(:,i) = q(point  :point+3)
            x(:,i) = q(point+4:point+6)
            point = point + 7
         end do

         if (this%fixed_pn == 1) then
            p(:,this%n1) = this%fixed_pn_orientation
         else
            p(:,this%n1) = q(point  :point+3)
            point = point + 4
         end if
         !
         if (this%fixed_xn == 1) then
            x(:,this%n1) = this%fixed_xn_position
         else
            x(:,this%n1) = q(point  :point+2)
         end if
      end subroutine px_from_q

      pure subroutine OmV_from_v(this,Om,cV,v)
         ! input
         class(myproblem),                intent(in   )  :: this
         real(8), dimension(:),           intent(in   )  :: v
         ! output
         real(8), dimension(3,this%n1),   intent(  out)  :: Om
         real(8), dimension(3,this%n1),   intent(  out)  :: cV
         ! internal
         integer                                         :: point
         integer                                         :: i
         !
         point = 1
         !
         if (this%fixed_p0 == 1) then
            Om(:,1) = 0.0_8
            point = point
         else
            Om(:,1) = v(1:3)
            point = point + 3
         end if
         !
         if (this%fixed_x0 == 1) then
            cV(:,1) = 0.0_8 ! Derivative of the constant this%fixed_x0_position
            point = point
         else
            cV(:,1) = v(4:6)
            point = point + 3
         end if

         do i = 2,this%n
            Om(:,i) = v(point  :point+2)
            cV(:,i) = v(point+3:point+5)
            point = point + 6
         end do

         if (this%fixed_pn == 1) then
            Om(:,this%n1) = 0.0_8
         else
            Om(:,this%n1) = v(point  :point+2)
            point = point + 3
         end if
         if (this%fixed_xn == 1) then
            cV(:,this%n1) = 0.0_8 ! Derivative of the constant this%fixed_xn_position
         else
            cV(:,this%n1) = v(point  :point+2)
         end if
      end subroutine OmV_from_v

      pure function lin_interpol_S3R3(t,q1,q2) result(rslt)
         ! input
         real(8), intent(in)  :: t
         real(8), intent(in)  :: q1(7)
         real(8), intent(in)  :: q2(7)
         ! result
         real(8)              :: rslt(7)
         !
         rslt = lp_s3sdr3(q1, expt_s3sdr3(t*logt_s3sdr3(lp_s3sdr3(inv_s3sdr3(q1),q2))))
      end function lin_interpol_S3R3

      pure function lin_interpol_Rn(t,v1,v2) result(rslt)
         ! input
         real(8), intent(in)  :: t
         real(8), intent(in)  :: v1(:)
         real(8), intent(in)  :: v2(:)
         ! result
         real(8)              :: rslt(size(v1))
         !
         rslt = (1.0_8-t)*v1 + t*v2
      end function lin_interpol_Rn

!      pure function quat_lin_interpol(t,q1,q2) result(rslt)
!         ! input
!         real(8), intent(in)  :: t
!         real(8), intent(in)  :: q1(0:3)
!         real(8), intent(in)  :: q2(0:3)
!         ! result
!         real(8)              :: rslt(0:3)
!         !
!         rslt = qp(q1,exptlogqS3(t,kqp(q1,q2)))
!      end function quat_lin_interpol

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
         ! internal
         integer                             :: i
         integer                             :: indq, indv
         !
         if (this%fixed_x0 == 1) then
            if (this%fixed_p0 == 1) then
               indq = 1
               indv = 1
            else
               rslt(1:4) = lp_s3(q(1:4),expt_s3(h*Dq(1:3)))
               indq = 5
               indv = 4
            end if
         else
            if (this%fixed_p0 == 1) then
               ERROR STOP "free x, fixed p not implemented"
            else
               rslt(1:7) = lp_s3sdr3(q(1:7),expt_s3sdr3(h*Dq(1:6)))
               indq = 8
               indv = 7
            end if
         end if

         do i=2,this%n
            associate (qk  => q(indq  :indq+6),    &
                       vk  => h*Dq(indv :indv+5)  )
               rslt(indq  :indq+6) = lp_s3sdr3(qk, expt_s3sdr3(vk))
            end associate
            indq = indq + 7
            indv = indv + 6
         end do

         if (this%fixed_xn == 1) then
            if (this%fixed_pn /= 1) then
               associate (pk  => q(indq  :indq+3),    &
                          Omk  => h*Dq(indv :indv+2)  )
                  rslt(indq  :indq+3) = lp_s3(pk, expt_s3(Omk))
               end associate
            end if
         else
            if (this%fixed_pn == 1) then
               ERROR STOP "free x, fixed p not implemented"
            else
               associate (qk  => q(indq  :indq+6),    &
                          vk  => h*Dq(indv :indv+5)  )
                  rslt(indq  :indq+6) = lp_s3sdr3(qk, expt_s3sdr3(vk))
               end associate
            end if
         end if

         !! DEBUG
         !if (abs(norm2(rslt(1:4)) - 1.0_8) > 1.0e-10_8) then
         !   print *, 'norm=', norm2(rslt(1:4))
         !   call print_vector(rslt,'rslt')
         !   pause
         !end if
         !! GUBED
      end function myqlpexphDqtilde

      pure function myitlbtvtw(this, v, w) result(rslt)
         ! input
         class(myproblem),      intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         real(8), dimension(:), intent(in)   :: w
         ! result
         real(8), dimension(size(v))         :: rslt
         ! internal
         integer                             :: i, r
         !
         !ERROR STOP "function myitlbtvtw is not implemented and should not be reached"

         if (this%fixed_x0 == 1) then
            if (this%fixed_p0 == 1) then
               r = 6
            else
               rslt(1:3) = cross(v(1:3),w(1:3))
               r = 3
            end if
         else
            if (this%fixed_p0 == 1) then
               ERROR STOP "free x, fixed p not implemented"
            else
               rslt(1:6) = lie_bracket_s3sdr3(v(1:6),w(1:6))
               r = 0
            end if
         end if

         do i=1,this%n-1
            rslt(6*i+1-r:6*i+6-r) = lie_bracket_s3sdr3(v(6*i+1-r:6*i+6-r),w(6*i+1-r:6*i+6-r))
         end do

         i = this%n
         if (this%fixed_xn == 1) then
            if (this%fixed_pn /= 1) then
               rslt(6*i+1-r:6*i+3-r) = cross(v(6*i+1-r:6*i+3-r),w(6*i+1-r:6*i+3-r))
            end if
         else
            if (this%fixed_pn == 1) then
               ERROR STOP "free x, fixed p not implemented"
            end if
            rslt(6*i+1-r:6*i+6-r) = lie_bracket_s3sdr3(v(6*i+1-r:6*i+6-r),w(6*i+1-r:6*i+6-r))
         end if
      end function myitlbtvtw

!      ! the following function gives the i-th row of the result of the
!      ! euler map of a quaternion p that is defined by
!      !     p*v*konj(p) = R(p) v
!      pure function pvkptoRi(p,i) result(rslt)
!         ! input
!         real(8), dimension(1:4), intent(in) :: p
!         integer,                 intent(in) :: i
!         ! result
!         real(8), dimension(3)               :: rslt
!         ! internal variables
!         real(8), dimension(3)               :: ei
!         !
!         ! probably wrong TODO
!         !if (i == 1) then
!         !   rslt = [ p(1)**2 + p(2)**2 - p(3)**2 - p(4)**2, &
!         !           -2*p(1)*p(4) + 2*p(2)*p(3),             &
!         !            2*p(1)*p(3) + 2*p(2)*p(4)               ]
!         !elseif (i == 2) then
!         !   rslt = [ 2*p(1)*p(4) + 2*p(2)*p(3),             &
!         !            p(1)**2 - p(2)**2 + p(3)**2 - p(4)**2, &
!         !           -2*p(1)*p(2) + 2*p(3)*p(4)               ]
!         !elseif (i == 3) then
!         !   rslt = [-2*p(1)*p(3) + 2*p(2)*p(4),             &
!         !            2*p(1)*p(2) + 2*p(3)*p(4),             &
!         !            p(1)**2 - p(2)**2 - p(3)**2 + p(4)**2   ]
!         !end if
!         !
!         ! probably right TODO
!         if (i == 1) then
!            rslt = [ p(1)**2 + p(2)**2 - p(3)**2 - p(4)**2, &
!                     2*p(1)*p(4) + 2*p(2)*p(3),             &
!                    -2*p(1)*p(3) + 2*p(2)*p(4)               ]
!         elseif (i == 2) then
!            rslt = [-2*p(1)*p(4) + 2*p(2)*p(3),             &
!                     p(1)**2 - p(2)**2 + p(3)**2 - p(4)**2, &
!                     2*p(1)*p(2) + 2*p(3)*p(4)               ]
!         elseif (i == 3) then
!            rslt = [ 2*p(1)*p(3) + 2*p(2)*p(4),             &
!                    -2*p(1)*p(2) + 2*p(3)*p(4),             &
!                     p(1)**2 - p(2)**2 - p(3)**2 + p(4)**2   ]
!         end if
!      end function pvkptoRi

      pure function myTg(this, h, dq) result(rslt)
         ! input
         class(myproblem),      intent(in)      :: this
         real(8),               intent(in)      :: h
         real(8), dimension(:), intent(in)      :: dq
         ! result
         real(8), dimension(size(dq),size(dq))  :: rslt
         ! internal
         integer                                :: i
         integer                                :: ind
         !
         rslt = 0.0_8

         if (this%fixed_x0 == 1) then
            if (this%fixed_p0 == 1) then
               ind = 1
            else
               rslt(1:3,1:3) = tan_op_s3(h*Dq(1:3))
               ind = 4
            end if
         else
            if (this%fixed_p0 == 1) then
               ERROR STOP "free x, fixed p not implemented"
            else
               rslt(1:6,1:6) = tan_op_s3sdr3(h*Dq(1:6))
               ind = 7
            end if
         end if

         do i=2,this%n
            rslt(ind  :ind+5,ind  :ind+5) = tan_op_s3sdr3(h*Dq(ind  :ind+5))
            ind = ind + 6
         end do


         if (this%fixed_xn == 1) then
            if (this%fixed_pn /= 1) then
               rslt(ind  :ind+2,ind  :ind+2) = tan_op_s3(h*Dq(ind  :ind+2))
            end if
         else
            if (this%fixed_pn == 1) then
               ERROR STOP "free x, fixed p not implemented"
            else
               rslt(ind  :ind+5,ind  :ind+5) = tan_op_s3sdr3(h*Dq(ind  :ind+5))
            end if
         end if
      end function myTg

      pure function my_Tg_inv_T(this, dq) result(rslt)
         ! input
         class(myproblem),      intent(in)         :: this
         real(8), dimension(:), intent(in)         :: dq
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
         ! internal
         integer                                   :: i, ind
         !
         rslt = 0.0_8

         if (this%fixed_x0 == 1) then
            if (this%fixed_p0 == 1) then
               ind = 1
            else
               rslt(1:3,1:3) = tan_tr_inv_s3(dq(1:3))
               ind = 4
            end if
         else
            if (this%fixed_p0 == 1) then
               ERROR STOP "free x, fixed p not implemented"
            else
               rslt(1:6,1:6) = tan_tr_inv_s3sdr3(dq(1:6))
               ind = 7
            end if
         end if

         do i=2,this%n
            rslt(ind  :ind+5,ind  :ind+5) = tan_tr_inv_s3sdr3(dq(ind  :ind+5))
            ind = ind + 6
         end do


         if (this%fixed_xn == 1) then
            if (this%fixed_pn /= 1) then
               rslt(ind  :ind+2,ind  :ind+2) = tan_tr_inv_s3(dq(ind  :ind+2))
            end if
         else
            if (this%fixed_pn == 1) then
               ERROR STOP "free x, fixed p not implemented"
            else
               rslt(ind  :ind+5,ind  :ind+5) = tan_tr_inv_s3sdr3(dq(ind  :ind+5))
            end if
         end if
      end function my_Tg_inv_T

      pure function my_d_Tg_inv_T(this, v, w) result(rslt)
         ! input
         class(myproblem),      intent(in)         :: this
         real(8), dimension(:), intent(in)         :: v, w
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
         ! internal
         integer                                   :: i, ind
         !
         rslt = 0.0_8

         if (this%fixed_x0 == 1) then
            if (this%fixed_p0 == 1) then
               ind = 1
            else
               rslt(1:3,1:3) = d_tan_tr_inv_s3(v(1:3),w(1:3))
               ind = 4
            end if
         else
            if (this%fixed_p0 == 1) then
               ERROR STOP "free x, fixed p not implemented"
            else
               rslt(1:6,1:6) = d_tan_tr_inv_s3sdr3(v(1:6),w(1:6))
               ind = 7
            end if
         end if

         do i=2,this%n
            rslt(ind  :ind+5,ind  :ind+5) = d_tan_tr_inv_s3sdr3(v(ind  :ind+5), w(ind  :ind+5))
            ind = ind + 6
         end do


         if (this%fixed_xn == 1) then
            if (this%fixed_pn /= 1) then
               rslt(ind  :ind+2,ind  :ind+2) = d_tan_tr_inv_s3(v(ind  :ind+2), w(ind  :ind+2))
            end if
         else
            if (this%fixed_pn == 1) then
               ERROR STOP "free x, fixed p not implemented"
            else
               rslt(ind  :ind+5,ind  :ind+5) = d_tan_tr_inv_s3sdr3(v(ind  :ind+5), w(ind  :ind+5))
            end if
         end if
      end function my_d_Tg_inv_T

      pure function myM(this, q) result(rslt)
         ! input
         class(myproblem),      intent(in)               :: this
         real(8), dimension(:), intent(in)               :: q
         ! result
         real(8), dimension(this%sizev,this%sizev)       :: rslt
         ! internal variables
         integer                                         :: i
         integer                                         :: j
         integer                                         :: ind
         !
         rslt = 0.0_8
         !
         ind = 1
         !
         if (.not. this%fixed_p0 == 1) then
            forall (j=0:2)
               rslt(ind+j,ind+j) = this%mI(1+j)/2
            end forall
            ind = ind + 3
         end if
         !
         if (.not. this%fixed_x0 == 1) then
            forall (j=0:2)
               rslt(ind+j,ind+j) = this%m/2
            end forall
            ind = ind + 3
         end if

         do i=2,this%n
            forall (j=0:2)
               rslt(ind+j,ind+j) = this%mI(1+j)
            end forall
            forall (j=ind+3:ind+5)
               rslt(j,j) = this%m
            end forall
            ind = ind + 6
         end do

         if (.not. this%fixed_p0 == 1) then
            forall (j=0:2)
               rslt(ind+j,ind+j) = this%mI(1+j)/2
            end forall
            ind = ind + 3
         end if
         !
         if (.not. this%fixed_xn == 1) then
            forall (j=ind+3:ind+5)
               rslt(j,j) = this%m/2
            end forall
         end if
         !
      end function myM

      pure function mydiagM(this, q) result(rslt)
         ! input
         class(myproblem),      intent(in)               :: this
         real(8), dimension(:), intent(in)               :: q
         ! result
         real(8), dimension(this%sizev)                  :: rslt
         ! internal variables
         integer                                         :: i
         integer                                         :: j
         integer                                         :: ind
         !
         ind = 1
         !
         if (.not. this%fixed_p0 == 1) then
            forall (j=1:3)
               rslt(j) = this%mI(j)/2
            end forall
            ind = ind + 3
         end if
         !
         if (.not. this%fixed_x0 == 1) then
            forall (i=0:2)
               rslt(ind+i) = this%m/2
            end forall
            ind = ind + 3
         end if

         do i=2,this%n
            forall (j=1:3)
               rslt(ind+j-1) = this%mI(j)
            end forall
            forall (j=ind+3:ind+5)
               rslt(j) = this%m
            end forall
            ind = ind + 6
         end do

         if (.not. this%fixed_pn == 1) then
            forall (j=1:3)
               rslt(ind+j-1) = this%mI(j)/2
            end forall
            ind = ind + 3
         end if
         !
         if (.not. this%fixed_xn == 1) then
            forall (j=0:2)
               rslt(ind+j) = this%m/2
            end forall
         end if
         !
      end function mydiagM

      pure function fill_initial_q(this, p, x) result(rslt)
         ! input
         class(myproblem), intent(in)  :: this
         real(8),          intent(in)  :: p(:,:)
         real(8),          intent(in)  :: x(:,:)
         ! result
         real(8)                       :: rslt(this%sizeq)
         ! internal
         integer                       :: i
         integer                       :: r
         !
         ! Fill in first discretization point and set offset r
         if (this%fixed_p0 == 1) then
            r = 4
         else
            rslt(1:4) = p(:,1)
            r = 0
         end if
         if (this%fixed_x0 == 1) then
            r = r + 3
         else
            rslt(5-r:7-r) = x(:,1)
         end if

         ! Fill in all inner discretization points
         forall (i=1:this%n-1)
            rslt(7*i-r+1:7*i-r+4) = p(:,i+1)
            rslt(7*i-r+5:7*i-r+7) = x(:,i+1)
         end forall

         ! Fill in last discretization point
         i = this%n
         if (.not. this%fixed_pn == 1) then
            rslt(7*i-r+1:7*i-r+4) = p(:,i+1)
         else
            r = r + 4
         end if
         if (.not. this%fixed_xn == 1) then
            rslt(7*i-r+5:7*i-r+7) = x(:,i+1)
         end if
      end function fill_initial_q

      pure function fill_initial_v(this, Om, V) result(rslt)
         ! input
         class(myproblem), intent(in)  :: this
         real(8),          intent(in)  :: Om(:,:)
         real(8),          intent(in)  :: V(:,:)
         ! result
         real(8)                       :: rslt(this%sizev)
         ! internal
         integer                       :: i
         integer                       :: r
         !
         ! Fill in first discretization point and set offset r
         if (this%fixed_p0 == 1) then
            r = 3
         else
            rslt(1:3) = Om(:,1)
            r = 0
         end if
         if (this%fixed_x0 == 1) then
            r = r + 3
         else
            rslt(4:6) = V(:,1)
         end if

         ! Fill in all inner discretization points
         forall (i=1:this%n-1)
            rslt(6*i-r+1:6*i-r+3) = Om(:,i+1)
            rslt(6*i-r+4:6*i-r+6) = V(:,i+1)
         end forall

         ! Fill in last discretization point
         i = this%n
         if (.not. this%fixed_pn == 1) then
            rslt(6*i-r+1:6*i-r+3) = Om(:,i+1)
         else
            r = r + 3
         end if
         if (.not. this%fixed_xn == 1) then
            rslt(6*i-r+4:6*i-r+6) = V(:,i+1)
         end if
      end function fill_initial_v

#if defined(INT_RATTLie) || defined(INT_SHAKELie) || defined(INT_varint4lie)
      pure function crm_f(this, q, v, t) result(rslt)
#else
      pure function crm_g(this, q, v, t) result(rslt)
#endif
         ! input
         class(myproblem),      intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8), dimension(:), intent(in)   :: v
         real(8),               intent(in)   :: t
         ! result
         real(8), dimension(size(v))         :: rslt
         ! internal variables
         integer                             :: i ! for counting
         integer                             :: r, rq
         real(8)                             :: tmpv(6)
         real(8)                             :: tmpq(7)
         real(8), dimension(6*this%n)        :: w
         real(8), dimension(6*this%n)        :: wd
         real(8), dimension(6,this%n1)       :: force_moment

         ! For debugging purposes, define a custom for-block
#define for(i,f,l) forall(i=f:l)
#define endfor end forall
#ifdef pure
#undef for
#undef endfor
#define for(i,f,l) do i=f,l
#define endfor end do
#endif

         !! DEBUG
         !if (.false.) then
         !   print*, '**************************************************************************************'
         !   print*, '**************************************************************************************'
         !   print*, 't=', t
         !   call print_vector(q,'q')
         !   call print_vector(v,'v')
         !end if
         !! GUBED

         ! Set negative offset r in case the node x0 is fixed
         if (this%fixed_x0 == 1) then
            r = 3
         else
            r = 0
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! Calculate the w_{k+1/2} an wd_{k+1/2}
         w = calculate_w(this, q)
         wd = calculate_wd(this, v, w)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! Calculate right hand sides

         !! For external forces and moments
         force_moment = external_force_and_moment(this, t)

         associate (dsscrM => [this%mI(1), this%mI(2), this%mI(3), this%m, this%m, this%m], &
                    scrC  => [this%CK(1), this%CK(2), this%CK(3), this%CGam(1), this%CGam(2), this%CGam(3)], &
                    scrCd => [this%CKd(1), this%CKd(2), this%CKd(3), this%CGamd(1), this%CGamd(2), this%CGamd(3)], &
                    majE => [0.0_8, 0.0_8, 0.0_8,    0.0_8, 0.0_8, 1.0_8], &
                    ds => this%ds )

            ! Right hand side for v_0
            if (this%fixed_x0 == 1) then
               if (this%fixed_p0 == 1) then
                  tmpv = 0.0_8
                  r = 6
                  rq = 7
               else
                  tmpv(1:3) = v(1:3)
                  tmpv(4:6) = 0.0_8
                  r = 3
                  rq = 3
               end if
            else
               if (this%fixed_p0 == 1) then
                  ERROR STOP "free x, fixed p not implemented"
               else
                  tmpv = v(1:6)
                  r = 0
                  rq = 0
               end if
            end if
            tmpv = &
#if !defined(INT_RATTLie) && !defined(INT_SHAKELie) && !defined(INT_varint4lie)
               hat_tr_mult_s3sdr3(tmpv, dsscrM*tmpv/2)  &
#endif
               + tan_tr_inv_mult_s3sdr3(-ds*w(1:6), scrC*(w(1:6) - majE) + 2*scrCd*wd(1:6))
            if (this%fixed_x0 == 1) then
               if (this%fixed_p0 == 1) then
                  ! nothing
               else
                  rslt(1:3) = tmpv(1:3) + apply_conj_quat(q(1:4), force_moment(1:3,1))
               end if
            else
               if (this%fixed_p0 == 1) then
                  ERROR STOP "free x, fixed p not implemented"
               else
                  rslt(1:3) = tmpv(1:3)  + apply_conj_quat(q(1:4), force_moment(1:3,1))
                  rslt(4:6) = tmpv(4:6)  + apply_conj_quat(q(1:4), force_moment(4:6,1))
               end if
            end if


            ! Right hand side for v_k
            for(i,1,this%n-1)
               rslt(6*i-r+1:6*i-r+6) = &
#if !defined(INT_RATTLie) && !defined(INT_SHAKELie) && !defined(INT_varint4lie)
                  hat_tr_mult_s3sdr3(v(6*i-r+1:6*i-r+6), dsscrM*v(6*i-r+1:6*i-r+6))   &
#endif
                   + tan_tr_inv_mult_s3sdr3(-ds*w(6*i+1:6*i+6), scrC*(w(6*i+1:6*i+6) - majE) + 2*scrCd*wd(6*i+1:6*i+6)) &
                   - tan_tr_inv_mult_s3sdr3( ds*w(6*i-5:6*i  ), scrC*(w(6*i-5:6*i  ) - majE) + 2*scrCd*wd(6*i-5:6*i  ))
               rslt(6*i-r+1:6*i-r+3) = rslt(6*i-r+1:6*i-r+3) + apply_conj_quat(q(7*i-rq+1:7*i-rq+4), force_moment(1:3,i+1))
               rslt(6*i-r+4:6*i-r+6) = rslt(6*i-r+4:6*i-r+6) + apply_conj_quat(q(7*i-rq+1:7*i-rq+4), force_moment(4:6,i+1))
            endfor

            ! Right hand side for v_n
            i = this%n
            if (this%fixed_xn == 1) then
               if (this%fixed_pn == 1) then
                  tmpv = 0.0_8
               else
                  tmpv(1:3) = v(6*i-r+1:6*i-r+3)
                  tmpv(4:6) = 0.0_8
               end if
            else
               if (this%fixed_pn == 1) then
                  ERROR STOP "free x, fixed p not implemented"
               else
                  tmpv = v(6*i-r+1:6*i-r+6)
               end if
            end if
            tmpv = &
#if !defined(INT_RATTLie) && !defined(INT_SHAKELie) && !defined(INT_varint4lie)
               hat_tr_mult_s3sdr3(tmpv, dsscrM*tmpv/2)   &
#endif
                - tan_tr_inv_mult_s3sdr3( ds*w(6*i-5:6*i ), scrC*(w(6*i-5:6*i  ) - majE) + 2*scrCd*wd(6*i-5:6*i  ))
            if (this%fixed_xn == 1) then
               if (this%fixed_pn == 1) then
                  ! nothing
               else
                  rslt(6*i-r+1:6*i-r+3) = tmpv(1:3) + apply_conj_quat(q(7*i-rq+1:7*i-rq+4), force_moment(1:3,i+1))
               end if
            else
               if (this%fixed_pn == 1) then
                  ERROR STOP "free x, fixed p not implemented"
               else
                  rslt(6*i-r+1:6*i-r+3) = tmpv(1:3) + apply_conj_quat(q(7*i-rq+1:7*i-rq+4), force_moment(1:3,i+1))
                  rslt(6*i-r+4:6*i-r+6) = tmpv(4:6) + apply_conj_quat(q(7*i-rq+1:7*i-rq+4), force_moment(4:6,i+1))
               end if
            end if

            !! DEBUG
            !call print_vector(dsscrM,'dsscrM')
            !call print_vector(scrC,'scrC')
            !! GUBED

         ! DEBUG end associate

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! Add external forces and moments



         !force_moment = 0.0_8
         !force_moment(4,:) = 1.0_8

         !do i=0,this%n
         !   !force_moment(6,i+1) = 10*8*(real(i,8)/this%n - 0.5_8)**3
         !   associate (ss => real(i,8)/this%n)
         !      force_moment(3+1,i+1) = ss
         !   end associate
         !end do
         !force_moment([1,2,3],1) = 10.0_8
         !force_moment([4,5,6],1) = 1.0_8
         !force_moment(4:6,1+4) = apply_quat( q(7*4+1:7*4+4),[0.0_8, 0.0_8, 5.0e-3_8])
         !force_moment(4,:) = 1.0_8

         !! DEBUG
         !!if (.true.) then
         !if (.false.) then
         !!if (t>=2.0_8) then
         !!if (norm2(rslt) > 100.0_8) then
         !   print*, '**************************************************************************************'
         !   print*, '**************************************************************************************'
         !   print*, 't=', t,';'
         !   !call print_vector(q,'q')
         !   !call print_vector(v,'v')
         !   !call print_vector(w,'w')
         !   !call print_matrix(force_moment,'force_moment')
         !   call print_m_vector(q,'qq')
         !   call print_m_vector(v,'vv')
         !   call print_m_vector(w,'ww')
         !   call print_m_vector(rslt,'rslt')
         !   i = 1
         !   call print_m_vector(  hat_tr_mult_s3sdr3(v(6*i-r+1:6*i-r+6), dsscrM*v(6*i-r+1:6*i-r+6)), 'hatvMv1')
         !   call print_m_vector(  scrC*(w(6*i+1:6*i+6) - majE), 'Cw1')
         !   call print_m_vector(  tan_tr_inv_mult_s3sdr3(-ds*w(6*i+1:6*i+6), scrC*(w(6*i+1:6*i+6) - majE)), 'TTmwCw1')
         !   call print_m_vector(- tan_tr_inv_mult_s3sdr3( ds*w(6*i-5:6*i  ), scrC*(w(6*i-5:6*i  ) - majE)), 'mTTwCW0')
         !   print *, 'ds =', this%ds
         !   print *,'nrmv =', norm2(v)
         !   print *,'nrmrslt =', norm2(rslt)
         !   !if (norm2(rslt)>=100.0_8) pause
         !   pause
         !end if
         !! GUBED


         !force_moment(:,1) = 2*force_moment(:,1)
         !force_moment(:,this%n1) = 2*force_moment(:,this%n1)
         !if (t<1.0_8) then
            !for(i,0,this%n)
            !   rslt(6*i+1:6*i+3) = rslt(6*i+1:6*i+3) + apply_conj_quat(q(7*i+1:7*i+4), force_moment(1:3,i+1))
            !   rslt(6*i+4:6*i+6) = rslt(6*i+4:6*i+6) + apply_conj_quat(q(7*i+1:7*i+4), force_moment(4:6,i+1))
            !endfor
         !end if

         ! DEBUG
         end associate

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! Don't forget the minus
         rslt = -rslt

         ! Undefine custom for blocks
#undef for
#undef endfor

#if defined(INT_RATTLie) || defined(INT_SHAKELie) || defined(INT_varint4lie)
      end function crm_f
#else
      end function crm_g
#endif

#if defined(INT_RATTLie)
      pure function crm_inertial(this, v) result(rslt)
         ! input
         class(myproblem), intent(in)  :: this
         real(8)         , intent(in)  :: v(:)
         ! result
         real(8)                       :: rslt(size(v))
         ! internal
         integer                       :: i, r
         real(8)                       :: tmpv(6)

         !
#define for(i,f,l) forall(i=f:l)
#define endfor end forall
#ifdef pure
#undef for
#undef endfor
#define for(i,f,l) do i=f,l
#define endfor end do
#endif

                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !! Calculate inertial forces

                  associate (dsscrM => [this%mI(1), this%mI(2), this%mI(3), this%m, this%m, this%m])

                     ! Right hand side for v_0
                     if (this%fixed_x0 == 1) then
                        if (this%fixed_p0 == 1) then
                           r = 6
                        else
                           tmpv(1:3) = v(1:3)
                           tmpv(4:6) = 0.0_8
                           tmpv = hat_tr_mult_s3sdr3(tmpv, dsscrM*tmpv/2)
                           rslt(1:3) = tmpv(1:3)
                           r = 3
                        end if
                     else
                        if (this%fixed_p0 == 1) then
                           ERROR STOP "free x, fixed p not implemented"
                        else
                           rslt(1:6) = hat_tr_mult_s3sdr3(v(1:6), dsscrM*v(1:6)/2)
                           r = 0
                        end if
                     end if

                     ! Right hand side for v_k
                     for(i,1,this%n-1)
                        rslt(6*i-r+1:6*i-r+6) = hat_tr_mult_s3sdr3(v(6*i-r+1:6*i-r+6), dsscrM*v(6*i-r+1:6*i-r+6))
                     endfor

                     ! Right hand side for v_n
                     i = this%n
                     if (this%fixed_xn == 1) then
                        if (this%fixed_pn == 1) then
                           ! nothing
                        else
                           tmpv(1:3) = v(6*i-r+1:6*i-r+3)
                           tmpv(4:6) = 0.0_8
                           tmpv = hat_tr_mult_s3sdr3(tmpv, dsscrM*tmpv/2)
                           rslt(6*i-r+1:6*i-r+3) = tmpv(1:3)
                        end if
                     else
                        if (this%fixed_pn == 1) then
                           ERROR STOP "free x, fixed p not implemented"
                        else
                           rslt(6*i-r+1:6*i-r+6) = hat_tr_mult_s3sdr3(v(6*i-r+1:6*i-r+6), dsscrM*v(6*i-r+1:6*i-r+6)/2)
                        end if
                     end if
                  end associate
                  ! Undefine custom for blocks
#undef for
#undef endfor
      end function crm_inertial
#endif

      pure function calculate_w(this, q) result(w)
         ! input
         class(myproblem), intent(in)  :: this
         real(8),          intent(in)  :: q(:)
         ! result
         real(8)                       :: w(6*this%n)
         ! internal
         integer                       :: i
         integer                       :: r
         real(8)                       :: tmpq(7)
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! Calculate the w_{k+1/2}

         ! Calculate w_{1/2}
         if (this%fixed_x0 == 1) then
            if (this%fixed_p0 == 1) then
               tmpq(1:4) = this%fixed_p0_orientation
               r = 7
            else
               tmpq(1:4) = q(1:4)
               r = 3
            end if
            tmpq(5:7) = this%fixed_x0_position
         else
            if (this%fixed_p0 == 1) then
               ERROR STOP "free x, fixed p not implemented"
            else
               tmpq = q(1:7)
               r = 0
            end if
         end if
         !
         i = 0
#ifdef LINEAR_W
         tmpq = lp_s3sdr3(inv_s3sdr3(tmpq),q(7*i-r+8:7*i-r+14))
         w(6*i+1:6*i+6) = 2/this%ds * tmpq(2:7)
#else
         w(6*i+1:6*i+6) = 1/this%ds &
            * logt_s3sdr3(lp_s3sdr3(inv_s3sdr3(tmpq),q(7*i-r+8:7*i-r+14)))
#endif

         ! For debugging purposes, define a custom for-block
#define for(i,f,l) forall(i=f:l)
#define endfor end forall
#ifdef pure
#undef for
#undef endfor
#define for(i,f,l) do i=f,l
#define endfor end do
#endif
         ! Calculate w_{k+1/2}
         for(i,1,this%n-2)
#ifdef LINEAR_W
            tmpq = lp_s3sdr3(inv_s3sdr3(q(7*i-r+1:7*i-r+7)),q(7*i-r+8:7*i-r+14))
            w(6*i+1:6*i+6) = 2/this%ds * tmpq(2:7)
#else
            w(6*i+1:6*i+6) = 1/this%ds &
               * logt_s3sdr3(lp_s3sdr3(inv_s3sdr3(q(7*i-r+1:7*i-r+7)),q(7*i-r+8:7*i-r+14)))
#endif
         endfor
#undef for
#undef endfor

         ! Calculate w_{n-1/2}
         i = this%n-1
         if (this%fixed_xn == 1) then
            if (this%fixed_pn == 1) then
               tmpq(1:4) = this%fixed_pn_orientation
            else
               tmpq(1:4) = q(7*i-r+8:7*i-r+11)
            end if
            tmpq(5:7) = this%fixed_xn_position
         else
            if (this%fixed_pn == 1) then
               ERROR STOP "free x, fixed p not implemented"
            else
               tmpq = q(7*i-r+8:7*i-r+14)
            end if
         end if
         !
#ifdef LINEAR_W
         tmpq = lp_s3sdr3(inv_s3sdr3(q(7*i-r+1:7*i-r+7)),tmpq)
         w(6*i+1:6*i+6) = 2/this%ds * tmpq(2:7)
#else
         w(6*i+1:6*i+6) = 1/this%ds &
            * logt_s3sdr3(lp_s3sdr3(inv_s3sdr3(q(7*i-r+1:7*i-r+7)),tmpq))
#endif

      end function calculate_w

      pure function calculate_wd(this, v, w) result(wd)
         ! input
         class(myproblem), intent(in)  :: this
         real(8),          intent(in)  :: v(:)
         real(8),          intent(in)  :: w(:)
         ! result
         real(8)                       :: wd(6*this%n)
         ! internal
         integer                       :: i
         integer                       :: r
         real(8)                       :: tmpv(6)
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! Calculate the w_{k+1/2}

         ! Calculate w_{1/2}
         if (this%fixed_x0 == 1) then
            if (this%fixed_p0 == 1) then
               tmpv(1:3) = 0.0_8
               r = 6
            else
               tmpv(1:3) = v(1:3)
               r = 3
            end if
            tmpv(4:6) = 0.0_8
         else
            if (this%fixed_p0 == 1) then
               ERROR STOP "free x, fixed p not implemented"
            else
               tmpv = v(1:6)
               r = 0
            end if
         end if
         !
         i = 0
         wd(6*i+1:6*i+6) = 1/this%ds &
            * (  tan_op_inv_mult_s3sdr3(this%ds*w(1:6), v(6*i-r+7:6*i-r+12)) &
               - tan_op_inv_mult_s3sdr3(-this%ds*w(1:6), tmpv) )


         ! For debugging purposes, define a custom for-block
#define for(i,f,l) forall(i=f:l)
#define endfor end forall
#ifdef pure
#undef for
#undef endfor
#define for(i,f,l) do i=f,l
#define endfor end do
#endif
         ! Calculate wd_{k+1/2}
         for(i,1,this%n-2)
            wd(6*i+1:6*i+6) = 1/this%ds &
               * (  tan_op_inv_mult_s3sdr3(this%ds*w(6*i+1:6*i+6), v(6*i-r+7:6*i-r+12)) &
                  - tan_op_inv_mult_s3sdr3(-this%ds*w(6*i+1:6*i+6), v(6*i-r+1:6*i-r+6)) )
         endfor
#undef for
#undef endfor

         ! Calculate w_{n-1/2}
         i = this%n-1
         if (this%fixed_xn == 1) then
            if (this%fixed_pn == 1) then
               tmpv(1:3) = 0.0_8
            else
               tmpv(1:3) = v(6*i-r+7:6*i-r+9)
            end if
            tmpv(4:6) = 0.0_8
         else
            if (this%fixed_pn == 1) then
               ERROR STOP "free x, fixed p not implemented"
            else
               tmpv = v(6*i-r+7:6*i-r+12)
            end if
         end if
         !
         wd(6*i+1:6*i+6) = 1/this%ds &
            * (  tan_op_inv_mult_s3sdr3(this%ds*w(6*i+1:6*i+6), tmpv) &
               - tan_op_inv_mult_s3sdr3(-this%ds*w(6*i+1:6*i+6), v(6*i-r+1:6*i-r+6)) )

      end function calculate_wd

      pure function myphi(this,q) result(rslt)
         ! input
         class(myproblem),                 intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         ! result
         real(8), dimension(this%sizel)               :: rslt
         ! internal
         real(8), dimension(6*this%n)                 :: w
         integer                                      :: i
         integer                                      :: start
         !

         !ERROR STOP "Function myphi is not implemented and should not be reached"

         ! For debugging purposes, define a custom for-block
#define for(i,f,l) forall(i=f:l)
#define endfor end forall
#ifdef pure
#undef for
#undef endfor
#define for(i,f,l) do i=f,l
#define endfor end do
#endif

         ! calculate w_{k+1/2}
         w = calculate_w(this, q)

         ! Check for Kirchhoff rod
         if (this%kirchhoff == 1) then
            for(i,0,this%n-1)
               rslt(i+1       ) = this%ds*w(6*i+3+1)
               rslt(i+1+this%n) = this%ds*w(6*i+3+2)
            endfor
!            start = 2*this%n1
!         else
!            start = 0
          end if
!
!         ! check for inextensible rod
!         if (this%inextensible == 1) then
!            for(i,1,this%n1)
!               rslt(start + i) = dot_product(Dx(:,i), Dx(:,i)) - this%ds**2
!            endfor
!         end if
!#undef for
!#undef endfor
      end function myphi

!! DEBUG
!   pure function my_num_B(this, q) result(rslt)
!      ! input
!      class(myproblem),   intent(in)            :: this
!      real(8), dimension(:), intent(in)         :: q
!      ! result
!      real(8), dimension(this%sizel,this%sizev) :: rslt
!      ! internal
!      integer                                   :: i
!      real(8), dimension(this%sizev)            :: w
!      real(8), dimension(this%sizev)            :: Phi0
!      !
!      !
!      ! Calculate Phi
!      Phi0 = this%GL(INTEGRATOR)_Phi(q)
!
!      ! set w to zero
!      w = 0.0_8
!      ! loop over the columns of Ct
!      do i=1,this%sizev
!         if (.not. i == 1) w(i-1) = 0.0_8
!         w(i) = max( abs(q(i))*1.0e-8_8,  1.0e-12_8)*1.0e8
!         rslt(:,i) = (this%GL(INTEGRATOR)_Phi(this%GL(INTEGRATOR)_qlpexphDqtilde(q,1.0_8, w)) - Phi0) / w(i)
!         !w    = 0.0_8
!         !w(i) = h
!         !rslt(:,i) = (this%GL(INTEGRATOR)_g(q+this%GL(INTEGRATOR)_tilde(w),v,t) - g0)/h
!      end do
!   end function my_num_B
!! GUBED

      pure function myB(this,q) result(rslt)
         ! input
         class(myproblem), intent(in)  :: this
         real(8),          intent(in)  :: q(:)
         ! result
         real(8)                       :: rslt(this%sizel,this%sizev)
         ! internal variables
         real(8)                       :: w(6*this%n)
         real(8)                       :: T(6,6)
         integer                       :: r
         integer                       :: i

         !ERROR STOP "function myB is not implemented and should not be reached"
!         ! get x and p
!         call xp_from_q(this, x, p, q)

         ! Result is zero in most places
         rslt = 0.0_8
!
!         ! We need Delta x
!         Dx = x(:,2:this%n) - x(:,1:this%n1)
!
!
         ! Calculate w_{k+1/2}
         w = calculate_w(this, q)

         ! Check for Kirchhoff rod
         if (this%kirchhoff == 1) then
            ! We need to take into account, that x1 and xn may be fixed,
            ! therefore the cases i=0 and i=n must be handled a little
            ! different.
            !i = 1
            T = -tan_op_inv_s3sdr3(-this%ds*w(1:6))
            if (this%fixed_x0 == 1) then
               if (this%fixed_p0 == 1) then
                  r = 6
               else
                  rslt(       1, 1:3) =  T(4,1:3)
                  rslt(this%n+1, 1:3) =  T(5,1:3)
                  r = 3
               end if
            else
               if (this%fixed_p0 == 1) then
                  ERROR STOP "free x, fixed p not implemented"
               else
                  rslt(       1, 1:6) =  T(4,:)
                  rslt(this%n+1, 1:6) =  T(5,:)
                  r = 0;
               end if
            end if
!

            do i=1,this%n-1
               T = tan_op_inv_s3sdr3(this%ds*w(6*(i-1)+1:6*(i-1)+6))
               rslt(       i, 6*i+1-r:6*i+6-r) = T(4,:)
               rslt(this%n+i, 6*i+1-r:6*i+6-r) = T(5,:)
               T = -tan_op_inv_s3sdr3(-this%ds*w(6*i+1:6*i+6))
               rslt(       i+1, 6*i+1-r:6*i+6-r) = T(4,:)
               rslt(this%n+i+1, 6*i+1-r:6*i+6-r) = T(5,:)
            end do

            i = this%n
            T = tan_op_inv_s3sdr3(this%ds*w(6*(i-1)+1:6*(i-1)+6))
            if (this%fixed_xn == 1) then
               if (this%fixed_pn == 1) then
                  ! nothing
               else
                  rslt(       i, 6*i+1-r:6*i+3-r) = T(4,1:3)
                  rslt(this%n+i, 6*i+1-r:6*i+3-r) = T(5,1:3)
               end if
            else
               if (this%fixed_p0 == 1) then
                  ERROR STOP "free x, fixed p not implemented"
               else
                  rslt(       i, 6*i+1-r:6*i+6-r) = T(4,:)
                  rslt(this%n+i, 6*i+1-r:6*i+6-r) = T(5,:)
               end if
            end if
!
!            i = this%n1
!            rslt(        i, ind+1:ind+3) = -pvkptoRi(p(:,i), 1)
!            rslt(        i, ind+4:ind+6) =  cross(e1, apply_conj_quat(p(:,i),Dx(:,i)))
!            rslt(this%n1+i, ind+1:ind+3) = -pvkptoRi(p(:,i), 2)
!            rslt(this%n1+i, ind+4:ind+6) =  cross(e2, apply_conj_quat(p(:,i),Dx(:,i)))
!            if (this%fixed_xn == 0) then
!               rslt(        i, ind+7:ind+9) = -rslt(        i, ind+1:ind+3)
!               rslt(this%n1+i, ind+7:ind+9) = -rslt(this%n1+i, ind+1:ind+3)
!            end if
!
!            start = 2*this%n1
!         else
!            start = 0
         end if !(this%kirchhoff == 1)
!
!
!         ! Check for inextensible rod
!         if (this%inextensible == 1) then
!            ! Handle fixed first position correctly
!            if (this%fixed_x1 == 0) then
!               rslt(start+1, 1:3) = -2*Dx(:,1)
!               ind = 6
!            else
!               ind = 3
!            end if
!            rslt(start+1, ind+1:ind+3) =  2*Dx(:,1)
!
!            ! Loop over the contraints except i=1 and i=this%n1
!            do i=2,this%n1-1
!               rslt(start+i, ind+1:ind+3) = -2*Dx(:,i)
!               rslt(start+i, ind+7:ind+9) =  2*Dx(:,i)
!               ind = ind + 6
!            end do
!
!            ! Handle fixed last position correctly
!            rslt(start+this%n1, ind+1:ind+ 3) = -2*Dx(:,this%n1)
!            if (this%fixed_xn == 0) then
!               rslt(start+this%n1, ind+7:ind+9) =  2*Dx(:,this%n1)
!            end if
!         end if !(this%inextensible == 1)
!
      !! DEBUG
      !call print_matrix(rslt,'B')
      !!pause
      !! GUBED

      end function myb

      pure function dtaninvSE3(v,dv) result(rslt)
         ! input
         real(8), intent(in)     :: v(6)
         real(8), intent(in)     :: dv(6)
         ! result
         real(8)                 :: rslt(6,6)
         ! internal
         real(8)                 :: nOm
         real(8)                 :: f, g, h
         real(8), dimension(3,3) :: EOm, EU, EdOm, EdU
         !
         associate(Om  =>  v(1:3),&
                   U   =>  v(4:6),&
                   dOm => dv(1:3),&
                   dU  => dv(4:6))
            nOm = norm2(Om)

            EOm  = skw(Om)
            EU   = skw(U)
            EdOm = skw(dOm)
            EdU  = skw(dU)

            if (nOm < 4.0e-3_8) then
               f = 1.0_8/12 + nOm**2/720 + nOm**4/30240
            else
               f = (2 - nOm/tan(nOm/2))/(2*nOm**2)
            end if
            if (nOm < 2.0e-1_8) then
               g = 1.0_8/360 + nOm**2/7560 + nOm**4/201600 + nOm**6/5987520
            else
               g = (-4 + nOm**2 + 4*cos(nOm) + nOm*sin(nOm))/(4*nOm**4*sin(nOm/2)**2)
            end if
            if (nOm < 2.0e-1_8) then
               h = 1.0_8/3780 + nOm**2/50400 + nOm**4/997920 + 691*nOm**6/16345929600.0_8 + nOm**8/622702080
            else
               h = (64 - 12*nOm/tan(nOm/2) + (4*nOm**2*(3 + nOm/tan(nOm/2)))/(-1 + cos(nOm)))/(8*nOm**6)
            end if

            rslt(1:3,1:3) = EdOm/2 &
              + g * dot_product(dOm,Om) * matmul(EOm,EOm) &
              + f * (matmul(EOm,EdOm) + matmul(EdOm,EOm))
            rslt(4:6,4:6) = rslt(1:3,1:3)
            rslt(1:3,4:6) = 0.0_8
            rslt(4:6,1:3) = EdU/2 &
              + g * dot_product(dOm,Om) * (matmul(EU,EOm) + matmul(EOm,EU)) &
              + f * (matmul(EdU,EOm) + matmul(EU,EdOm) + matmul(EdOm,EU) + matmul(EOm,EdU)) &
              + h * dot_product(dOm,Om) * dot_product(Om,U) * matmul(EOm,EOm) &
              + g * (dot_product(dOm,U) + dot_product(Om,dU)) * matmul(EOm,EOm) &
              + g * dot_product(Om,U) * (matmul(EdOm,EOm) + matmul(EOm,EdOm))
         end associate
      end function dtaninvSE3

      pure function myZ(this,q,v) result(rslt)
         ! input
         class(myproblem), intent(in)  :: this
         real(8),          intent(in)  :: q(:)
         real(8),          intent(in)  :: v(:)
         ! result
         real(8)                       :: rslt(this%sizel)
         ! internal variables
         real(8)                       :: w(6*this%n)
         real(8)                       :: tmpv(6), tmp(6)
         integer                       :: i, r

         !ERROR STOP "function myZ is not implemented and should not be reached"
!         ! get x and p
!         call xp_from_q(this, x, p, q)
!
!         ! get xd and Om
!         call xdom_from_v(this, xd, Om, v)
!
!         ! set result to zero, because we will be adding up on it
!         rslt = 0.0_8
!
!         ! We need Delta x
!         Dx = x(:,2:this%n) - x(:,1:this%n1)
!
         ! TODO TODO TODO

         ! Check for Kirchhoff rod
         if (this%kirchhoff == 1) then
            ! Calculate w_{k+1/2}
            w = calculate_w(this, q)

            if (this%fixed_x0 == 1) then
               if (this%fixed_p0 == 1) then
                  tmpv(1:3) = 0.0_8
                  r = 6
               else
                  tmpv(1:3) = v(1:3)
                  r = 3
               end if
               tmpv(4:6) = 0.0_8
            else
               if (this%fixed_p0 == 1) then
                  ERROR STOP "free x, fixed p not implemented"
               else
                  tmpv = v(1:6)
                  r = 0
               end if
            end if
            tmp =   tan_op_inv_mult_s3sdr3( this%ds*w(1:6),v(6+1-r:6+6-r)) &
                  - tan_op_inv_mult_s3sdr3(-this%ds*w(1:6),tmpv)
            tmp =   matmul(dtaninvSE3(-this%ds*w(1:6),tmp),tmpv) &
                  + matmul(dtaninvSE3( this%ds*w(1:6),tmp),v(6+1-r:6+6-r))
            rslt(       1) = tmp(4)
            rslt(this%n+1) = tmp(5)

            do i=1,this%n-2
               tmp =   tan_op_inv_mult_s3sdr3( this%ds*w(6*i+1:6*i+6),v(6*(i+1)+1-r:6*(i+1)+6-r)) &
                     - tan_op_inv_mult_s3sdr3(-this%ds*w(6*i+1:6*i+6),v(6* i   +1-r:6* i   +6-r))
               tmp =   matmul(dtaninvSE3(-this%ds*w(6*i+1:6*i+6),tmp),v(6* i   +1-r:6* i   +6-r)) &
                     + matmul(dtaninvSE3( this%ds*w(6*i+1:6*i+6),tmp),v(6*(i+1)+1-r:6*(i+1)+6-r))
               rslt(       i+1) = tmp(4)
               rslt(this%n+i+1) = tmp(5)
            end do

            i=this%n-1
            if (this%fixed_xn == 1) then
               if (this%fixed_pn == 1) then
                  tmpv(1:3) = 0.0_8
               else
                  tmpv(1:3) = v(6*(i+1)+1-r:6*(i+1)+3-r)
               end if
               tmpv(4:6) = 0.0_8
            else
               if (this%fixed_pn == 1) then
                  ERROR STOP "free x, fixed p not implemented"
               else
                  tmpv = v(6*(i+1)+1-r:6*(i+1)+6-r)
               end if
            end if
            tmp =   tan_op_inv_mult_s3sdr3( this%ds*w(6*i+1:6*i+6),tmpv) &
                  - tan_op_inv_mult_s3sdr3(-this%ds*w(6*i+1:6*i+6),v(6* i   +1-r:6* i   +6-r))
            tmp =   matmul(dtaninvSE3(-this%ds*w(6*i+1:6*i+6),tmp),v(6* i   +1-r:6* i   +6-r)) &
                  + matmul(dtaninvSE3( this%ds*w(6*i+1:6*i+6),tmp),tmpv)
            rslt(       i+1) = tmp(4)
            rslt(this%n+i+1) = tmp(5)
         end if ! (this%kirchhoff == 1)


!            zerothisn1 = [0, this%n1]
!            ! Loop over the contraints
!            do j=1,2
!               do i=1,this%n1
!                  ! Handle fixed x1 and fixed xn correctly
!                  if (i > 1 .or. this%fixed_x1 == 0) then
!                     ! Derivative wrt. $p_{i+1/2}$ of the derivative wrt. $x_i$
!                     !  is the same as the derivative wrt. $x_i$ of the derivative wrt. $p_{i+1/2}$
!                     rslt(zerothisn1(j)+i) = rslt(zerothisn1(j)+i) -  &
!                        2*dot_product(ej(:,j), cross(apply_conj_quat(p(:,i), xd(:,i)), Om(:,i)))
!                  end if
!                  if (i < this%n1 .or. this%fixed_xn == 0) then
!                     ! Derivative wrt. $p_{i+1/2}$ of the derivative wrt. $x_{i+1}$
!                     !  is the same as the derivative wrt. $x_{i+1}$ of the derivative wrt. $p_{i+1/2}$
!                     rslt(zerothisn1(j)+i) = rslt(zerothisn1(j)+i) +  &
!                        2*dot_product(ej(:,j), cross(apply_conj_quat(p(:,i), xd(:,i+1)), Om(:,i)))
!                  end if
!                  ! Derivative wrt $p_{i+1/2}$ of the derivatives wrt. $p_{i+1/2}$
!                  rslt(zerothisn1(j)+i) = rslt(zerothisn1(j)+i) +  &
!                     dot_product(ej(:,j), cross(cross(apply_conj_quat(p(:,i), Dx(:,i)), Om(:,i)), Om(:,i)))
!                  ! All other derivatives vanish
!               end do
!            end do
!            start = 2*this%n1
!         else
!            start = 0
!         end if !(this%kirchhoff == 1)
!
!
!         ! Check for inextensible rod
!         if (this%inextensible == 1) then
!            ! Loop over the contraints
!            do i=1,this%n1
!               ! Derivative wrt. $x_i$ of the derivative of $x_i$
!               rslt(start+i) = rslt(start+i) + 2*dot_product(xd(:,i),xd(:,i))
!               ! Derivative wrt. $x_{i+1}$ of the derivative of $x_i$
!               rslt(start+i) = rslt(start+i) - 2*dot_product(xd(:,i+1),xd(:,i))
!               ! Derivative wrt. $x_i$ of the derivative of $x_{i+1}$
!               rslt(start+i) = rslt(start+i) - 2*dot_product(xd(:,i),xd(:,i+1))
!               ! Derivative wrt. $x_{i+1}$ of the derivative of $x_{i+1}$
!               rslt(start+i) = rslt(start+i) + 2*dot_product(xd(:,i+1),xd(:,i+1))
!            end do
!         end if !(this%inextensible == 1)
!
      end function myZ

      pure function mymatZ(this,q,v,T) result(rslt)
         ! input
         class(myproblem),                 intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         real(8), dimension(:),            intent(in) :: v
         real(8), dimension(:,:),          intent(in) :: T
         ! result
         real(8), dimension(this%sizel, this%sizev)   :: rslt
         ! internal
         real(8), dimension(3,3)                      :: tildeOm
         real(8), dimension(3,3)                      :: RT
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
         !! the standard 2-norm
         !rslt = sqrt(sum(v**2))
         !!! the standard maximum norm
         !!rslt = maxval(abs(v))
      end function mynorm

      pure function myCt(this, q, v, t) result(rslt)
         ! input
         class(myproblem),      intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8), dimension(:), intent(in)   :: v
         real(8),               intent(in)   :: t
         ! result
         real(8), dimension(size(v),size(v)) :: rslt
         ! internal
         integer                             :: i
         !
         ERROR STOP "function myCt is not implemented and should not be reached"
         !DUMMY FUNCTION, returns identity matrix
         !rslt = 0.0_8
         !forall (i=1:size(v)) rslt(i,i) = 1.0_8
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
         ! internal
         integer                             :: i
         !
         ERROR STOP "function myKt is not implemented and should not be reached"
         !! DUMMY FUNCTION, returns identity matrix
         !rslt = 0.0_8
         !forall (i=1:size(v)) rslt(i,i) = 1.0_8
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
         ! internal
         integer                             :: i
         !
         ERROR STOP "function myKtl is not implemented and should not be reached"
         !! DUMMY FUNCTION, returns identity matrix
         !rslt = 0.0_8
         !forall (i=1:size(v)) rslt(i,i) = 1.0_8
      end function myKtl

      subroutine print_at_custom_s(fid,prob)
         ! input/output
         integer,          intent(in)  :: fid
         class(myproblem), intent(in)  :: prob
         ! internal
         integer                       :: i
!         integer                       :: startind
!         real(8)                       :: x(3,prob%n)
!         real(8)                       :: p(4,prob%n1)
!         real(8)                       :: xd(3,prob%n)
!         real(8)                       :: Om(3,prob%n1)
!         real(8)                       :: xdd(3,prob%n)
!         real(8)                       :: Omd(3,prob%n1)
!         !

         if (prob%fixed_x0 == 1 .or. &
             prob%fixed_p0 == 1 .or. &
             prob%fixed_xn == 1 .or. &
             prob%fixed_pn == 1) then
            ERROR STOP "print at custom s works only for completely free beams"
         end if

         ! Print all q
         ! Loop over the custom s
         do i=1,size(prob%output_s)
            associate (q   => prob%q,              &
                       idx => prob%output_s_id(i), &
                       ds  => prob%output_s_del(i) )
               if (idx == prob%n) then
                  write(fid) q(7*idx+1:7*idx+7)
               else
                  write(fid) lin_interpol_S3R3(ds, q(7*idx+1:7*idx+7), q(7*(idx+1)+1:7*(idx+1)+7))
               end if
            end associate
         end do

         ! Print all v
         ! Loop over the custom s
         do i=1,size(prob%output_s)
            associate (v   => prob%v,              &
                       idx => prob%output_s_id(i), &
                       ds  => prob%output_s_del(i) )
               if (idx == prob%n) then
                  write(fid) v(6*idx+1:6*idx+6)
               else
                  write(fid) lin_interpol_Rn(ds, v(6*idx+1:6*idx+6), v(6*(idx+1)+1:6*(idx+1)+6))
               end if
            end associate
         end do

#if !defined(INT_RATTLie) && !defined(INT_SHAKELie) && !defined(INT_varint4lie)
         ! Print all vd
         ! Loop over the custom s
         do i=1,size(prob%output_s)
            associate (vd  => prob%vd,             &
                       idx => prob%output_s_id(i), &
                       ds  => prob%output_s_del(i) )
               if (idx == prob%n) then
                  write(fid) vd(6*idx+1:6*idx+6)
               else
                  write(fid) lin_interpol_Rn(ds, vd(6*idx+1:6*idx+6), vd(6*(idx+1)+1:6*(idx+1)+6))
               end if
            end associate
         end do
#endif

      ! Print all lambda
      ! Loop over the custom s
      if (prob%kirchhoff == 1) then
         do i=1,size(prob%output_s)
            associate (l   => prob%l,              &
                       idx => prob%output_s_id_l(i), &
                       ds  => prob%output_s_del_l(i) )
               write(fid) lin_interpol_Rn(ds, l(2*idx+1:2*idx+2), l(2*(idx+1)+1:2*(idx+1)+2))
            end associate
         end do
      end if

#if defined(INT_RATTLie)
      ! Print all lambda
      ! Loop over the custom s
      if (prob%kirchhoff == 1) then
         do i=1,size(prob%output_s)
            associate (l   => prob%lm,               &
                       idx => prob%output_s_id_l(i), &
                       ds  => prob%output_s_del_l(i) )
               write(fid) lin_interpol_Rn(ds, l(2*idx+1:2*idx+2), l(2*(idx+1)+1:2*(idx+1)+2))
            end associate
         end do
      end if
#endif

!
!         if (prob%kirchhoff == 1) then
!            ! Loop over the custom p
!            do i=1,size(prob%s_output_p)
!               associate (idp => prob%s_output_p_id(i), &
!                          ds  => prob%ds_output_p(i)    )
!                  write(fid) ds * prob%l(idp) + (1-ds) * prob%l(idp+1)
!               end associate
!            end do
!            ! Loop over the custom p
!            do i=1,size(prob%s_output_p)
!               associate (idp => prob%s_output_p_id(i), &
!                          ds  => prob%ds_output_p(i)    )
!                  write(fid) ds * prob%l(prob%n1+idp) + (1-ds) * prob%l(prob%n1+idp+1)
!               end associate
!            end do
!            startind = prob%n1*2
!         else
!            startind = 0
!         end if
!
!         if (prob%inextensible == 1) then
!            ! Loop over the custom p
!            do i=1,size(prob%s_output_p)
!               associate (idp => prob%s_output_p_id(i), &
!                          ds  => prob%ds_output_p(i)    )
!                  write(fid) ds * prob%l(startind+idp) + (1-ds) * prob%l(startind+idp+1)
!               end associate
!            end do
!         end if
!
!         if (prob%opts%stab2 == 1) then
!            if (prob%kirchhoff == 1) then
!               ! Loop over the custom p
!               do i=1,size(prob%s_output_p)
!                  associate (idp => prob%s_output_p_id(i), &
!                             ds  => prob%ds_output_p(i)    )
!                     write(fid) ds * prob%eta(idp) + (1-ds) * prob%eta(idp+1)
!                  end associate
!               end do
!               ! Loop over the custom p
!               do i=1,size(prob%s_output_p)
!                  associate (idp => prob%s_output_p_id(i), &
!                             ds  => prob%ds_output_p(i)    )
!                     write(fid) ds * prob%eta(prob%n1+idp) + (1-ds) * prob%eta(prob%n1+idp+1)
!                  end associate
!               end do
!               startind = prob%n1*2
!            else
!               startind = 0
!            end if
!
!            if (prob%inextensible == 1) then
!               ! Loop over the custom p
!               do i=1,size(prob%s_output_p)
!                  associate (idp => prob%s_output_p_id(i), &
!                             ds  => prob%ds_output_p(i)    )
!                     write(fid) ds * prob%eta(startind+idp) + (1-ds) * prob%eta(startind+idp+1)
!                  end associate
!               end do
!            end if
!         end if
!
      end subroutine print_at_custom_s

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
               !open (newunit=126, file="out/testprobcrm_misc")
            case (1)
               if (this%output_t_at == 0 .or. &
                   abs(modmod(this%t,this%t_output_at_multiples_of))  < 1.0e-9_8) then ! TODO: Besserer Wert? DEBUG

                  ! normal output:
                  ! write time to file
                  write (this%out_bin_lun) this%t
                  ! write desired values to file
                  if (this%output_s_at == 1) then
                     call print_at_custom_s(this%out_bin_lun, this)
                  else
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
                  end if

                  ! write misc output (DEBUG)
                  write (this%out_misc_lun,*) this%t , this%GL(INTEGRATOR)_stats%newt_steps_curr, this%GL(INTEGRATOR)_stats%newt_steps_max, this%GL(INTEGRATOR)_stats%newt_steps_avg
                  flush (this%out_misc_lun)

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
         ! internal
         integer                             :: kirch_ind
         integer                             :: inext_ind
         integer                             :: v_ind
         integer                             :: ind
         integer                             :: i
         !

         if (this%fixed_x0 == 1 .or. &
             this%fixed_p0 == 1 .or. &
             this%fixed_xn == 1 .or. &
             this%fixed_pn == 1) then
            ERROR STOP "calculate jour works only for completely free beams"
         end if

         ! initialize the index for kirchhoff rod
         kirch_ind = 1

         !! initialize the index for inextensible case
         !if (this%inextensible == 1 .and. this%kirchhoff == 1) then
         !   inext_ind = 2*this%n1 + 1
         !else
         !   inext_ind = 1
         !end if

         !! check for fixed first position and initialize index of jour
         !! and index for the v
         !if (this%fixed_x1 == 0) then
         !   ! write indices for first position
         !   this%opts%jour(1:3) = [1, 2, 3]
         !   ind = 4
         !   v_ind = 4
         !else
            ind = 1
            v_ind = 1
         !end if

         ! loop
         do i=1,this%n1
            ! write indices for angular velocity and velocity
            this%opts%jour(ind:ind+5) = v_ind + [0, 1, 2, 3, 4, 5]
            ind = ind + 6
            v_ind = v_ind + 6

            if (i < this%n1) then
               ! write indices for kirchhoff
               if (this%kirchhoff == 1) then
                  this%opts%jour(ind  ) = this%sizev + kirch_ind
                  this%opts%jour(ind+1) = this%sizev + this%n + kirch_ind
                  ind = ind + 2
#ifndef INT_RATTLie
                  if (this%opts%stab2 == 1) then
                     this%opts%jour(ind  ) = this%sizev + this%sizel + kirch_ind
                     this%opts%jour(ind+1) = this%sizev + this%sizel + this%n + kirch_ind
                     ind = ind + 2
                  end if
#endif
                  kirch_ind = kirch_ind + 1
               end if
            end if

            !! write indices for inextensible
            !if (this%inextensible == 1) then
            !   this%opts%jour(ind) = this%sizev + inext_ind
            !   ind = ind + 1
            !   if (this%opts%stab2 == 1) then
            !      this%opts%jour(ind) = this%sizev + this%sizel + inext_ind
            !      ind = ind + 1
            !   end if
            !   inext_ind = inext_ind + 1
            !end if

            !! write indices for angular velocity
            !this%opts%jour(ind:ind+2) = v_ind + [0, 1, 2]
            !
            !! write indices for next position
            !if (i < this%n1 .or. this%fixed_xn == 0) then !notice lazy evaluation
            !   this%opts%jour(ind+3:ind+5) = v_ind + [3, 4, 5]
            !   v_ind = v_ind + 6
            !   ind = ind + 6
            !end if
         end do

         !! DEBUG
         !print *, 'ind = ', ind
         !call print_vector_int(this%opts%jour, 'jour')
         !stop
         !! GUBED
      end subroutine calculate_jour

      pure function external_force_and_moment(this,t) result(rslt)
         ! input
         class(myproblem), intent(in   )  :: this
         real(8),          intent(in   )  :: t
         !real(8),          intent(in   )  :: x(:,:)
         !real(8),          intent(in   )  :: p(:,:)
         ! result
         real(8)                          :: rslt(6,this%n+1)
         !
         select case (this%externalfm_type)
            ! gravity
            case (1)
               ! No moments
               rslt(1:3,:) = 0.0_8

               ! Only downward forces
               rslt(4:5,:) = 0.0_8
               rslt(6,2:this%n) = this%externalfm_gravity_g * this%m
               ! Gravitational force at the ends is only half as big,
               ! because the mass is only half as big.
               rslt(6,[1, this%n+1]) = this%externalfm_gravity_g * this%m/2

            ! flying spaghetti
            case (2)
               !if (.false.) then ! DEBUG
               ! Torque is applied only at the end
               rslt(1:3,2:this%n+1) = 0.0_8
               if (t < this%externalfm_flying_spaghetti_increasing_time) then
                  rslt(1:3,1) = this%externalfm_flying_spaghetti_moment_factors       &
                           * this%externalfm_flying_spaghetti_maximum_height          &
                           * t / this%externalfm_flying_spaghetti_increasing_time
               elseif (t <   this%externalfm_flying_spaghetti_increasing_time &
                           + this%externalfm_flying_spaghetti_decreasing_time) then
                  rslt(1:3,1) = this%externalfm_flying_spaghetti_moment_factors        &
                           * this%externalfm_flying_spaghetti_maximum_height           &
                           * ((  this%externalfm_flying_spaghetti_increasing_time      &
                               + this%externalfm_flying_spaghetti_decreasing_time) - t)&
                           / this%externalfm_flying_spaghetti_decreasing_time
               else
                  rslt(1:3,:) = 0.0_8
               end if

               ! Force is applied only at the end
               rslt(4:6,2:this%n+1) = 0.0_8
               if (t < this%externalfm_flying_spaghetti_increasing_time) then
                  rslt(4:6,1) = this%externalfm_flying_spaghetti_force_factors      &
                          * this%externalfm_flying_spaghetti_maximum_height         &
                          * t / this%externalfm_flying_spaghetti_increasing_time
               elseif (t <   this%externalfm_flying_spaghetti_increasing_time &
                           + this%externalfm_flying_spaghetti_decreasing_time) then
                  rslt(4:6,1) = this%externalfm_flying_spaghetti_force_factors        &
                          * this%externalfm_flying_spaghetti_maximum_height           &
                          * ((  this%externalfm_flying_spaghetti_increasing_time      &
                              + this%externalfm_flying_spaghetti_decreasing_time) - t)&
                          / this%externalfm_flying_spaghetti_decreasing_time
               else
                  rslt(4:6,:) = 0.0_8
               end if
               !else
               !   rslt(:,2:this%n+1) = 0.0_8
               !   if (t < 5) then
               !      rslt(4:6,1) = this%externalfm_flying_spaghetti_force_factors*200*(256*(-5 + t)**4 * t**4)/390625
               !      rslt(1:3,1) = this%externalfm_flying_spaghetti_moment_factors*200*(256*(-5 + t)**4 * t**4)/390625
               !   else
               !      rslt(:,1) = 0.0_8
               !   end if
               !end if ! DEBUG

            ! roll-up into a full circle
            case (3)
               ! no forces, only moments at the end
               rslt = 0.0_8
               rslt(2,1) = -this%externalfm_roll_up_factor
               rslt(2,this%n+1) = this%externalfm_roll_up_factor
            case default
               ! This case should never be reached
               ERROR STOP 'Unknown externalfm_type'
         end select
      end function external_force_and_moment

!      pure function external_moment(this,t,x,p) result(rslt)
!         ! input
!         class(myproblem), intent(in   )  :: this
!         real(8),          intent(in   )  :: t
!         real(8),          intent(in   )  :: x(:,:)
!         real(8),          intent(in   )  :: p(:,:)
!         ! result
!         real(8)                          :: rslt(3,this%N1)
!         !
!         select case (this%externalfm_type)
!            ! gravity
!            case (1)
!               rslt = 0.0_8
!
!            ! flying spaghetti
!            case (2)
!               ! Torque is applied only at the end
!               rslt(:,2:) = 0.0_8
!               if (t < this%externalfm_flying_spaghetti_increasing_time) then
!                  rslt(:,1) = this%externalfm_flying_spaghetti_moment_factors       &
!                         * this%externalfm_flying_spaghetti_maximum_height          &
!                         * t / this%externalfm_flying_spaghetti_increasing_time
!               elseif (t <   this%externalfm_flying_spaghetti_increasing_time &
!                           + this%externalfm_flying_spaghetti_decreasing_time) then
!                  rslt(:,1) = this%externalfm_flying_spaghetti_moment_factors        &
!                         * this%externalfm_flying_spaghetti_maximum_height           &
!                         * ((  this%externalfm_flying_spaghetti_increasing_time      &
!                             + this%externalfm_flying_spaghetti_decreasing_time) - t)&
!                         / this%externalfm_flying_spaghetti_decreasing_time
!               else
!                  rslt(:,:) = 0.0_8
!               end if
!            case default
!               ! This case should never be reached
!               rslt = 37.7777777777_8
!         end select
!      end function external_moment

      subroutine myinit(this)
         ! input/output
         class(myproblem), intent(inout)        :: this
         ! internal
         real(8), allocatable :: x(:,:)
         real(8), allocatable :: xd(:,:)
         real(8), allocatable :: p(:,:)
         real(8), allocatable :: Om(:,:)
         real(8), allocatable :: Mq(:,:)
         real(8), allocatable :: vdl(:)
         integer, allocatable :: ipiv(:)
         integer              :: info
         integer              :: i
         integer              :: kirch_dummy, inext_dummy
         real(8)              :: testdummy(6,6)

         !! DEBUG
         !do i=1,6
         !   testdummy(:,i) = 0.0_8
         !   testdummy(i,i) = 1.0_8
         !   testdummy(:,i) = tan_tr_inv_mult_s3sdr3(1.25_8*[-0.000330197_8, -1.46033e-8_8, 8.59721e-11_8, 3.99133e-6_8, 0.0180207_8, 0.999983_8],testdummy(:,i))
         !end do
         !call print_m_matrix(testdummy,'tanTf')
         !pause
         !! GUBED

         !! DEBUG
         !call print_vector(hat_tr_mult_s3sdr3([1.0_8,2._8,3._8,4._8,5._8,6._8],[7.0_8,8._8,9._8,10._8,11._8,12._8]),'rslt')
         !pause
         !! GUBED

         !! DEBUG
         !call !print_vector(inv_s3sdr3([1.0_8,2.0_8,3.0_8,4.0_8,5.0_8,6.0_8,7.0_8]/norm2([1.0_8,2.0_8,3.0_8,4.0_8,5.0_8,6.0_8!,7.0_8])),'inv')
         !pause
         !! GUBED

         !! DEBUG
         !call !print_vector(logt_s3sdr3([1.0_8,2.0_8,3.0_8,4.0_8,5.0_8,6.0_8,7.0_8]/norm2([1.0_8,2.0_8,3.0_8,4.0_8,5.0_8,6.0!_8,7.0_8])),'log')
         !pause
         !! GUBED

         !! DEBUG
         !call print_m_vector([-1.2345678901234567e100_8, 1.0_8],'asd')
         !pause
         !! GUBED

         ! Set current time point (t0)
         this%t  = this%opts%t0

         ! Set the sizes of q and v
         this%sizeq = (4+3)*this%n1
         this%sizev = (3+3)*this%n1
         if (this%fixed_x0 == 1) then
            this%sizeq = this%sizeq - 3
            this%sizev = this%sizev - 3
         end if
         if (this%fixed_p0 == 1) then
            this%sizeq = this%sizeq - 4
            this%sizev = this%sizev - 3
         end if
         if (this%fixed_xn == 1) then
            this%sizeq = this%sizeq - 3
            this%sizev = this%sizev - 3
         end if
         if (this%fixed_pn == 1) then
            this%sizeq = this%sizeq - 4
            this%sizev = this%sizev - 3
         end if

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
#if defined(INT_RATTLie) || defined(INT_varint4lie)
         allocate(this%p(this%sizev))
#endif
#if !defined(INT_RATTLie) && !defined(INT_SHAKELie) && !defined(INT_varint4lie)
         allocate(this%vd(this%sizev))
#endif
#ifdef INT_gena
         allocate(this%a(this%sizev))
#endif

         ! Fill in the initial values
         this%q = fill_initial_q(this, this%p0,  this%x0)
         this%v = fill_initial_v(this, this%Om0, this%V0)

         ! In the constrained case
         if (this%opts%constrained == 1) then
            ! Set the length of lambda (and eta)
            ! The Kirchhoff case needs two times the number of segments
            ! as contraints (the cross sections stays orthogonal to
            ! the centerline).
            ! In the inextensible case, we have the
            ! constraint, that no segment changes length.
            if (this%kirchhoff == 1) then
               this%sizel = 2*this%n
            else
               this%sizel = 0
            end if
!            if (this%inextensible == 1) then
!               this%sizel = this%sizel + this%n1
!            end if

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
            ! allocate this%l
            if (allocated(this%l)) then
               deallocate(this%l)
            end if
            allocate(this%l(this%sizel))
            ! Initial values for lambda are calculated by integrator

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

            ! allocate and calculate ordering vector for the iteration matrix
            if (allocated(this%opts%jour)) then
               deallocate(this%opts%jour)
            end if
#ifdef INT_RATTLie
            allocate(this%opts%jour(this%sizev+this%sizel))
#else
            if (this%opts%stab2 == 1) then
               allocate(this%opts%jour(this%sizev+2*this%sizel))
            else
               allocate(this%opts%jour(this%sizev+this%sizel))
            end if
#endif
            call this%calculate_jour()
         end if

         ! Allocate and create variables needed for output at custom
         ! points of the rod
         if (this%output_s_at == 1) then
            if (allocated(this%output_s_id)) then
               deallocate(this%output_s_id)
            end if
            if (allocated(this%output_s_del)) then
               deallocate(this%output_s_del)
            end if

            ! Allocate private variables for outputting at different s
            allocate(this%output_s_id(size(this%output_s)))
            allocate(this%output_s_del(size(this%output_s)))
            ! Fill them. output_s_id is the index of the lower
            ! adjoint q and output_s_del is the normalized distance
            ! between the lower adjoint and the upper adjoint,
            ! where 0<=del<1
            do concurrent (i=1:size(this%output_s))
               this%output_s_id(i) = floor(this%n * this%output_s(i))
               this%output_s_del(i) = this%n * this%output_s(i) - this%output_s_id(i)
            end do

            if (this%opts%constrained == 1) then
               ! The Lagrange multipliers are defined in the middle between the
               ! configurations
               if (allocated(this%output_s_id_l)) then
                  deallocate(this%output_s_id_l)
               end if
               if (allocated(this%output_s_del_l)) then
                  deallocate(this%output_s_del_l)
               end if

               ! Allocate private variables for outputting at different s
               allocate(this%output_s_id_l(size(this%output_s)))
               allocate(this%output_s_del_l(size(this%output_s)))
               ! Fill them. output_s_id is the index of the lower
               ! adjoint q and output_s_del is the normalized distance
               ! between the lower adjoint and the upper adjoint,
               ! where 0<=del<1
               do concurrent (i=1:size(this%output_s))
                   this%output_s_id_l(i) = floor(this%n * this%output_s(i) - 0.5_8)
                   ! in the case of very small s, we are left of the actual first lambda
                   ! so we have to extrapolate
                   if (this%output_s_id_l(i) < 0) this%output_s_id_l(i) = 0
                   ! we could do the same for the right end, but this will only ever happen
                   ! when the requested s is bigger than one
                   if (this%output_s_id_l(i) > this%n1) this%output_s_id_l(i) = this%n1
                   !
                   this%output_s_del_l(i) = this%n * this%output_s(i) - 0.5_8 - this%output_s_id_l(i)
               end do
            end if
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

end module testprobcrm
