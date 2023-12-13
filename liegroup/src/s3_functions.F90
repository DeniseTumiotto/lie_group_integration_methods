module s3_functions

use cross_functions
use quaternion_functions
use singular_functions

implicit none

contains

   pure function lp_s3(p1, p2) result(rslt)
      ! input
      real(8), intent(in)  :: p1(4), p2(4)
      ! result
      real(8)              :: rslt(4)
      !
      rslt = qp(p1,p2)
   end function lp_s3

   pure function inv_s3(p) result(rslt)
      ! input
      real(8), intent(in)  :: p(4)
      ! result
      real(8)              :: rslt(4)
      !
      rslt = conj_quat(p)
   end function inv_s3

   pure function lie_bracket_s3(v1, v2) result(rslt)
      ! input
      real(8), intent(in)  :: v1(3), v2(3)
      ! result
      real(8)              :: rslt(3)
      !
      rslt = cross(v1, v2)
   end function lie_bracket_s3

   pure function expt_s3(v) result(rslt)
      ! input
      real(8), intent(in)  :: v(3)
      ! result
      real(8)              :: rslt(0:3)
      ! internal
      real(8)              :: nv
      !
      nv = norm2(v)
      rslt(0) = cos(nv/2)
      rslt(1:3) = v * sinx2_x(nv)
   end function expt_s3

   pure function tan_op_s3(v) result(rslt)
      ! input
      real(8), intent(in)  :: v(3)
      ! result
      real(8)              :: rslt(3,3)
      ! internal
      real(8)              :: nv
      integer              :: i
      !
      nv = norm2(v)
      !
      rslt = skw(v)
      !
      rslt = cosx_1_x2(nv)*rslt + x_sinx_x3(nv)*matmul(rslt,rslt)
      !
      forall (i=1:3)
         rslt(i,i) = rslt(i,i) + 1
      end forall
   end function tan_op_s3

   pure function tan_tr_mult_s3(v,w) result(rslt)
      ! input
      real(8), intent(in)  :: v(3)
      real(8), intent(in)  :: w(3)
      ! result
      real(8)              :: rslt(3)
      ! internal
      real(8)              :: nv
      !
      nv = norm2(v)
      !
      rslt = cross(v,w)
      !
      rslt = w - cosx_1_x2(nv)*rslt + x_sinx_x3(nv)*cross(v,rslt)
   end function tan_tr_mult_s3

   pure function tan_op_inv_s3(v) result(rslt)
      ! input
      real(8), intent(in)  :: v(3)
      ! result
      real(8)              :: rslt(3,3)
      ! internal
      real(8)              :: nv
      integer              :: i
      !
      nv = norm2(v)
      !
      rslt = skw(v)
      !
      rslt = rslt/2 + two_xcotx2_2x2(nv)*matmul(rslt,rslt)
      !
      forall (i=1:3)
         rslt(i,i) = rslt(i,i) + 1
      end forall
   end function tan_op_inv_s3

   pure function tan_tr_inv_s3(v) result(rslt)
      ! input
      real(8), intent(in)  :: v(3)
      ! result
      real(8)              :: rslt(3,3)
      ! internal
      real(8)              :: nv
      integer              :: i
      !
      nv = norm2(v)
      !
      rslt = skw(v)
      !
      rslt = -rslt/2 + two_xcotx2_2x2(nv)*matmul(rslt,rslt)
      !
      forall (i=1:3)
         rslt(i,i) = rslt(i,i) + 1
      end forall
   end function tan_tr_inv_s3

   pure function tan_op_inv_mult_s3(v, w) result(rslt)
      ! input
      real(8), intent(in)  :: v(3), w(3)
      ! result
      real(8)              :: rslt(3)
      ! internal
      real(8)              :: nv
      !
      nv = norm2(v)
      !
      rslt = w + cross(v,w)/2 + two_xcotx2_2x2(nv)*cross(v,cross(v,w))
   end function tan_op_inv_mult_s3

   pure function tan_tr_inv_mult_s3(v, w) result(rslt)
      ! input
      real(8), intent(in)  :: v(3), w(3)
      ! result
      real(8)              :: rslt(3)
      ! internal
      real(8)              :: nv
      !
      nv = norm2(v)
      !
      rslt = w - cross(v,w)/2 + two_xcotx2_2x2(nv)*cross(v,cross(v,w))
   end function tan_tr_inv_mult_s3

   pure function logt_s3(p) result(rslt)
      ! input
      real(8), intent(in)  :: p(0:3)
      ! result
      real(8)              :: rslt(3)
      !
      rslt = two_acosx_sqrt_1_x2(p(0))*p(1:3)
   end function logt_s3

   pure function d_tan_tr_inv_s3(Om, A) result(rslt)
      ! input
      real(8), intent(in)  :: Om(3), A(3)
      ! result
      real(8)              :: rslt(3,3)
      ! internal
      real(8)              :: nOm
      !
      nOm = norm2(Om)
      !
      rslt = skw(a)/2                                                               &
        - xsinx_4cosx_x2_4_4sinx22_x4(nOm)*dyadic_product(cross(Om,cross(A,Om)),Om) &
        - two_xcotx2_2x2(nOm)*(matmul(skw(Om),skw(A)) + skw(cross(Om,A)))
   end function d_tan_tr_inv_s3

end module s3_functions
