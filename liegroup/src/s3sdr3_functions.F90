module s3sdr3_functions

use cross_functions
use singular_functions
use s3_functions

implicit none

contains

   pure function lp_s3sdr3(q1, q2) result(rslt)
      ! input
      real(8), intent(in)  :: q1(7), q2(7)
      ! result
      real(8)              :: rslt(7)
      !
      rslt(1:4) = qp(q1(1:4), q2(1:4))
      rslt(5:7) = q1(5:7) + apply_quat(q1(1:4), q2(5:7))
   end function lp_s3sdr3

   pure function inv_s3sdr3(q) result(rslt)
      ! input
      real(8), intent(in)  :: q(7)
      ! result
      real(8)              :: rslt(7)
      !
      rslt(1:4) = conj_quat(q(1:4))
      rslt(5:7) = - apply_quat(rslt(1:4), q(5:7))
   end function inv_s3sdr3

   pure function lie_bracket_s3sdr3(v1, v2) result(rslt)
      ! input
      real(8), intent(in)  :: v1(6), v2(6)
      ! result
      real(8)              :: rslt(6)
      !
      associate (Om1 => v1(1:3), Om2 => v2(1:3), &
                 Vv1 => v1(4:6), Vv2 => v2(4:6))
         rslt(1:3) = cross(Om1, Om2)
         rslt(4:6) = cross(Om1, Vv2) - cross(Om2, Vv1)
      end associate
   end function lie_bracket_s3sdr3

   pure function hat_tr_s3sdr3(v) result(rslt)
      ! input
      real(8), intent(in)  :: v(6)
      ! result
      real(8)              :: rslt(6,6)
      !
      rslt(1:3,1:3) = skw(-v(1:3))
      rslt(4:6,4:6) = rslt(1:3,1:3)
      rslt(1:3,4:6) = skw(-v(4:6))
      rslt(4:6,1:3) = 0.0_8
   end function hat_tr_s3sdr3

   pure function hat_tr_mult_s3sdr3(v, w) result(rslt)
      ! input
      real(8), intent(in)  :: v(6), w(6)
      ! result
      real(8)              :: rslt(6)
      !
      rslt(1:3) = cross(w(1:3), v(1:3)) + cross(w(4:6),v(4:6))
      rslt(4:6) = cross(w(4:6), v(1:3))
   end function hat_tr_mult_s3sdr3

   pure function expt_s3sdr3(v) result(rslt)
      ! input
      real(8), intent(in)  :: v(6)
      ! result
      real(8)              :: rslt(7)
      !
      rslt(1:4) = expt_s3(v(1:3))
      rslt(5:7) = tan_tr_mult_s3(v(1:3),v(4:6))
   end function expt_s3sdr3

   pure function tan_op_s3sdr3(v) result(rslt)
      ! input
      real(8), intent(in)  :: v(6)
      ! result
      real(8)              :: rslt(6,6)
      ! internal
      real(8)              :: nOm
      real(8)              :: skwV(3,3), skwOm(3,3)
      !
      rslt(1:3,1:3) = tan_op_s3(v(1:3))
      rslt(4:6,4:6) = rslt(1:3,1:3)
      rslt(1:3,4:6) = 0.0_8
      !
      nOm = norm2(v(1:3))
      !
      skwOm = skw(v(1:3))
      skwV  = skw(v(4:6))
      !
      rslt(4:6,1:3) = cosx_1_x2(nOm) * skwV &
         + x_sinx_x3(nOm) * (matmul(skwV,skwOm) + matmul(skwOm,skwV)) &
         + dot_product(v(1:3),v(4:6)) * ( &
               two_2cosx_xsinx_x4(nOm) * skwOm &
             - x_2_cosx_3sinx_x5(nOm) * matmul(skwOm,skwOm) )
   end function tan_op_s3sdr3

   pure function tan_op_inv_s3sdr3(v) result(rslt)
      ! input
      real(8), intent(in)  :: v(6)
      ! result
      real(8)              :: rslt(6,6)
      ! internal
      real(8)              :: nOm
      real(8)              :: skwOm(3,3), skwV(3,3)
      !
      rslt(1:3,1:3) = tan_op_inv_s3(v(1:3))
      rslt(4:6,4:6) = rslt(1:3,1:3)
      rslt(1:3,4:6) = 0.0_8
      !
      nOm = norm2(v(1:3))
      !
      skwOm = skw(v(1:3))
      skwV  = skw(v(4:6))
      !
      rslt(4:6,1:3) = skwV/2 &
         + two_xcotx2_2x2(nOm) * (matmul(skwV,skwOm) + matmul(skwOm,skwV)) &
         + xsinx_4cosx_x2_4_4sinx22_x4(nOm) * dot_product(v(1:3),v(4:6)) * matmul(skwOm,skwOm)
   end function tan_op_inv_s3sdr3

   pure function tan_tr_inv_s3sdr3(v) result(rslt)
      ! input
      real(8), intent(in)  :: v(6)
      ! result
      real(8)              :: rslt(6,6)
      ! internal
      real(8)              :: nOm
      real(8)              :: skwOm(3,3), skwV(3,3)
      !
      rslt(1:3,1:3) = tan_tr_inv_s3(v(1:3))
      rslt(4:6,4:6) = rslt(1:3,1:3)
      rslt(4:6,1:3) = 0.0_8
      !
      nOm = norm2(v(1:3))
      !
      skwOm = skw(v(1:3))
      skwV  = skw(v(4:6))
      !
      rslt(1:3,4:6) = -skwV/2 &
         + two_xcotx2_2x2(nOm) * (matmul(skwV,skwOm) + matmul(skwOm,skwV)) &
         + xsinx_4cosx_x2_4_4sinx22_x4(nOm) * dot_product(v(1:3),v(4:6)) * matmul(skwOm,skwOm)
   end function tan_tr_inv_s3sdr3

   pure function tan_op_inv_mult_s3sdr3(v, w) result(rslt)
      ! input
      real(8), intent(in)  :: v(6), w(6)
      ! result
      real(8)              :: rslt(6)
      ! internal
      real(8)              :: nOm
      real(8)              :: skwOm(3,3), skwV(3,3)
      !
      rslt(1:3) = tan_op_inv_mult_s3(v(1:3),w(1:3))
      rslt(4:6) = tan_op_inv_mult_s3(v(1:3),w(4:6))
      !
      nOm = norm2(v(1:3))
      !
      skwOm = skw(v(1:3))
      skwV  = skw(v(4:6))
      !
      rslt(4:6) = rslt(4:6) + cross(v(4:6), w(1:3))/2 &
         + two_xcotx2_2x2(nOm) * (cross(v(4:6),cross(v(1:3),w(1:3))) + cross(v(1:3),cross(v(4:6),w(1:3)))) &
         + xsinx_4cosx_x2_4_4sinx22_x4(nOm) * dot_product(v(1:3),v(4:6)) * cross(v(1:3),cross(v(1:3),w(1:3)))
   end function tan_op_inv_mult_s3sdr3

   pure function tan_tr_inv_mult_s3sdr3(v, w) result(rslt)
      ! input
      real(8), intent(in)  :: v(6), w(6)
      ! result
      real(8)              :: rslt(6)
      ! internal
      real(8)              :: nOm
      real(8)              :: skwOm(3,3), skwV(3,3)
      !
      rslt(1:3) = tan_tr_inv_mult_s3(v(1:3),w(1:3))
      rslt(4:6) = tan_tr_inv_mult_s3(v(1:3),w(4:6))
      !
      nOm = norm2(v(1:3))
      !
      skwOm = skw(v(1:3))
      skwV  = skw(v(4:6))
      !
      rslt(1:3) = rslt(1:3) - cross(v(4:6), w(4:6))/2 &
         + two_xcotx2_2x2(nOm) * (cross(v(4:6),cross(v(1:3),w(4:6))) + cross(v(1:3),cross(v(4:6),w(4:6)))) &
         + xsinx_4cosx_x2_4_4sinx22_x4(nOm) * dot_product(v(1:3),v(4:6)) * cross(v(1:3),cross(v(1:3),w(4:6)))
   end function tan_tr_inv_mult_s3sdr3

   pure function logt_s3sdr3(q) result(rslt)
      ! input
      real(8), intent(in)  :: q(7)
      ! result
      real(8)              :: rslt(6)
      !
      rslt(1:3) = logt_s3(q(1:4))
      rslt(4:6) = tan_op_inv_mult_s3(-rslt(1:3), q(5:7))
   end function logt_s3sdr3

   pure function d_tan_tr_inv_s3sdr3(lower_v, lower_w) result(rslt)
      ! input
      real(8), intent(in)  :: lower_v(6), lower_w(6)
      ! result
      real(8)              :: rslt(6,6)
      ! internal
      real(8)              :: nOm
      !
      associate (Om => lower_v(1:3), &
                 V  => lower_v(4:6), &
                 A  => lower_w(1:3), &
                 W  => lower_w(4:6) )
         nOm = norm2(Om)
         rslt(1:3,1:3) = d_tan_tr_inv_s3(Om, A)
         rslt(1:3,4:6) = d_tan_tr_inv_s3(Om, W)
         rslt(4:6,1:3) = rslt(1:3,4:6)
         rslt(4:6,4:6) = 0.0_8
         !
         rslt(1:3,1:3) = rslt(1:3,1:3)                                              &
            + two_xcotx2_2x2(nOm) * (skw(cross(W,V)) - matmul(skw(V), skw(W)))      &
            + xsinx_4cosx_x2_4_4sinx22_x4(nOm) * (                                  &
                - dot_product(V,Om)*( matmul(skw(Om),skw(W))                        &
                                      + skw(cross(Om,W)) )                          &
                + dyadic_product(cross(Om,cross(Om,W)),V)                           &
                + dyadic_product(cross(V,cross(Om,W)) + cross(Om,cross(V,W)), Om) ) &
            + sixtyfour_12xcotx2_4x2_3_xcotx2_cosx_1_8x6(nOm) * dot_product(Om,V)   &
                * dyadic_product(cross(Om,cross(Om,W)), Om)
      end associate
   end function d_tan_tr_inv_s3sdr3

end module s3sdr3_functions
