module lie_group_functions
   use s3sdr3_functions

contains
   ! calculate the imaginary part of the multiplication of the conjugated p with q
   pure function qikpq(p,q) result(rslt)
      ! input
      real(8), intent(in)  :: p(0:3)
      real(8), intent(in)  :: q(0:3)
      ! result
      real(8)              :: rslt(3)
      !
      rslt = [- p(1)*q(0) + p(0)*q(1) + p(3)*q(2) - p(2)*q(3), &
              - p(2)*q(0) - p(3)*q(1) + p(0)*q(2) + p(1)*q(3), &
              - p(3)*q(0) + p(2)*q(1) - p(1)*q(2) + p(0)*q(3)]
   end function qikpq

   ! calculate skw but save the resulting matrix as a 9-vector
   pure function qEv(v) result(rslt)
      ! input
      real(8), intent(in)  :: v(3)
      ! result
      real(8)              :: rslt(9)
      !
      ! Note that we enter the result COLUMNWISE
      rslt(1:3) = [0.0_8,  v(3), -v(2)]
      rslt(4:6) = [-v(3), 0.0_8,  v(1)]
      rslt(7:9) = [ v(2), -v(1), 0.0_8]
   end function qEv

   ! calculate rotation matrix associated to a quaternion
   pure function qR(p) result(rslt)
      ! input
      real(8), intent(in)  :: p(0:3)
      ! result
      real(8)              :: rslt(3,3)
      !
      ! Note that we enter the result COLUMNWISE
      rslt(1:3,1) = [p(0)*p(0) + p(1)*p(1) - p(2)*p(2) - p(3)*p(3), 2*(p(1)*p(2) + p(0)*p(3)), 2*(p(1)*p(3) - p(0)*p(2))]
      rslt(1:3,2) = [2*(p(1)*p(2) - p(0)*p(3)), p(0)*p(0) - p(1)*p(1) + p(2)*p(2) - p(3)*p(3), 2*(p(0)*p(1) + p(2)*p(3))]
      rslt(1:3,3) = [2*(p(0)*p(2) + p(1)*p(3)), 2*(p(2)*p(3) - p(0)*p(1)), p(0)*p(0) - p(1)*p(1) - p(2)*p(2) + p(3)*p(3)]
   end function qR

   ! see above, but save the matrix as a 9-vector
   pure function qRv(p) result(rslt)
      ! input
      real(8), intent(in)  :: p(0:3)
      ! result
      real(8)              :: rslt(9)
      !
      ! Note that we enter the result COLUMNWISE
      rslt(1:3) = [p(0)*p(0) + p(1)*p(1) - p(2)*p(2) - p(3)*p(3), 2*(p(1)*p(2) + p(0)*p(3)), 2*(p(1)*p(3) - p(0)*p(2))]
      rslt(4:6) = [2*(p(1)*p(2) - p(0)*p(3)), p(0)*p(0) - p(1)*p(1) + p(2)*p(2) - p(3)*p(3), 2*(p(0)*p(1) + p(2)*p(3))]
      rslt(7:9) = [2*(p(0)*p(2) + p(1)*p(3)), 2*(p(2)*p(3) - p(0)*p(1)), p(0)*p(0) - p(1)*p(1) - p(2)*p(2) + p(3)*p(3)]
   end function qRv

   ! dual quaternion product
   pure function dqp(a,b) result(rslt)
      ! input
      real(8), intent(in)  :: a(8)
      real(8), intent(in)  :: b(8)
      ! result
      real(8)              :: rslt(8)
      !
      rslt(1:4) = qp(a(1:4), b(1:4))
      rslt(5:8) = qp(a(1:4), b(5:8)) + qp(a(5:8), b(1:4))
   end function dqp

   ! exponential function  for SO(3)
   pure function so3expt(w) result(rslt)
      ! input
      real(8), intent(in)  :: w(3)
      ! result
      real(8)              :: rslt(3,3)
      ! internal
      real(8)              :: nw
      real(8)              :: f1, f2
      !
      nw = norm2(w)
      !
      rslt(1:3,1) = [1.0_8, 0.0_8, 0.0_8]
      rslt(1:3,2) = [0.0_8, 1.0_8, 0.0_8]
      rslt(1:3,3) = [0.0_8, 0.0_8, 1.0_8]
      !
      rslt = rslt + sinx_x(nw) * skw(w) - cosx_1_x2(nw) * matmul(skw(w),skw(w))
   end function so3expt

   ! see above, but save 3x3-matrix in 9-vector format
   pure function so3exptv(w) result(rslt)
      ! input
      real(8), intent(in)  :: w(3)
      ! result
      real(8)              :: rslt(9)
      ! internal
      real(8)              :: nw
      real(8)              :: f1, f2
      !
      nw = norm2(w)
      !
      rslt(1:3) = [1.0_8, 0.0_8, 0.0_8]
      rslt(4:6) = [0.0_8, 1.0_8, 0.0_8]
      rslt(7:9) = [0.0_8, 0.0_8, 1.0_8]
      !
      rslt = rslt + sinx_x(nw) * qEv(w) - cosx_1_x2(nw) * matmatv(qEv(w),qEv(w))
   end function so3exptv

   ! multiply two 3x3 matrices in the 9-vector format
   pure function matmatv(A, B) result(rslt)
      ! input
      real(8), intent(in)  :: A(9)
      real(8), intent(in)  :: B(9)
      ! result
      real(8)              :: rslt(9)
      ! internal
      integer              :: i, j
      !
      forall (j=1:3)
         forall (i=1:3)
            rslt(i + 3*(j-1)) = dot_product(A(i:i+6:3), B(1+3*(j-1):3+3*(j-1)))
         end forall
      end forall
      !print *, sum((matmul(reshape(A,[3,3]),reshape(B,[3,3]))-reshape(rslt,[3,3]))**2)!DEBUG
   end function matmatv

   ! see above, but result is actually a 3x3-matrix
   pure function matvmat2mat(A, B) result(rslt)
      ! input
      real(8), intent(in)  :: A(9)
      real(8), intent(in)  :: B(3,3)
      ! result
      real(8)              :: rslt(3,3)
      ! internal
      integer              :: i, j
      !
      forall (j=1:3)
         forall (i=1:3)
            rslt(i,j) = dot_product(A(i:i+6:3), B(:,j))
         end forall
      end forall
      !print *, sum((matmul(reshape(A,[3,3]),B)-rslt)**2)!DEBUG
   end function matvmat2mat

   ! multiply 3x3-matrix in 9-vector format with 3-vector
   pure function matvecv(A, B) result(rslt)
      ! input
      real(8), intent(in)  :: A(9)
      real(8), intent(in)  :: b(3)
      ! result
      real(8)              :: rslt(3)
      ! internal
      integer              :: i
      !
      forall (i=1:3)
         rslt(i) = dot_product(A(i:i+6:3), b)
      end forall
      !print *, sum((matmul(reshape(A,[3,3]),b)-rslt)**2)!DEBUG
   end function matvecv

   ! see above, but matrix is transposed
   pure function matTvecv(A, B) result(rslt)
      ! input
      real(8), intent(in)  :: A(9)
      real(8), intent(in)  :: b(3)
      ! result
      real(8)              :: rslt(3)
      ! internal
      integer              :: i
      !
      forall (i=1:3)
         rslt(i) = dot_product(A(1+3*(i-1):3+3*(i-1)), b)
      end forall
      !print *, sum((matmul(transpose(reshape(A,[3,3])),b)-rslt)**2)!DEBUG
   end function matTvecv

   ! multiply two elements of SE(3) in 12-vector format
   pure function se3pv(a, b) result(rslt)
      ! input
      real(8), intent(in)  :: a(12)
      real(8), intent(in)  :: b(12)
      ! result
      real(8)              :: rslt(12)
      !
      rslt( 1: 9) = matmatv(a(1:9), b(1:9))
      rslt(10:12) = a(10:12) + matvecv(a(1:9), b(10:12))
   end function se3pv

   ! calculate diagonal block part of tangent operator of SE(3)
   pure function se3A(w,U) result(rslt)
      ! input
      real(8), intent(in)  :: w(3)
      real(8), intent(in)  :: U(3)
      ! result
      real(8)              :: rslt(3)
      ! internal
      real(8)              :: nw
      real(8)              :: f1, f2
      !
      nw = norm2(w)
      !
      rslt = U - cosx_1_x2(nw) * cross(w,U) + x_sinx_x3(nw) * cross(w,cross(w,U))
   end function se3A

   ! calculate off-diagonal block part of tangent operator of SE(3)
   pure function se3C(w,U) result(rslt)
      ! input
      real(8), intent(in)  :: w(3)
      real(8), intent(in)  :: U(3)
      ! result
      real(8)              :: rslt(3,3)
      ! internal
      real(8)              :: nw, wU
      !
      nw = norm2(w)
      !
      wU = dot_product(w,u)
      !
      rslt =  cosx_1_x2(nw) * skw(U) &
               + x_sinx_x3(nw) * (matmul(skw(U),skw(w)) + matmul(skw(w),skw(U))) &
               + wU * two_2cosx_xsinx_x4(nw) * skw(w) &
               + wU * x_2_cosx_3sinx_x5(nw) * matmul(skw(w),skw(w))
   end function se3C

end module lie_group_functions
