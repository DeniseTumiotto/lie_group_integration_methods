module cross_functions

implicit none

contains

   pure function cross(v,w) result(rslt)
      ! input
      real(8), intent(in)  :: v(3)
      real(8), intent(in)  :: w(3)
      ! result
      real(8)              :: rslt(3)
      !
      rslt = [v(2)*w(3) - v(3)*w(2), &
              v(3)*w(1) - v(1)*w(3), &
              v(1)*w(2) - v(2)*w(1)]
   end function cross

   pure function skw(v) result(rslt)
      ! input
      real(8), intent(in)  :: v(3)
      ! result
      real(8)              :: rslt(3,3)
      !
      ! Note that we enter the result COLUMNWISE
      rslt(1:3,1) = [0.0_8,  v(3), -v(2)]
      rslt(1:3,2) = [-v(3), 0.0_8,  v(1)]
      rslt(1:3,3) = [ v(2), -v(1), 0.0_8]
   end function skw

   pure function dyadic_product(x,y) result(rslt)
      ! input
      real(8), intent(in)  :: x(:), y(:)
      ! result
      real(8)              :: rslt(size(x),size(y))
      ! internal
      integer              :: i, j
      !
      forall(i=1:size(x), j=1:size(y)) rslt(i,j) = x(i)*y(j)
   end function dyadic_product

end module cross_functions
