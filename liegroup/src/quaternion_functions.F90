module quaternion_functions

implicit none

contains

   pure function qp(p,q) result(rslt)
      ! input
      real(8), intent(in)  :: p(0:3)
      real(8), intent(in)  :: q(0:3)
      ! result
      real(8)              :: rslt(4)
      !
      rslt = [p(0)*q(0) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3), &
              p(1)*q(0) + p(0)*q(1) - p(3)*q(2) + p(2)*q(3), &
              p(2)*q(0) + p(3)*q(1) + p(0)*q(2) - p(1)*q(3), &
              p(3)*q(0) - p(2)*q(1) + p(1)*q(2) + p(0)*q(3)]
   end function qp

   pure function apply_quat(p,v) result(rslt)
      ! input
      real(8), intent(in)  :: p(0:3)
      real(8), intent(in)  :: v(1:3)
      ! result
      real(8)              :: rslt(3)
      ! internal
      real(8)              :: r(0:3)
      !
      ! first, calculate p*[0\\v], keeping in mind that v(0)=0
      r =  [- p(1)*v(1) - p(2)*v(2) - p(3)*v(3), &
              p(0)*v(1) - p(3)*v(2) + p(2)*v(3), &
              p(3)*v(1) + p(0)*v(2) - p(1)*v(3), &
            - p(2)*v(1) + p(1)*v(2) + p(0)*v(3)]
      ! now calculate Im(r*conj(p)), keeping in mind that p*v*conj(p) will
      ! always yield a vector and directly applying the conj
      rslt = [r(1)*p(0) - r(0)*p(1) + r(3)*p(2) - r(2)*p(3), &
              r(2)*p(0) - r(3)*p(1) - r(0)*p(2) + r(1)*p(3), &
              r(3)*p(0) + r(2)*p(1) - r(1)*p(2) - r(0)*p(3)]
   end function apply_quat

   pure function apply_conj_quat(p,v) result(rslt)
      ! input
      real(8), intent(in)  :: p(0:3)
      real(8), intent(in)  :: v(1:3)
      ! result
      real(8)              :: rslt(3)
      ! internal
      real(8)              :: r(0:3)
      !
      ! first, calculate conj(p)*[0\\v], keeping in mind that v(0)=0 and
      ! directly applying the conj
      r =  [  p(1)*v(1) + p(2)*v(2) + p(3)*v(3), &
              p(0)*v(1) + p(3)*v(2) - p(2)*v(3), &
            - p(3)*v(1) + p(0)*v(2) + p(1)*v(3), &
              p(2)*v(1) - p(1)*v(2) + p(0)*v(3)]
      ! now calculate Im(r*p), keeping in mind that conj(p)*v*p will
      ! always yield a vector
      rslt = [r(1)*p(0) + r(0)*p(1) - r(3)*p(2) + r(2)*p(3), &
              r(2)*p(0) + r(3)*p(1) + r(0)*p(2) - r(1)*p(3), &
              r(3)*p(0) - r(2)*p(1) + r(1)*p(2) + r(0)*p(3)]
   end function apply_conj_quat

   pure function conj_quat(p) result(rslt)
      ! input
      real(8), intent(in)  :: p(0:3)
      ! result
      real(8)              :: rslt(0:3)
      !
      rslt(0) = p(0)
      rslt(1:3) = - p(1:3)
   end function conj_quat

end module quaternion_functions
