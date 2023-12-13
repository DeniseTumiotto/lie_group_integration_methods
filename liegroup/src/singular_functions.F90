module singular_functions

implicit none

contains

   pure function sinx2_x(x) result(rslt)
      ! input
      real(8), intent(in)  :: x
      ! result
      real(8)              :: rslt
      !
      if (x < 1.0e-4_8) then
         rslt = 1.0_8/2 - x**2/48
      else
         rslt = sin(x/2)/x
      end if
   end function sinx2_x

   pure function sinx_x(x) result(rslt)
      ! input
      real(8), intent(in)  :: x
      ! result
      real(8)              :: rslt
      !
      if (x < 1.0e-4_8) then
         rslt = 1.0_8 - x**2/6
      else
         rslt = sin(x)/x
      end if
   end function sinx_x

   pure function cosx_1_x2(x) result(rslt)
      ! input
      real(8), intent(in)  :: x
      ! result
      real(8)              :: rslt
      !
      if (x < 1.0e-2_8) then
         rslt = -1.0_8/2 + x**2/24 - x**4/720
      else
         rslt = (cos(x) - 1.0_8)/x**2
      end if
   end function cosx_1_x2

   pure function x_sinx_x3(x) result(rslt)
      ! input
      real(8), intent(in)  :: x
      ! result
      real(8)              :: rslt
      !
      if (x < 1.0e-4_8) then
         rslt = 1.0_8/6 - x**2/120 + x**4/5040
      else
         rslt = (x - sin(x))/x**3
      end if
   end function x_sinx_x3

   pure function two_2cosx_xsinx_x4(x) result(rslt)
      ! input
      real(8), intent(in)  :: x
      ! result
      real(8)              :: rslt
      !
      if (x < 1.0e-1_8) then
         rslt = 1.0_8/12 - x**2/180 + x**4/6720 - x**6/453600
      else
         rslt = (2.0_8 - 2*cos(x) - x*sin(x))/x**4
      end if
   end function two_2cosx_xsinx_x4

   pure function x_2_cosx_3sinx_x5(x) result(rslt)
      ! input
      real(8), intent(in)  :: x
      ! result
      real(8)              :: rslt
      !
      if (x < 1.0e-1_8) then
         rslt = 1.0_8/60 - x**2/1260 + x**4/60480 - x**6/4989600
      else
         rslt = (x*(2.0_8 + cos(x)) - 3*sin(x))/x**5
      end if
   end function x_2_cosx_3sinx_x5

   pure function two_xcotx2_2x2(x) result(rslt)
      ! input
      real(8), intent(in)  :: x
      ! result
      real(8)              :: rslt
      !
      if (x < 1.0e-2_8) then
         rslt = 1.0_8/12 + x**2/720 + x**4/30240
      else
         rslt = (2.0_8 - x/tan(x/2))/(2*x**2)
      end if
   end function two_xcotx2_2x2

   pure function two_acosx_sqrt_1_x2(x) result(rslt)
      ! input
      real(8), intent(in)  :: x
      ! result
      real(8)              :: rslt
      !
      if (x > 1.0_8 - 5.0e-3_8) then
         rslt = 2.0_8 - 2*(x-1)/3 + 4*(x-1)**2/15 - 4*(x-1)**3/35 + 16*(x-1)**4/315 - 16*(x-1)**5/693
      else
         rslt = 2*acos(x)/sqrt(1 - x**2)
      end if
   end function two_acosx_sqrt_1_x2

   pure function xsinx_4cosx_x2_4_4sinx22_x4(x) result(rslt)
      ! input
      real(8), intent(in)  :: x
      ! result
      real(8)              :: rslt
      !
      if (x < 2.0e-1_8) then
         rslt = 1.0_8/360 + x**2/7560 + x**4/201600 + x**6/5987520 + x**8 * 5.284190138687493e-9_8
      else
         rslt = (x*sin(x) + 4*cos(x) + x**2 - 4.0_8)/(4 * sin(x/2)**2 * x**4)
      end if
   end function xsinx_4cosx_x2_4_4sinx22_x4

   pure function sixtyfour_12xcotx2_4x2_3_xcotx2_cosx_1_8x6(x) result(rslt)
      ! input
      real(8), intent(in)  :: x
      ! resul
      real(8)              :: rslt
      !
      if (x < 4.0e-1_8) then
         rslt = 1.0_8/3780 + x**2/50400 + x**4/997920 + x**6 * 4.227352110949995e-8_8
      else
         rslt = (64.0_8 - 12*x/tan(x/2) + (4*x**2*(3.0_8 + x/tan(x/2)))/(cos(x) - 1.0_8))/(8*x**6)
      end if
   end function sixtyfour_12xcotx2_4x2_3_xcotx2_cosx_1_8x6
end module singular_functions
