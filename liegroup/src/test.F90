program test

   use cross_functions
   use quaternion_functions
   use singular_functions
   use s3_functions
   use s3sdr3_functions

   implicit none

   logical, parameter  :: explicit = .true. ! DEBUG

   integer  :: failures
   integer  :: i
   real(8)  :: p(4), v(3)
   real(8)  :: rsltp(4), rsltv(3)


   failures = 0

   ! Testing for module cross_functions
   write (*, '("########### Testing cross_functions")')
   print *, 'dyadic_product'
   call compare_matrix_matrix(dyadic_product([1.0_8, 2.0_8, 3.0_8],[4.0_8, 5.0_8]), &
      reshape([4.0_8,  8.0_8, 12.0_8, &
               5.0_8, 10.0_8, 15.0_8], [3,2]))
   call resume

   ! Testing for module quaternion_functions
   write (*, '("########### Testing quaternion_functions")')
   do i=1,10
      call random_number(p)
      p = 2*p-1
      p = p/norm2(p)
      call random_number(v)

      rsltv = apply_quat(p,v)
      rsltp = qp(qp(p,[0.0_8,v(1),v(2),v(3)]),conj_quat(p))

      call compare_scalar(rsltp(1),0.0_8)
      call compare_vector(rsltv, rsltp(2:4))
   end do
   call resume

   ! Testing for module singular_functions
   write (*, '("########### Testing singular_functions")')
   print *, 'sinx2_x'
   call compare_scalar(sinx2_x(1.0e-10_8), 0.5_8                )
   call compare_scalar(sinx2_x(1.0e-01_8), 0.4997916927067833_8 )
   call compare_scalar(sinx2_x(1.0e+00_8), 0.479425538604203_8  )
   call compare_scalar(sinx2_x(1.0e+01_8),-0.09589242746631385_8)
   print *, 'sinx_x'
   call compare_scalar(sinx_x(1.0e-10_8), 1.0_8                )
   call compare_scalar(sinx_x(1.0e-01_8), 0.9983341664682815_8 )
   call compare_scalar(sinx_x(1.0e+00_8), 0.8414709848078965_8 )
   call compare_scalar(sinx_x(1.0e+01_8),-0.05440211108893698_8)
   print *, 'cosx_1_x2'
   call compare_scalar(cosx_1_x2(1.0e-10_8),-0.5_8                )
   call compare_scalar(cosx_1_x2(1.0e-01_8),-0.4995834721974234_8 )
   call compare_scalar(cosx_1_x2(1.0e+00_8),-0.4596976941318603_8 )
   call compare_scalar(cosx_1_x2(1.0e+01_8),-0.01839071529076452_8)
   print *, 'x_sinx_x3'
   call compare_scalar(x_sinx_x3(1.0e-10_8),0.1666666666666667_8 )
   call compare_scalar(x_sinx_x3(1.0e-01_8),0.1665833531718477_8 )
   call compare_scalar(x_sinx_x3(1.0e+00_8),0.1585290151921035_8 )
   call compare_scalar(x_sinx_x3(1.0e+01_8),0.01054402111088937_8)
   print *, 'two_2cosx_xsinx_x4'
   call compare_scalar(two_2cosx_xsinx_x4(1.0e-10_8),0.08333333333333333_8  )
   call compare_scalar(two_2cosx_xsinx_x4(1.0e-01_8),0.08327779265652578_8  )
   call compare_scalar(two_2cosx_xsinx_x4(1.0e+00_8),0.07792440345582406_8  )
   call compare_scalar(two_2cosx_xsinx_x4(1.0e+01_8),0.0009118354167046603_8)
   print *, 'x_2_cosx_3sinx_x5'
   call compare_scalar(x_2_cosx_3sinx_x5(1.0e-10_8),0.01666666666666667_8  )
   call compare_scalar(x_2_cosx_3sinx_x5(1.0e-01_8),0.01665873181196891_8  )
   call compare_scalar(x_2_cosx_3sinx_x5(1.0e+00_8),0.0158893514444502_8   )
   call compare_scalar(x_2_cosx_3sinx_x5(1.0e+01_8),0.0001324134804190358_8)
   print *, 'two_xcotx2_2x2'
   call compare_scalar(two_xcotx2_2x2(1.0e-10_8),0.08333333333333333_8)
   call compare_scalar(two_xcotx2_2x2(1.0e-01_8),0.08334722552992746_8)
   call compare_scalar(two_xcotx2_2x2(1.0e+00_8),0.08475613914377404_8)
   call compare_scalar(two_xcotx2_2x2(1.0e+01_8),0.02479064577663728_8)
   print *, 'xsinx_4cosx_x2_4_4sinx22_x4'
   call compare_scalar(xsinx_4cosx_x2_4_4sinx22_x4(1.0e-10_8),0.002777777777777778_8)
   call compare_scalar(xsinx_4cosx_x2_4_4sinx22_x4(1.0e-01_8),0.002779101025299342_8)
   call compare_scalar(xsinx_4cosx_x2_4_4sinx22_x4(1.0e+00_8),0.002915185691236665_8)
   call compare_scalar(xsinx_4cosx_x2_4_4sinx22_x4(1.0e+01_8),0.002370856744723585_8)
   print *, 'two_acosx_sqrt_1_x2'
   call compare_scalar(two_acosx_sqrt_1_x2(1.0_8 - 1.0e-10_8),2.000000000066667_8)
   call compare_scalar(two_acosx_sqrt_1_x2(1.0_8 - 1.0e-01_8),2.069452940470786_8)
   call compare_scalar(two_acosx_sqrt_1_x2(1.0_8 - 1.0e+00_8),3.141592653589793_8)
   print *, 'sixtyfour_12xcotx2_4x2_3_xcotx2_cosx_1_8x6'
   call compare_scalar(sixtyfour_12xcotx2_4x2_3_xcotx2_cosx_1_8x6(1.0e-10_8),0.0002645502645502646_8)
   call compare_scalar(sixtyfour_12xcotx2_4x2_3_xcotx2_cosx_1_8x6(1.0e-01_8),0.0002647487774994004_8)
   call compare_scalar(sixtyfour_12xcotx2_4x2_3_xcotx2_cosx_1_8x6(1.0e+00_8),0.0002854375570870349_8)
   call resume

   ! Testing for module s3_functions
   write (*, '("########### Testing s3_functions")')
   print *, 'lp_s3'
   call compare_vector(lp_s3([1.0_8, 2.0_8, 3.0_8, 4.0_8],[5.0_8, 6.0_8, 7.0_8, 8.0_8]), &
      qp([1.0_8, 2.0_8, 3.0_8, 4.0_8],[5.0_8, 6.0_8, 7.0_8, 8.0_8]))
   print *, 'inv_s3'
   call compare_vector(inv_s3([1.0_8, 2.0_8, 3.0_8, 4.0_8]), conj_quat([1.0_8, 2.0_8, 3.0_8, 4.0_8]))
   print *, 'expt_s3'
   call compare_vector(expt_s3([1.0_8, 2.0_8, 3.0_8]), &
      [-0.2955511274929782_8, 0.2553218600452643_8, 0.5106437200905286_8, 0.7659655801357929_8])
   print *, 'tan_op_s3'
   call compare_matrix(tan_op_s3([1.0_8, 2.0_8, 3.0_8]), &
      [-0.06871266098996704_8, -0.2267181808418462_8, 0.5073830075578865_8, &
        0.5555528457618361_8, 0.1779133377000254_8, 0.3628734929460377_8, &
       -0.0141310101779017_8, 0.6236305018139318_8, 0.5889566688500127_8])
   print *, 'tan_tr_mult_s3'
   call compare_vector(tan_tr_mult_s3([1.0_8, 2.0_8, 3.0_8], [4.0_8, 5.0_8, 6.0_8]), &
      [1.63585649717822_8,5.289019029223697_8,6.595368481458128_8])
   print *, 'tan_op_inv_s3'
   call compare_matrix(tan_op_inv_s3([1.0_8, 2.0_8, 3.0_8]), &
      [-0.4660113732424363_8, 1.725540211268067_8, -0.6616896830978993_8, &
       -1.274459788731933_8, -0.1277010563403356_8, 1.176620633804201_8, &
        1.338310316902101_8, 0.1766206338042014_8, 0.4361494718298322_8])
   print *, 'tan_tr_inv_s3'
   call compare_matrix_matrix(tan_tr_inv_s3([1.0_8, 2.0_8, 3.0_8]), &
         transpose(tan_op_inv_s3([1.0_8, 2.0_8, 3.0_8])))
   print *, 'tan_op_inv_mult_s3'
   call compare_vector(tan_op_inv_mult_s3([1.0_8, 2.0_8, 3.0_8], [-2.0_8, 3.0_8, -4.0_8]), &
      matmul(tan_op_inv_s3([1.0_8, 2.0_8, 3.0_8]), [-2.0_8, 3.0_8, -4.0_8]))
   print *, 'tan_tr_inv_mult_s3'
   call compare_vector(tan_tr_inv_mult_s3([1.0_8, 2.0_8, 3.0_8], [-2.0_8, 3.0_8, -4.0_8]), &
      matmul(tan_tr_inv_s3([1.0_8, 2.0_8, 3.0_8]), [-2.0_8, 3.0_8, -4.0_8]))
   print *, 'logt_s3'
   call compare_vector(logt_s3([0.1825741858350554_8, 0.3651483716701107_8, 0.5477225575051661_8, 0.7302967433402215_8]), &
      [1.03038058532817_8, 1.545570877992255_8, 2.06076117065634_8])
   print *, 'd_tan_tr_inv_s3'
   call compare_matrix(d_tan_tr_inv_s3([1.0_8, 2.0_8, 3.0_8],[4.0_8, 5.0_8, 6.0_8]), &
         [ 3.003741751049189_8, 2.736004487055995_8, -2.423089396648125_8, &
          -4.548113575381870_8, 2.404031720596863_8,  1.138890255997449_8, &
           0.008674478476143939_8, -4.145227806440417_8, 1.809513288932096_8])
   call resume

   ! Testing for module s3sdr3_functions
   write (*, '("########### Testing s3sdr3_functions")')
   print *, 'lp_s3sdr3'
   call compare_vector(lp_s3sdr3(&
         [ 0.1825741858350554_8, 0.3651483716701107_8, 0.5477225575051661_8, 0.7302967433402215_8, &
           1.0_8, 2.0_8, 3.0_8], &
         [-0.6804138174397717_8, 0.5443310539518174_8,-0.408248290463863_8, 0.2721655269759087_8, &
          -3.0_8, 2.0_8, -1.0_8]), &
      [-0.298142396999972_8, 0.298142396999972_8, -0.149071198499986_8, -0.8944271909999159_8, &
        2.533333333333333_8, -1.333333333333333_8, 3.733333333333333_8])
   print *, 'inv_s3sdr3'
   call compare_vector(lp_s3sdr3(&
         [ 0.1825741858350554_8, 0.3651483716701107_8, 0.5477225575051661_8, 0.7302967433402215_8, &
           1.0_8, 2.0_8, 3.0_8], &
         inv_s3sdr3([ 0.1825741858350554_8, 0.3651483716701107_8, 0.5477225575051661_8, 0.7302967433402215_8, &
           1.0_8, 2.0_8, 3.0_8])), &
      [1.0_8, 0.0_8, 0.0_8, 0.0_8,  0.0_8, 0.0_8, 0.0_8])
   print *, 'lie_bracket_s3sdr3'
   call compare_vector(lie_bracket_s3sdr3([ 1.0_8, 2.0_8,  3.0_8, 4.0_8,  5.0_8, 6.0_8], &
                                          [-2.0_8, 3.0_8, -4.0_8, 5.0_8, -6.0_8, 7.0_8]), &
      [-17.0_8, -2.0_8, 7.0_8, -6.0_8, 12.0_8, 6.0_8])
   print *, 'hat_tr_s3sdr3'
   call compare_vector(matmul(transpose(hat_tr_s3sdr3([ 1.0_8, 2.0_8,  3.0_8, 4.0_8,  5.0_8, 6.0_8])), &
         [-2.0_8, 3.0_8, -4.0_8, 5.0_8, -6.0_8, 7.0_8]), &
      [-17.0_8, -2.0_8, 7.0_8, -6.0_8, 12.0_8, 6.0_8])
   print *, 'hat_tr_mult_s3sdr3'
   call compare_vector(matmul(hat_tr_s3sdr3([ 1.0_8, 2.0_8,  3.0_8, 4.0_8,  5.0_8, 6.0_8]), &
         [-2.0_8, 3.0_8, -4.0_8, 5.0_8, -6.0_8, 7.0_8]), &
      hat_tr_mult_s3sdr3([ 1.0_8, 2.0_8,  3.0_8, 4.0_8,  5.0_8, 6.0_8], &
         [-2.0_8, 3.0_8, -4.0_8, 5.0_8, -6.0_8, 7.0_8]))
   print *, 'expt_s3sdr3'
   call compare_vector(expt_s3sdr3([1.0_8, 2.0_8, 3.0_8, 4.0_8, 5.0_8, 6.0_8]), &
      [-0.2955511274929782_8, 0.2553218600452643_8, 0.5106437200905286_8, 0.7659655801357929_8, &
        1.63585649717822_8, 5.289019029223697_8, 6.595368481458128_8])
   print *, 'tan_op_s3sdr3'
   call compare_matrix(tan_op_s3sdr3([1.0_8, 2.0_8, 3.0_8, 4.0_8, 5.0_8, 6.0_8]), &
      [-0.06871266098996704_8, -0.2267181808418462_8, 0.5073830075578865_8, &
       -1.149474050985939_8, 2.577961679532536_8, -0.5474352684191178_8, &
        0.5555528457618361_8, 0.1779133377000254_8, 0.3628734929460377_8, &
       -1.503370590750731_8, -0.9600957311245662_8, 1.044847674592055_8, &
       -0.0141310101779017_8, 0.6236305018139318_8, 0.5889566688500127_8, &
        1.912695902901833_8, 0.205917602233421_8, -0.9732998629422679_8, &
        0.0_8, 0.0_8, 0.0_8, &
       -0.06871266098996704_8, -0.2267181808418462_8, 0.5073830075578865_8, &
        0.0_8, 0.0_8, 0.0_8, &
        0.5555528457618361_8, 0.1779133377000254_8, 0.3628734929460377_8, &
        0.0_8, 0.0_8, 0.0_8, &
       -0.0141310101779017_8, 0.6236305018139318_8, 0.5889566688500127_8])
   print *, 'tan_op_inv_s3sdr3'
   call compare_matrix(tan_op_inv_s3sdr3([1.0_8, 2.0_8, 3.0_8, 4.0_8, 5.0_8, 6.0_8]), &
      [-0.4660113732424363_8, 1.725540211268067_8, -0.6616896830978993_8, &
       -8.981360165037557_8, 4.876201257785771_8, 0.1451467282276067_8, &
       -1.274459788731933_8, -0.1277010563403356_8, 1.176620633804201_8, &
       -1.123798742214229_8, -7.012834070614152_8, 6.275362505748911_8, &
        1.338310316902101_8, 0.1766206338042014_8, 0.4361494718298322_8, &
        5.145146728227607_8, 2.275362505748911_8, -4.183037669111277_8, &
        0.0_8, 0.0_8, 0.0_8, &
       -0.4660113732424363_8, 1.725540211268067_8, -0.6616896830978993_8, &
        0.0_8, 0.0_8, 0.0_8, &
       -1.274459788731933_8, -0.1277010563403356_8, 1.176620633804201_8, &
        0.0_8, 0.0_8, 0.0_8, &
        1.338310316902101_8, 0.1766206338042014_8, 0.4361494718298322_8])
   print *, 'tan_tr_inv_s3sdr3'
   call compare_matrix_matrix(tan_tr_inv_s3sdr3([1.0_8, 2.0_8, 3.0_8, 4.0_8, 5.0_8, 6.0_8]), &
      transpose(tan_op_inv_s3sdr3([1.0_8, 2.0_8, 3.0_8, 4.0_8, 5.0_8, 6.0_8])) )
   print *, 'tan_op_inv_mult_s3sdr3'
   call compare_vector(tan_op_inv_mult_s3sdr3([1.0_8, 2.0_8, 3.0_8, 4.0_8, 5.0_8, 6.0_8], &
                                              [-2.0_8, 3.0_8, -4.0_8, 5.0_8, -6.0_8, 7.0_8]), &
     matmul(tan_op_inv_s3sdr3([1.0_8, 2.0_8, 3.0_8, 4.0_8, 5.0_8, 6.0_8]), &
            [-2.0_8, 3.0_8, -4.0_8, 5.0_8, -6.0_8, 7.0_8]))
   print *, 'tan_tr_inv_mult_s3sdr3'
   call compare_vector(tan_tr_inv_mult_s3sdr3([1.0_8, 2.0_8, 3.0_8, 4.0_8, 5.0_8, 6.0_8], &
                                              [-2.0_8, 3.0_8, -4.0_8, 5.0_8, -6.0_8, 7.0_8]), &
     matmul(tan_tr_inv_s3sdr3([1.0_8, 2.0_8, 3.0_8, 4.0_8, 5.0_8, 6.0_8]), &
            [-2.0_8, 3.0_8, -4.0_8, 5.0_8, -6.0_8, 7.0_8]))
   print *, 'logt_s3sdr3'
   call compare_vector(logt_s3sdr3(&
         [0.1825741858350554_8, 0.3651483716701107_8, 0.5477225575051661_8, 0.7302967433402215_8, &
          1.0_8, 2.0_8, 3.0_8]), &
      [1.03038058532817_8, 1.545570877992255_8, 2.06076117065634_8, &
       1.024006694714424_8, 2.566390627399806_8, 2.563203682092933_8])
   print *, 'd_tan_tr_inv_s3sdr3'
   call compare_matrix(d_tan_tr_inv_s3sdr3([1.0_8, 2.0_8, 3.0_8,  4.0_8,  5.0_8,  6.0_8],  &
                                           [7.0_8, 8.0_8, 9.0_8, 10.0_8, 11.0_8, 12.0_8]), &
      [ 25.38228822722371_8, -1.712241562328065_8, -3.389202354821169_8, &
        6.079202506662695_8,  5.659093883704120_8, -4.592647556140173_8, &
       -29.00216325836857_8,  23.66240557730227_8, -1.508163581144240_8, &
       -10.19326030360948_8,  4.956693049109918_8,  3.769912035600749_8, &
       -31.96904574690218_8, -26.79808527718474_8,  27.79044418959137_8, &
       -1.297355930767367_8, -10.08244215171285_8,  4.300838810455952_8, &
        6.079202506662695_8,  5.659093883704120_8, -4.592647556140173_8, &
                      0.0_8,                0.0_8,                0.0_8, &
       -10.19326030360948_8,  4.956693049109918_8,  3.769912035600749_8, &
                      0.0_8,                0.0_8,                0.0_8, &
       -1.297355930767367_8, -10.08244215171285_8,  4.300838810455952_8, &
                      0.0_8,                0.0_8,                0.0_8 ])
   call resume


contains

   subroutine compare_scalar(in1,in2)
      ! input
      real(8), intent(in)  :: in1, in2
      !
      call print_ok_error_increment_failure( abs(in1-in2) )
   end subroutine compare_scalar

   subroutine compare_vector(in1,in2)
      ! input
      real(8), intent(in)  :: in1(:), in2(:)
      !
      call print_ok_error_increment_failure( norm2(in1-in2) )
   end subroutine compare_vector

   subroutine compare_matrix(in1,in2)
      ! input
      real(8), intent(in)  :: in1(:,:), in2(:)
      !
      call print_ok_error_increment_failure( norm2(in1-reshape(in2,shape(in1))) )
   end subroutine compare_matrix

   subroutine compare_matrix_matrix(in1,in2)
      ! input
      real(8), intent(in)  :: in1(:,:), in2(:,:)
      !
      call print_ok_error_increment_failure( norm2(in1-in2) )
   end subroutine compare_matrix_matrix

   subroutine print_ok_error_increment_failure(res)
      ! input
      real(8), intent(in)  :: res
      !
      if (res < 1.0e-14_8) then
         if (explicit) write(*,'("ok (" ES8.1 ")")') res
      else
         write(*,'("ERROR (" ES8.1 ")")') res
         failures = failures + 1
      end if
   end subroutine print_ok_error_increment_failure

   subroutine resume
      if (failures == 0) then
         write (*, '("No errors so far")')
      else
         write(*,'("Attention! There were " I0 " errors so far!")') failures
      end if
   end subroutine resume

end program test
