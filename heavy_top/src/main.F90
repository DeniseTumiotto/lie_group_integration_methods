#define GL(x) x

program main
   !!!!!! Use modules
   ! Use modul to be able to identify the output unit
   use iso_fortran_env
   ! Use the module with the actual crm problem
   use heavy_top
   ! For reading the config file
   use aotus_module, only: flu_State, open_config_file, close_config, aot_get_val, &
                           aot_top_get_val, &
                           aoterr_Fatal, aoterr_WrongType, aoterr_NonExistent
   use aot_table_module, only:   aot_table_open, aot_table_close, aot_table_length
   use aot_fun_module, only:  aot_fun_type, aot_fun_open, aot_fun_put, &
                              aot_fun_do, aot_fun_close
   ! For reading lines of variable length
   use get_line_of_variable_length

   !!!!!! No implicit variables
   implicit none

   !!!!!! Variables
   ! crm problem object
   type(heavy_top_t)   :: prob
   ! For looping
   integer                 :: i, j
   ! Aotus handles and error variables, etc.
   character(len=256)      :: conf_fname
   type(flu_State)         :: conf
   integer                 :: iError
   integer                 :: iv3Error(3)
   integer                 :: iv4Error(4)
   character(len=256)      :: cError
   integer, parameter      :: max_length = 256
   integer                 :: thandle
   !
   integer                       :: out_lua_lun
   integer                       :: conf_lun
   character(len=:), allocatable :: conf_line
   character(len=128)            :: iomsg_str
   integer                       :: iostat_number


   !!!!!! Get name of the config file
   call get_command_argument(1, conf_fname)
   if (len_trim(conf_fname) == 0) then
      print *, 'FATAL Error: Please provide a config file.'
      ERROR STOP
   end if

   !!!!!! Read config file
   ! Open file
   call open_config_file(L = conf, filename = conf_fname, &
                         ErrCode = iError, ErrString = cError)
   if (iError /= 0) then
      print *, 'FATAL Error ', iError, ' when opening the Lua config file:', cError
      ERROR STOP
   end if

#ifdef INT_gena
   ! Algorithmic parameters
   call aot_get_val(L = conf, key = 'alpha_m', val = prob%GL(INTEGRATOR)_alpha_m, ErrCode = iError)
   call error_check(conf, iError, 'alpha_m')
   print *, 'alpha_m = ', prob%GL(INTEGRATOR)_alpha_m

   call aot_get_val(L = conf, key = 'alpha_f', val = prob%GL(INTEGRATOR)_alpha_f, ErrCode = iError)
   call error_check(conf, iError, 'alpha_f')
   print *, 'alpha_f = ', prob%GL(INTEGRATOR)_alpha_f

   call aot_get_val(L = conf, key = 'beta', val = prob%GL(INTEGRATOR)_beta, ErrCode = iError)
   call error_check(conf, iError, 'beta')
   print *, 'beta = ', prob%GL(INTEGRATOR)_beta

   call aot_get_val(L = conf, key = 'gamma', val = prob%GL(INTEGRATOR)_gamma, ErrCode = iError)
   call error_check(conf, iError, 'gamma')
   print *, 'gamma = ', prob%GL(INTEGRATOR)_gamma
#endif

#ifdef INT_BLieDF
   ! Algorithmic parameters
   call aot_get_val(L = conf, key = 'k_bdf', val = prob%k_bdf, ErrCode = iError)
   call error_check(conf, iError, 'alpha_m')
   print *, 'k_bdf = ', prob%k_bdf
#endif

   ! Integrator options
   call aot_get_val(L = conf, key = 'const_mass_matrix', val = prob%opts%const_mass_matrix, ErrCode = iError)
   call error_check(conf, iError, 'const_mass_matrix')
   print *, 'const_mass_matrix = ', prob%opts%const_mass_matrix

   call aot_get_val(L = conf, key = 'diag_mass_matrix', val = prob%opts%diag_mass_matrix, ErrCode = iError)
   call error_check(conf, iError, 'diag_mass_matrix')
   print *, 'diag_mass_matrix = ', prob%opts%diag_mass_matrix

   call aot_get_val(L = conf, key = 'banded_iteration_matrix', val = prob%opts%banded_iteration_matrix, ErrCode = iError)
   call error_check(conf, iError, 'banded_iteration_matrix')
   print *, 'banded_iteration_matrix = ', prob%opts%banded_iteration_matrix

   if (prob%opts%banded_iteration_matrix == 1) then
      call aot_get_val(L = conf, key = 'nr_subdiag', val = prob%opts%nr_subdiag, ErrCode = iError)
      call error_check(conf, iError, 'nr_subdiag')
      print *, 'nr_subdiag = ', prob%opts%nr_subdiag

      call aot_get_val(L = conf, key = 'nr_superdiag', val = prob%opts%nr_superdiag, ErrCode = iError)
      call error_check(conf, iError, 'nr_superdiag')
      print *, 'nr_superdiag = ', prob%opts%nr_superdiag
   else
      print *, 'nr_superdiag = n/a'
      print *, 'nr_subdiag = n/a'
   end if

   call aot_get_val(L = conf, key = 'recalc_iteration_matrix', val = prob%opts%recalc_iteration_matrix, ErrCode = iError)
   call error_check(conf, iError, 'recalc_iteration_matrix')
   print *, 'recalc_iteration_matrix = ', prob%opts%recalc_iteration_matrix

#ifdef INT_gena
   call aot_get_val(L = conf, key = 'perturb', val = prob%opts%pertube, ErrCode = iError)
   call error_check(conf, iError, 'perturb')
   print *, 'perturb = ', prob%opts%pertube

   call aot_get_val(L = conf, key = 'perturb_s', val = prob%opts%pertube_s, ErrCode = iError)
   call error_check(conf, iError, 'perturb_s')
   print *, 'perturb_s = ', prob%opts%pertube_s
#endif

   call aot_get_val(L = conf, key = 'use_num_K', val = prob%opts%use_num_Kt, ErrCode = iError)
   call error_check(conf, iError, 'use_num_K')
   print *, 'use_num_K = ', prob%opts%use_num_Kt

   call aot_get_val(L = conf, key = 'use_num_D', val = prob%opts%use_num_Ct, ErrCode = iError)
   call error_check(conf, iError, 'use_num_D')
   print *, 'use_num_D = ', prob%opts%use_num_Ct

   call aot_get_val(L = conf, key = 'no_K', val = prob%opts%no_Kt, ErrCode = iError)
   call error_check(conf, iError, 'no_K')
   print *, 'no_K = ', prob%opts%no_Kt

   call aot_get_val(L = conf, key = 'no_D', val = prob%opts%no_Ct, ErrCode = iError)
   call error_check(conf, iError, 'no_D')
   print *, 'no_D = ', prob%opts%no_Ct

   call aot_get_val(L = conf, key = 'rtol', val = prob%opts%rtol, ErrCode = iError)
   call error_check(conf, iError, 'rtol')
   print *, 'rtol = ', prob%opts%rtol

   call aot_get_val(L = conf, key = 'atol', val = prob%opts%atol, ErrCode = iError)
   call error_check(conf, iError, 'atol')
   print *, 'atol = ', prob%opts%atol

   call aot_get_val(L = conf, key = 'imax', val = prob%opts%imax, ErrCode = iError)
   call error_check(conf, iError, 'imax')
   print *, 'imax = ', prob%opts%imax

   call aot_get_val(L = conf, key = 'stab2', val = prob%opts%stab2, ErrCode = iError)
   call error_check(conf, iError, 'imax')
   print *, 'stab2 = ', prob%opts%stab2

   ! Integration interval and step size
   call aot_get_val(L = conf, key = 't0', val = prob%opts%t0, ErrCode = iError)
   call error_check(conf, iError, 'imax')
   print *, 't0 = ', prob%opts%t0

   call aot_get_val(L = conf, key = 'te', val = prob%opts%te, ErrCode = iError)
   call error_check(conf, iError, 'imax')
   print *, 'te = ', prob%opts%te

   call aot_get_val(L = conf, key = 'steps', val = prob%opts%nsteps, ErrCode = iError)
   call error_check(conf, iError, 'imax')
   print *, 'steps = ', prob%opts%nsteps

   ! Problem options
   ! Mass of the pendulum
   call aot_get_val(L = conf, key = 'mass', val = prob%mass, ErrCode = iError)
   call error_check(conf, iError, 'mass')
   print *, 'mass = ', prob%mass
   if (prob%mass <= 0) then
      print *, 'FATAL Error: mass must be positive'
      ERROR STOP
   end if

   ! Lie group formulation
   call aot_get_val(L = conf, key = 'liegroup', val = prob%liegroup, ErrCode = iError)
   call error_check(conf, iError, 'liegroup')
   print *, 'liegroup = ', prob%liegroup
   if (prob%liegroup < 1 .or. prob%liegroup > 7) then
      print *, 'FATAL Error: liegroup must be in the range of 1,2,..,7'
      ERROR STOP
   end if
   if (prob%liegroup > 2) then
      prob%opts%constrained = 1
   else
      prob%opts%constrained = 0
   end if

   ! Gravity
   call aot_get_val(val = prob%gravity, ErrCode = iv3Error, L = conf, key = 'gravity')
   do i=1,3
      call error_check(conf, iv3Error(i), 'gravity')
   end do
   print *, 'gravity = [', prob%gravity, ']'

   ! Diagonal elements of inertial tensor wrt. the center of mass
   call aot_get_val(val = prob%inerJ, ErrCode = iv3Error, L = conf, key = 'inerJ')
   do i=1,3
      call error_check(conf, iv3Error(i), 'inerJ')
   end do
   print *, 'inerJ = [', prob%inerJ, ']'

   ! Initial configuration: x0
   call aot_get_val(val = prob%x0, ErrCode = iv3Error, L = conf, key = 'x0')
   do i=1,3
      call error_check(conf, iv3Error(i), 'x0')
   end do
   print *, 'x0 = [', prob%x0, ']'

   ! Initial configuration: p0
   call aot_get_val(val = prob%p0, ErrCode = iv4Error, L = conf, key = 'p0')
   do i=1,4
      call error_check(conf, iv4Error(i), 'p0')
   end do
   print *, 'p0 = [', prob%p0, ']'

   ! Initial configuration: Om0
   call aot_get_val(val = prob%Om0, ErrCode = iv3Error, L = conf, key = 'Om0')
   do i=1,3
      call error_check(conf, iv3Error(i), 'Om0')
   end do
   print *, 'Om0 = [', prob%Om0, ']'

   ! Output type
   call aot_get_val(L = conf, key = 'output_type', val = prob%output_type, ErrCode = iError)
   call error_check(conf, iError, 'output_type')
   print *, 'output_type = ', prob%output_type
   if (prob%output_type < 0 .or. prob%output_type > 2) then
      print *, 'FATAL Error: output_type must be in the range of 0,1,2'
      ERROR STOP
   end if

   ! Output only at certain times
   call aot_get_val(L = conf, key = 'output_t_at', val = prob%output_t_at, ErrCode = iError)
   call error_check(conf, iError, 'output_t_at')
   print *, 'output_t_at = ', prob%output_t_at

   if (prob%output_t_at == 1) then
      call aot_get_val(L = conf, key = 't_output_at_multiples_of', val = prob%t_output_at_multiples_of, ErrCode = iError)
      call error_check(conf, iError, 't_output_at_multiples_of')
      print *, 't_output_at_multiples_of = ', prob%t_output_at_multiples_of
   else
      print *, 't_output_at_multiples_of = n/a'
   end if

   call close_config(conf)

   ! flush stdout (output_unit is defined in the module iso_fortran_env)
   flush(output_unit)

   !!!!!! Prepare output file
   ! Get name of the output file
   call get_command_argument(2, prob%out_fname)
   if (len_trim(prob%out_fname) == 0) then
      print *, 'FATAL Error: Please provide an output file.'
      ERROR STOP
   elseif (len_trim(prob%out_fname) == 256) then
      print *, 'FATAL Error: Output filename is 256, characters long.'
      print *, '             Most likely we lost some characters of the'
      print *, '             filename. Please use a shorter file name.'
      ERROR STOP
   end if

   ! Open output file
   open(newunit = out_lua_lun,                 &
        file    = trim(prob%out_fname)//'.lua',&
        status  = 'new',                       &
        iostat  = iostat_number                )
   if (iostat_number /= 0) then
      print *, 'FATAL Error ', iostat_number
      print *, '  While creating the lua output file ', trim(prob%out_fname)//'.lua', '(does it exist already?)'
      ERROR STOP
   end if

   ! Header
   write (out_lua_lun, *) '-- ######################################################'
   write (out_lua_lun, *) '-- Original configuration file:'
   write (out_lua_lun, *) '-- ######################################################'

   ! Open config file and write it to the output file
   open(newunit=conf_lun, file=conf_fname)
   do
      call get_line(unit   = conf_lun,      &
                    line   = conf_line,     &
                    iostat = iostat_number, &
                    iomsg  = iomsg_str      )
      if (is_iostat_end(iostat_number)) then
         exit
      elseif (iostat_number /= 0) then
         print *, 'FATAL Error reading config file:', iostat_number, iomsg_str
         close(out_lua_lun)
         close(conf_lun)
         ERROR STOP
      end if
      write (out_lua_lun, *) conf_line
   end do
   ! Close config file
   close(conf_lun)

   ! Header for the results:
   write (out_lua_lun, *) ''
   write (out_lua_lun, *) ''
   write (out_lua_lun, *) '-- ######################################################'
   write (out_lua_lun, *) '-- Results of the integration:'
   write (out_lua_lun, *) '-- ######################################################'
   write (out_lua_lun, *) '-- Using heavy_top which was compiled on ', __DATE__, ' at ', __TIME__
   write (out_lua_lun, *) '-- '
   write (out_lua_lun, *) '-- We used the following integrator'
   ! This is FIXME, a really ugly way of stringification
   write (out_lua_lun, '(AAAA)') ' integrator = ', '"', '&
   INTEGRATOR', '"'
   ! EMXIF
   write (out_lua_lun, *) ''

   ! Flush lua file
   flush(out_lua_lun)

   ! Open binary output file
   open(newunit = prob%out_bin_lun,            &
        file    = trim(prob%out_fname)//'.bin',&
        form    = 'unformatted',               &
        access  = 'stream',                    &
        status  = 'new',                       &
        iostat  = iostat_number                )
   if (iostat_number /= 0) then
      print *, 'FATAL Error ', iostat_number
      print *, '  While creating the binary output file ', trim(prob%out_fname)//'.bin', '(does it exist already?)'
      close(out_lua_lun)
      ERROR STOP
   end if

   ! Open misc output file
   open(newunit = prob%out_misc_lun,            &
        file    = trim(prob%out_fname)//'.misc',&
        status  = 'new',                        &
        iostat  = iostat_number                 )
   if (iostat_number /= 0) then
      print *, 'FATAL Error ', iostat_number
      print *, '  While creating the misc output file ', trim(prob%out_fname)//'.lua', '(does it exist already?)'
      ERROR STOP
   end if

   !!!!!! Start the actual integration
   call prob%GL(INTEGRATOR)_integrate()

   ! Close binary and misc output files
   close(prob%out_bin_lun)
   close(prob%out_misc_lun)

   ! Write statistics
   write (out_lua_lun, *) 'cpu_time = ', prob%GL(INTEGRATOR)_stats%time
   write (out_lua_lun, *) 'newt_steps_sum = ', prob%GL(INTEGRATOR)_stats%newt_steps_sum
   write (out_lua_lun, *) 'newt_steps_max = ', prob%GL(INTEGRATOR)_stats%newt_steps_max
   write (out_lua_lun, *) 'newt_steps_avg = ', prob%GL(INTEGRATOR)_stats%newt_steps_avg
   write (out_lua_lun, *) 'n_g_calls = ', prob%GL(INTEGRATOR)_stats%ngcalls
   write (out_lua_lun, *) 'n_B_calls = ', prob%GL(INTEGRATOR)_stats%nBcalls

   ! Show statistics
   print *, ''
   call prob%GL(INTEGRATOR)_print_stats

   ! Clean up
   call prob%GL(INTEGRATOR)_cleanup()

   ! Close files
   close(out_lua_lun)

   ! Done
   print *, ''
   print *, 'Done'

contains

   subroutine error_check(conf, iError, key)
      use aotus_module, only: aoterr_Fatal, aoterr_WrongType, aoterr_NonExistent
      implicit none
      type(flu_State),  intent(inout)  :: conf
      integer,          intent(in   )  :: iError
      character(len=*), intent(in   )  :: key
      !
      if (btest(iError, aoterr_Fatal)) then
         write(*,*) 'FATAL Error occured, while retrieving variable ', key, ':', iError
         if (btest(iError, aoterr_NonExistent)) write(*,*) 'Variable not existent!'
         if (btest(iError, aoterr_WrongType)) write(*,*) 'Variable has wrong type!'
         ERROR STOP
      end if
   end subroutine error_check

end program main
