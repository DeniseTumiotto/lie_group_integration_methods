#define GL(x) x

program main
   !!!!!! Use modules
   ! Use modul to be able to identify the output unit
   use iso_fortran_env
   ! Use the module with the actual rolling_disk problem
   use testprobrolling_disk
   ! For reading the config file
   use aotus_module, only: flu_State, open_config_file, close_config, aot_get_val, &
                           aot_top_get_val, &
                           aoterr_Fatal, aoterr_WrongType, aoterr_NonExistent
   use aot_table_module, only:   aot_table_open, aot_table_close, aot_table_length
   use aot_fun_module, only:  aot_fun_type, aot_fun_open, aot_fun_put, &
                              aot_fun_do, aot_fun_close
   ! For reading lines of variable length
   use get_line_of_variable_length
   ! For quaternion functions
   use quaternion_functions

   !!!!!! No implicit variables
   implicit none

   !!!!!! Variables
   ! rolling_disk problem object
   type(myproblem)      :: prob
   ! For looping
   integer              :: i, j
   ! Aotus handles and error variables, etc.
   character(len=256)   :: conf_fname
   type(flu_State)      :: conf
   integer              :: iError
   integer              :: ivError(3)
   integer              :: ipError(4)
   character(len=256)   :: cError
   integer, parameter   :: max_length = 256
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
      errorstop
   end if

   !!!!!! Read config file
   ! Open file
   call open_config_file(L = conf, filename = conf_fname, &
                         ErrCode = iError, ErrString = cError)
   if (iError /= 0) then
      write(*,*) 'FATAL Error ', iError, ' when opening the Lua config file:', cError
      STOP
   end if

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

   call aot_get_val(L = conf, key = 'nr_subdiag', val = prob%opts%nr_subdiag, ErrCode = iError)
   call error_check(conf, iError, 'nr_subdiag')
   print *, 'nr_subdiag = ', prob%opts%nr_subdiag

   call aot_get_val(L = conf, key = 'nr_superdiag', val = prob%opts%nr_superdiag, ErrCode = iError)
   call error_check(conf, iError, 'nr_superdiag')
   print *, 'nr_superdiag = ', prob%opts%nr_superdiag

   call aot_get_val(L = conf, key = 'recalc_iteration_matrix', val = prob%opts%recalc_iteration_matrix, ErrCode = iError)
   call error_check(conf, iError, 'recalc_iteration_matrix')
   print *, 'recalc_iteration_matrix = ', prob%opts%recalc_iteration_matrix

   call aot_get_val(L = conf, key = 'rtol', val = prob%opts%rtol, ErrCode = iError)
   call error_check(conf, iError, 'rtol')
   print *, 'rtol = ', prob%opts%rtol

   call aot_get_val(L = conf, key = 'atol', val = prob%opts%atol, ErrCode = iError)
   call error_check(conf, iError, 'atol')
   print *, 'atol = ', prob%opts%atol

   call aot_get_val(L = conf, key = 'imax', val = prob%opts%imax, ErrCode = iError)
   call error_check(conf, iError, 'imax')
   print *, 'imax = ', prob%opts%imax

   prob%opts%stab2 = 1

   ! Integration interval and step size
   call aot_get_val(L = conf, key = 't0', val = prob%opts%t0, ErrCode = iError)
   call error_check(conf, iError, 't0')
   print *, 't0 = ', prob%opts%t0

   call aot_get_val(L = conf, key = 'te', val = prob%opts%te, ErrCode = iError)
   call error_check(conf, iError, 'te')
   print *, 'te = ', prob%opts%te

   call aot_get_val(L = conf, key = 'steps', val = prob%opts%nsteps, ErrCode = iError)
   call error_check(conf, iError, 'steps')
   print *, 'steps = ', prob%opts%nsteps

   ! Problem options and initial configuration
   call aot_get_val(L = conf, key = 'problem_name', val = prob%probname, ErrCode = iError)
   call error_check(conf, iError, 'problem_name')
   print *, 'problem_name = ', trim(prob%probname)

   prob%opts%constrained = 1

   call aot_get_val(L = conf, key = 'mass', val = prob%mass, ErrCode = iError)
   call error_check(conf, iError, 'mass')
   print *, 'mass = ', prob%mass

   call aot_get_val(L = conf, key = 'gravity', val = prob%gravity, ErrCode = ivError)
   do i=1,3; call error_check(conf, ivError(i), 'gravity'); end do
   print *, 'gravity = ', prob%gravity

   call aot_get_val(L = conf, key = 'inerJ', val = prob%inerJ, ErrCode = ivError)
   do i=1,3; call error_check(conf, ivError(i), 'inerJ'); end do
   print *, 'inerJ = ', prob%inerJ

   call aot_get_val(L = conf, key = 'p0', val = prob%p0, ErrCode = ipError)
   do i=1,4; call error_check(conf, ipError(i), 'p0'); end do
   print *, 'p0 = ', prob%p0

   call aot_get_val(L = conf, key = 'x0', val = prob%x0, ErrCode = ivError)
   do i=1,3; call error_check(conf, ivError(i), 'x0'); end do
   print *, 'x0 = ', prob%x0

   call aot_get_val(L = conf, key = 'Om0', val = prob%Om0, ErrCode = ivError)
   do i=1,3; call error_check(conf, ivError(i), 'Om0'); end do
   print *, 'Om0 = ', prob%Om0

   ! Output options
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
      errorstop
   elseif (len_trim(prob%out_fname) == 256) then
      print *, 'FATAL Error: Output filename is 256, characters long.'
      print *, '             Most likely we lost some characters of the'
      print *, '             filename. Please use a shorter file name.'
      errorstop
   end if

   ! Open output file
   open(newunit = out_lua_lun,                 &
        file    = trim(prob%out_fname)//'.lua',&
        status  = 'new',                       &
        iostat  = iostat_number                )
   if (iostat_number /= 0) then
      print *, 'FATAL Error ', iostat_number
      print *, '  While creating the lua output file ', trim(prob%out_fname)//'.lua', '(does it exist already?)'
      STOP
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
         STOP
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
   write (out_lua_lun, *) '-- Using rolling_disk which was compiled on ', __DATE__, ' at ', __TIME__
   write (out_lua_lun, *) '-- '
   write (out_lua_lun, *) '-- We used the following integrator'
   ! This is FIXME, a really ugly way of stringification
   write (out_lua_lun, '(AAAA)') ' integrator = ', '"', '&
   INTEGRATOR', '"'
   ! EMXIF

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
      STOP
   end if

   ! Open misc output file
   open(newunit = prob%out_misc_lun,            &
        file    = trim(prob%out_fname)//'.misc',&
        status  = 'new',                        &
        iostat  = iostat_number                 )
   if (iostat_number /= 0) then
      print *, 'FATAL Error ', iostat_number
      print *, '  While creating the misc output file ', trim(prob%out_fname)//'.lua', '(does it exist already?)'
      STOP
   end if

   !!!!!! Start the actual integration
   call prob%GL(INTEGRATOR)_integrate()

   ! Close binary and misc output files
   close(prob%out_bin_lun)
   close(prob%out_misc_lun)

   ! Write statistics
   write (out_lua_lun, *) 'cpu_time = ', prob%GL(INTEGRATOR)_stats%time
   write (out_lua_lun, *) 'newt_steps_max = ', prob%GL(INTEGRATOR)_stats%newt_steps_max
   write (out_lua_lun, *) 'newt_steps_avg = ', prob%GL(INTEGRATOR)_stats%newt_steps_avg
   write (out_lua_lun, *) 'n_g_calls = ', prob%GL(INTEGRATOR)_stats%ngcalls
   write (out_lua_lun, *) 'n_B_calls = ', prob%GL(INTEGRATOR)_stats%nBcalls

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
         errorstop
      end if
   end subroutine error_check

end program main
