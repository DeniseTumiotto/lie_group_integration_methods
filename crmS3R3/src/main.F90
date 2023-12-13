! 2016-01-11: Created by moving program main form crm.f90
!             Removed testdeom
! 2016-01-13: Cleaned up and completed main.F90

#define GL(x) x

program main
   !!!!!! Use modules
   ! Use modul to be able to identify the output unit
   use iso_fortran_env
   ! Use the module with the actual crm problem
   use testprobcrm
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
   type(myproblem)      :: prob
   ! For looping
   integer              :: i, j
   ! Aotus handles and error variables, etc.
   character(len=256)   :: conf_fname
   type(flu_State)      :: conf
   integer              :: iError
   integer              :: ivError(3)
   integer              :: iqError(4)
   integer, allocatable :: ivvError(:)
   character(len=256)   :: cError
   integer, parameter   :: max_length = 256
   type(aot_fun_type)   :: x0_fun
   type(aot_fun_type)   :: V0_fun
   type(aot_fun_type)   :: p0_fun
   type(aot_fun_type)   :: Om0_fun
   integer              :: external_parameters ! handle to a table
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

#ifdef INT_gena
   ! Algorithmic parameters
   call aot_get_val(L = conf, key = 'alpha_m', val = prob%gena_alpha_m, ErrCode = iError)
   call error_check(conf, iError, 'alpha_m')
   print *, 'alpha_m = ', prob%gena_alpha_m

   call aot_get_val(L = conf, key = 'alpha_f', val = prob%gena_alpha_f, ErrCode = iError)
   call error_check(conf, iError, 'alpha_f')
   print *, 'alpha_f = ', prob%gena_alpha_f

   call aot_get_val(L = conf, key = 'beta', val = prob%gena_beta, ErrCode = iError)
   call error_check(conf, iError, 'beta')
   print *, 'beta = ', prob%gena_beta

   call aot_get_val(L = conf, key = 'gamma', val = prob%gena_gamma, ErrCode = iError)
   call error_check(conf, iError, 'gamma')
   print *, 'gamma = ', prob%gena_gamma
#endif

#ifdef INT_BLieDF
   ! Algorithmic parameters
   call aot_get_val(L = conf, key = 'k_bdf', val = prob%k_bdf, ErrCode = iError)
   call error_check(conf, iError, 'k_bdf')
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

   call aot_get_val(L = conf, key = 'nr_subdiag', val = prob%opts%nr_subdiag, ErrCode = iError)
   call error_check(conf, iError, 'nr_subdiag')
   print *, 'nr_subdiag = ', prob%opts%nr_subdiag

   call aot_get_val(L = conf, key = 'nr_superdiag', val = prob%opts%nr_superdiag, ErrCode = iError)
   call error_check(conf, iError, 'nr_superdiag')
   print *, 'nr_superdiag = ', prob%opts%nr_superdiag

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

#ifdef VARIABLE_STEPS
   if (allocated(ivvError)) deallocate(ivvError)
   call aot_get_val(L = conf, key = 'tspan', val = prob%opts%tspan, MaxLength = 4*prob%opts%nsteps, ErrCode = ivvError)
   do i=1,size(ivvError); call error_check(conf, ivvError(i), 'tspan'); end do
   if ( size(prob%opts%tspan) == 4*prob%opts%nsteps ) ERROR STOP "tspan is not the right length"
   print *, 'tspan = ', prob%opts%tspan(1), prob%opts%tspan(2), '...', prob%opts%tspan(size(prob%opts%tspan)-2), prob%opts%tspan(size(prob%opts%tspan)-1), prob%opts%tspan(size(prob%opts%tspan))
#endif

   ! Problem options
   call aot_get_val(L = conf, key = 'problem_name', val = prob%probname, ErrCode = iError)
   call error_check(conf, iError, 'problem_name')
   print *, 'problem_name = ', trim(prob%probname)

   call aot_get_val(L = conf, key = 'kirchhoff', val = prob%kirchhoff, ErrCode = iError)
   call error_check(conf, iError, 'kirchhoff')
   print *, 'kirchhoff = ', prob%kirchhoff

!   call aot_get_val(L = conf, key = 'inextensible', val = prob%inextensible, ErrCode = iError)
!   call error_check(conf, iError, 'inextensible')
!   print *, 'inextensible = ', prob%inextensible
!
   ! Constrained problem if we use a Kirchhoff or inextensible beam model
!   if (prob%inextensible == 1 .or. prob%kirchhoff == 1) then
   if (prob%kirchhoff == 1) then
      prob%opts%constrained = 1
   else
      prob%opts%constrained = 0
   end if

!   call aot_get_val(L = conf, key = 'use_SMA', val = prob%use_SMA, ErrCode = iError)
!   call error_check(conf, iError, 'use_SMA')
!   print *, 'use_SMA = ', prob%use_SMA

   call aot_get_val(L = conf, key = 'n', val = prob%n, ErrCode = iError)
   call error_check(conf, iError, 'n')
   print *, 'n = ', prob%n

   prob%n1 = prob%n + 1

   call aot_get_val(L = conf, key = 'fixed_x0', val = prob%fixed_x0, ErrCode = iError)
   call error_check(conf, iError, 'fixed_x0')
   print *, 'fixed_x0 = ', prob%fixed_x0

   if (prob%fixed_x0 == 1) then
      call aot_get_val(L = conf, key = 'fixed_x0_position', val = prob%fixed_x0_position, ErrCode = ivError)
      do i=1,3; call error_check(conf, ivError(i), 'fixed_x0_position'); end do
      print *, 'fixed_x0_position = ', prob%fixed_x0_position
   end if

   call aot_get_val(L = conf, key = 'fixed_xn', val = prob%fixed_xN, ErrCode = iError)
   call error_check(conf, iError, 'fixed_xn')
   print *, 'fixed_xn = ', prob%fixed_xn

   if (prob%fixed_xn == 1) then
      call aot_get_val(L = conf, key = 'fixed_xn_position', val = prob%fixed_xn_position, ErrCode = ivError)
      do i=1,3; call error_check(conf, ivError(i), 'fixed_xn_position'); end do
      print *, 'fixed_xn_position = ', prob%fixed_xn_position
   end if

   call aot_get_val(L = conf, key = 'fixed_p0', val = prob%fixed_p0, ErrCode = iError)
   call error_check(conf, iError, 'fixed_p0')
   print *, 'fixed_p0 = ', prob%fixed_p0

   if (prob%fixed_p0 == 1) then
      call aot_get_val(L = conf, key = 'fixed_p0_orientation', val = prob%fixed_p0_orientation, ErrCode = iqError)
      do i=1,4; call error_check(conf, iqError(i), 'fixed_p0_orientation'); end do
      print *, 'fixed_p0_orientation = ', prob%fixed_p0_orientation
   end if

   call aot_get_val(L = conf, key = 'fixed_pn', val = prob%fixed_pn, ErrCode = iError)
   call error_check(conf, iError, 'fixed_pn')
   print *, 'fixed_pn = ', prob%fixed_pn

   if (prob%fixed_pn == 1) then
      call aot_get_val(L = conf, key = 'fixed_pn_orientation', val = prob%fixed_pn_orientation, ErrCode = iqError)
      do i=1,4; call error_check(conf, iqError(i), 'fixed_pn_orientation'); end do
      print *, 'fixed_pn_orientation = ', prob%fixed_pn_orientation
   end if

   call aot_get_val(L = conf, key = 'ds', val = prob%ds, ErrCode = iError)
   call error_check(conf, iError, 'ds')
   print *, 'ds = ', prob%ds

   call aot_get_val(L = conf, key = 'mI', val = prob%mI, ErrCode = ivError)
   do i=1,3; call error_check(conf, ivError(i), 'mI'); end do
   print *, 'mI = ', prob%mI

   call aot_get_val(L = conf, key = 'm', val = prob%m, ErrCode = iError)
   call error_check(conf, iError, 'm')
   print *, 'm = ', prob%m

   call aot_get_val(L = conf, key = 'CGam', val = prob%CGam, ErrCode = ivError)
   do i=1,3; call error_check(conf, ivError(i), 'CGam'); end do
   print *, 'CGam = ', prob%CGam

   call aot_get_val(L = conf, key = 'CGamd', val = prob%CGamd, ErrCode = ivError)
   do i=1,3; call error_check(conf, ivError(i), 'CGamd'); end do
   print *, 'CGamd = ', prob%CGamd

   call aot_get_val(L = conf, key = 'CK', val = prob%CK, ErrCode = ivError)
   do i=1,3; call error_check(conf, ivError(i), 'CK'); end do
   print *, 'CK = ', prob%CK

   call aot_get_val(L = conf, key = 'CKd', val = prob%CKd, ErrCode = ivError)
   do i=1,3; call error_check(conf, ivError(i), 'CKd'); end do
   print *, 'CKd = ', prob%CKd

   ! Initial configuration
   call aot_fun_open(L = conf, fun = x0_fun, key = 'x0')
   allocate(prob%x0(3,prob%n1))
   print *, 'x0 = '
   do i=1,prob%n1
      call aot_fun_put(L = conf, fun = x0_fun, arg = (i-1.0_8)/(prob%N1-1))
      call aot_fun_do(L = conf, fun = x0_fun, nresults = 1)
      call aot_top_get_val(L = conf, val = prob%x0(:,i), ErrCode = ivError)
      do j=1,3; call error_check(conf, ivError(j), 'x0'); end do
      print *, prob%x0(:,i)
   end do
   call aot_fun_close(L = conf, fun = x0_fun)

   call aot_fun_open(L = conf, fun = V0_fun, key = 'V0')
   allocate(prob%V0(3,prob%n1))
   print *, 'V0 = '
   do i=1,prob%n1
      call aot_fun_put(L = conf, fun = V0_fun, arg = (i-1.0_8)/(prob%n1-1))
      call aot_fun_do(L = conf, fun = V0_fun, nresults = 1)
      call aot_top_get_val(L = conf, val = prob%V0(:,i), ErrCode = ivError)
      do j=1,3; call error_check(conf, ivError(j), 'V0'); end do
      print *, prob%V0(:,i)
   end do
   call aot_fun_close(L = conf, fun = V0_fun)

   call aot_fun_open(L = conf, fun = p0_fun, key = 'p0')
   allocate(prob%p0(4,prob%n1))
   print *, 'p0 = '
   do i=1,prob%n1
      call aot_fun_put(L = conf, fun = p0_fun, arg = (i-1.0_8)/(prob%n1-1))
      call aot_fun_do(L = conf, fun = p0_fun, nresults = 1)
      call aot_top_get_val(L = conf, val = prob%p0(:,i), ErrCode = iqError)
      do j=1,4; call error_check(conf, iqError(j), 'p0'); end do
      print *, prob%p0(:,i)
   end do
   call aot_fun_close(L = conf, fun = p0_fun)

   call aot_fun_open(L = conf, fun = Om0_fun, key = 'Om0')
   allocate(prob%Om0(3,prob%n1))
   print *, 'Om0 = '
   do i=1,prob%N1
      call aot_fun_put(L = conf, fun = Om0_fun, arg = (i-1.0_8)/(prob%n1-1))
      call aot_fun_do(L = conf, fun = Om0_fun, nresults = 1)
      call aot_top_get_val(L = conf, val = prob%Om0(:,i), ErrCode = ivError)
      do j=1,3; call error_check(conf, ivError(j), 'Om0'); end do
      print *, prob%Om0(:,i)
   end do
   call aot_fun_close(L = conf, fun = Om0_fun)

   ! Output options
   call aot_get_val(L = conf, key = 'output_s_at', val = prob%output_s_at, ErrCode = iError)
   call error_check(conf, iError, 'output_s_at')
   print *, 'output_s_at = ', prob%output_s_at

   if (prob%output_s_at == 1) then
      if (allocated(ivvError)) deallocate(ivvError)
      call aot_get_val(L = conf, key = 'output_s', val = prob%output_s, MaxLength = max_length, ErrCode = ivvError)
      do i=1,size(ivvError); call error_check(conf, ivvError(i), 'output_s'); end do
      if (size(prob%output_s) == max_length) print *, 'WARNING Variable output_s may be truncated'
      print *, 'output_s = ', prob%output_s
   end if
!
!   if (allocated(ivvError))  deallocate(ivvError)
!   call aot_get_val(L = conf, key = 's_output_p', val = prob%s_output_p, MaxLength = max_length, ErrCode = ivvError)
!   do i=1,size(ivvError); call error_check(conf, ivvError(i), 's_output_p'); end do
!   if (size(prob%s_output_p) == max_length) print *, 'WARNING Variable s_output_p may be truncated'
!   print *, 's_output_p = ', prob%s_output_p

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

   ! External forces and moments
   call aot_get_val(L = conf, key = 'external', val = prob%externalfm_name, ErrCode = iError)
   call error_check(conf, iError, 'external')
   print *, 'external = ', trim(prob%externalfm_name)

   call aot_table_open(L = conf, key = 'external_parameters', thandle = external_parameters)
   select case (trim(prob%externalfm_name))
      case ('gravity')
         prob%externalfm_type = 1
         call aot_get_val(L = conf, thandle = external_parameters, key = 'g', val = prob%externalfm_gravity_g, ErrCode = iError)
         call error_check(conf, iError, 'gravity.g')
         print *, 'gravity.g = ', prob%externalfm_gravity_g
      case ('flying_spaghetti')
         prob%externalfm_type = 2

         call aot_get_val(L = conf, thandle = external_parameters, key = 'increasing_time', val = prob%externalfm_flying_spaghetti_increasing_time, ErrCode = iError)
         call error_check(conf, iError, 'flying_spaghetti.increasing_time')
         print *, 'flying_spaghetti.increasing_time = ', prob%externalfm_flying_spaghetti_increasing_time

         call aot_get_val(L = conf, thandle = external_parameters, key = 'decreasing_time', val = prob%externalfm_flying_spaghetti_decreasing_time, ErrCode = iError)
         call error_check(conf, iError, 'flying_spaghetti.decreasing_time')
         print *, 'flying_spaghetti.decreasing_time = ', prob%externalfm_flying_spaghetti_decreasing_time

         call aot_get_val(L = conf, thandle = external_parameters, key = 'maximum_height', val = prob%externalfm_flying_spaghetti_maximum_height, ErrCode = iError)
         call error_check(conf, iError, 'flying_spaghetti.maximum_height')
         print *, 'flying_spaghetti.maximum_height = ', prob%externalfm_flying_spaghetti_maximum_height

         call aot_get_val(L = conf, thandle = external_parameters, key = 'force_factors', val = prob%externalfm_flying_spaghetti_force_factors, ErrCode = ivError)
         do i=1,3; call error_check(conf, ivError(i), 'flying_spaghetti.force_factors'); end do
         print *, 'flying_spaghetti.force_factors = ', prob%externalfm_flying_spaghetti_force_factors

         call aot_get_val(L = conf, thandle = external_parameters, key = 'moment_factors', val = prob%externalfm_flying_spaghetti_moment_factors, ErrCode = ivError)
         do i=1,3; call error_check(conf, ivError(i), 'flying_spaghetti.moment_factors'); end do
         print *, 'flying_spaghetti.moment_factors = ', prob%externalfm_flying_spaghetti_moment_factors
      case ('roll-up')
         prob%externalfm_type = 3
         call aot_get_val(L = conf, thandle = external_parameters, key = 'factor', val = prob%externalfm_roll_up_factor, ErrCode = iError)
         call error_check(conf, iError, 'roll_up.moment_factors')
         print *, 'roll_up.factor = ', prob%externalfm_roll_up_factor

      case default
         print *, 'FATAL Error: unknown external forces and moments'
         call close_config(conf)
         errorstop
   end select

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
   write (out_lua_lun, *) '-- Using crmS3R3 which was compiled on ', __DATE__, ' at ', __TIME__
   write (out_lua_lun, *) '-- '
   write (out_lua_lun, *) '-- We used the following integrator'
   ! This is FIXME, a really ugly way of stringification
   write (out_lua_lun, '(AAAA)') ' integrator = ', '"', '&
   INTEGRATOR', '"'
   ! EMXIF
   write (out_lua_lun, *) ''
#ifdef VARIABLE_STEPS
   write (out_lua_lun, *) 'variable_steps = 1'
#else
   write (out_lua_lun, *) 'variable_steps = 0'
#endif

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

   ! write everything to disk
   flush(prob%out_misc_lun)
   flush(prob%out_bin_lun)
   flush(out_lua_lun)

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
