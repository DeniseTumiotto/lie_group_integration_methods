module get_line_of_variable_length

contains

   subroutine get_line(unit, line, iostat, iomsg)
      ! Input/output
      integer,          intent(in   )              :: unit
      character(len=:), intent(  out), allocatable :: line
      integer,          intent(  out)              :: iostat
      character(*),     intent(inout)              :: iomsg
      ! Internal variables
      integer, parameter            :: buffer_len = 64
      character(len=buffer_len)     :: buffer
      integer                       :: size_read
      
      ! Line is initially empty (notice implicit allocation)
      line = ''
      do 
         ! Read chunks of size buffer_len of the line
         read( unit    = unit,     &
               fmt     = '(A)',    &
               iostat  = iostat,   &
               iomsg   = iomsg,    &
               advance = 'no',     &
               size    = size_read &
              ) buffer
         ! Check if we reached the end of the line
         if (is_iostat_eor(iostat)) then
            ! Update line (notice implicit reallocation
            line = line // buffer(:size_read)
            ! Set iostat = 0, because we finished without errors
            iostat = 0
            exit
         else if (iostat == 0) then
            ! Successfully read a full buffer, update line
            ! (notice implicit reallocation)
            line = line // buffer
         else
            ! Something went wrong, exit, where iomsg is the
            ! iomsg that read returned
            exit
         end if
      end do
   end subroutine get_line
   
end module get_line_of_variable_length
