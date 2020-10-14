   subroutine read_menu_input(file_name, file_unit, input_section_in, input_found)

      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(*), intent(in)  :: file_name
      integer,      intent(in)  :: file_unit
      character(*), intent(in)  :: input_section_in
      logical,      intent(out) :: input_found
!     --------------------------------------------------------------------------
      character(kw_length)      :: input_section
      character(kw_length)      :: word
      character(kw_length)      :: kw_section
      logical                   :: read_all_input
      logical                   :: file_is_open
!     --------------------------------------------------------------------------

!     this is to catch keywords that appear before some section starts
      kw_section = '       '

!     check for the input section we would like to read
      input_section = uppercase(input_section_in) ! truncate or extend to kw_length

!     if applicable process the complete input (all modules)
      if (input_section(1:3) == 'ALL') then
         read_all_input = .true.
      else
         read_all_input = .false.
      end if

      inquire(file=file_name, opened=file_is_open)
      if (.not. file_is_open) then
         open(file_unit, file=file_name)
      end if

      rewind(file_unit)
      call set_file_unit(file_unit)

!     initialize
      input_found = .false.

      do while (.true.)

        read(file_unit, '(a7)', end=1) word

        if (word == '       ') then
!           blank line
        else
           select case (word(1:1))

             case ('!', '#')
!               comment

             case ('*')
!               section
                kw_section = uppercase(word)
                if (word == '*END OF' .or. word == '**END O') go to 1

             case default
!               keyword
                if (read_all_input .or. kw_section == input_section) then
                   input_found = .true.
                   call read_input_sections(word, kw_section)
                end if

           end select
        end if

      end do

1     if (.not. file_is_open) then
         close(file_unit, status='keep')
      end if

   end subroutine

   subroutine move_to_next_star(word_io, file_unit)

      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(inout) :: word_io
      integer,              intent(in)    :: file_unit
!     --------------------------------------------------------------------------
      character(kw_length)                :: word
!     --------------------------------------------------------------------------

      do while (.true.)
         read(file_unit, '(a7)') word
         if (word(1:1) == '*') then
            word_io = word
            return
         end if
      end do

   end subroutine
