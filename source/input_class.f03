! input class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module input_class

      use perrors_class
      use parallel_class
      use parallel_pipe_class
      use spect3d_class
      use mpi
      use json_module
         
      implicit none

      private

      public :: input_json
      
      type input_json

         class(spect3d), pointer :: sp => null()
         class(perrors), pointer :: err => null()
         class(parallel), pointer :: p => null()
         class(parallel_pipe), pointer :: pp => null()
         type(json_file), private, pointer :: input => null()

         contains
         
         generic :: new => read_input_json
         generic :: get => json_file_get_object,json_file_get_integer,&
         &json_file_get_double, json_file_get_logical,&
         &json_file_get_string, json_file_get_integer_vec,&
         &json_file_get_double_vec,json_file_get_logical_vec,&
         &json_file_get_string_vec,json_file_get_alloc_string_vec,&
         &json_file_get_root
         generic :: info => json_file_variable_info

         procedure :: found
         procedure, private :: read_input_json
         procedure, private :: initialize, set_json_core_in_file
         procedure, private :: load_file, print_to_string
         procedure, private :: load_from_string
         procedure, private :: json_file_get_object,json_file_get_integer,&
         &json_file_get_double, json_file_get_logical,&
         &json_file_get_string, json_file_get_integer_vec,&
         &json_file_get_double_vec,json_file_get_logical_vec,&
         &json_file_get_string_vec,json_file_get_alloc_string_vec,&
         &json_file_get_root
         procedure, private :: json_file_variable_info
         
      end type
      
      character(len=10), save :: class = 'input:'
      character(len=128), save :: erstr
      type(spect3d), save, target  :: sp
      type(perrors), save, target :: err
      type(parallel), save, target :: p
      type(parallel_pipe), save, target :: pp
     
      contains
!      
      subroutine read_input_json(this)
    
         implicit none
         
         class(input_json), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'read_input_json:'
         logical :: found, stat
         character(len=:), allocatable :: ff, boundary, error_msg
         integer :: length, num_stages, verbose, indx, indy, indz, psolve, ierr
         
         call p%new()
         
         this%p => p

         call err%new(this%p,2,monitor=0)
         
         this%err => err
         
         call this%err%werrfl0(class//sname//' started')
         
         call this%initialize(comment_char='!')
         
         if (p%getidproc() == 0) then
            inquire(FILE='./qpinput.json', EXIST=found)
            if(found) then
! read the file
               call this%load_file(filename = './qpinput.json')
               call this%input%check_for_errors(stat,error_msg)
               if (.not. stat) then
                  call this%err%equit(class//sname//trim(error_msg))
                  return
               end if
            else
                write (erstr,*) 'error: cannot find the input file'
                call this%err%equit(class//sname//erstr)
                return
            end if
            call this%print_to_string(ff)
            length = len(ff)
         end if

         call MPI_BCAST(length, 1, this%p%getmint(), 0, this%p%getlworld(), ierr)
         
         if (.not. allocated(ff)) allocate(character(len=length) :: ff)
         
         call MPI_BCAST(ff, length, this%p%getmchar(), 0, this%p%getlworld(), ierr)
         
         call this%load_from_string(ff)
         
         call this%get('simulation.nodes(2)',num_stages)
         
         call pp%new(nst=num_stages)
         
         this%pp => pp
                     
         call this%get('simulation.indx',indx)

         call this%get('simulation.indy',indy)

         call this%get('simulation.indz',indz)

         call this%get('simulation.boundary',boundary)
         
         select case (trim(boundary))
         case ("conducting")
            psolve = 1
         case default
            psolve = 1
         end select
         
         call sp%new(this%pp,this%err,indx,indy,indz,psolve,1)

         this%sp => sp

         call this%err%werrfl0(ff)

         call this%err%werrfl0(class//sname//' ended')
         
         call this%get('simulation.verbose',verbose)

         call err%setmonitor(verbose)

      end subroutine read_input_json
!
      subroutine initialize(this,verbose,compact_reals,print_signs,&
      &real_format,spaces_per_tab,strict_type_checking,trailing_spaces_significant,&
      &case_sensitive_keys,no_whitespace,unescape_strings,comment_char,path_mode,&
      &path_separator,compress_vectors,allow_duplicate_keys)
    
         implicit none
         
         class(input_json), intent(inout) :: this
         logical,intent(in),optional :: verbose
         logical,intent(in),optional :: compact_reals
         logical,intent(in),optional :: print_signs
         character(len=*),intent(in),optional :: real_format
         integer,intent(in),optional :: spaces_per_tab
         logical,intent(in),optional :: strict_type_checking
         logical,intent(in),optional :: trailing_spaces_significant
         logical,intent(in),optional :: case_sensitive_keys
         logical,intent(in),optional :: no_whitespace
         logical,intent(in),optional :: unescape_strings
         character(len=1),intent(in),optional :: comment_char
         integer,intent(in),optional :: path_mode
         character(len=1),intent(in),optional :: path_separator
         logical,intent(in),optional :: compress_vectors
         logical,intent(in),optional :: allow_duplicate_keys
! local data
         character(len=38), save :: sname = 'initialize_json_core_in_file:'
         call this%err%werrfl2(class//sname//' started')
         
         allocate(this%input)
         call this%input%initialize(verbose,compact_reals,print_signs,&
         &real_format,spaces_per_tab,strict_type_checking,trailing_spaces_significant,&
         &case_sensitive_keys,no_whitespace,unescape_strings,comment_char,path_mode,&
         &path_separator,compress_vectors,allow_duplicate_keys)

         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine initialize
!
      subroutine set_json_core_in_file(this,core)

         implicit none

         class(input_json),intent(inout) :: this
         type(json_core),intent(in) :: core
! local data
         character(len=38), save :: sname = 'set_json_core_in_file:'
         call this%err%werrfl2(class//sname//' started')

         call this%input%initialize(core)

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine set_json_core_in_file   
!
      subroutine load_file(this, filename, unit)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: filename
         integer,intent(in),optional :: unit
! local data
         character(len=18), save :: sname = 'json_file_load:'
         call this%err%werrfl2(class//sname//' started')

         call this%input%load_file(filename, unit)

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine load_file
!
      subroutine print_to_string(this, str)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=:),allocatable,intent(out) :: str
! local data
         character(len=38), save :: sname = 'json_file_print_to_string:'
         call this%err%werrfl2(class//sname//' started')

         call this%input%print_to_string(str)

         call this%err%werrfl2(class//sname//' ended')

      end subroutine print_to_string     
!
      subroutine load_from_string(this, str)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: str
! local data
         character(len=38), save :: sname = 'json_file_load_from_string:'
         call this%err%werrfl2(class//sname//' started')

         call this%input%load_from_string(str)

         call this%err%werrfl2(class//sname//' ended')

      end subroutine load_from_string
!
      function found(this, path)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         logical :: found
! local data
         character(len=38), save :: sname = 'json_file_found:'

         call this%err%werrfl2(class//sname//' started')

         call this%input%info(path, found=found)

         call this%err%werrfl2(class//sname//' ended')

      end function found
!
      subroutine json_file_get_object(this, path, p)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         type(json_value),pointer,intent(out) :: p
! local data
         character(len=38), save :: sname = 'json_file_get_object:'
         character(len=:), allocatable :: error
         logical :: st

         call this%err%werrfl2(class//sname//' started')

         call this%input%get(path, p)
         if (this%input%failed()) then
            call this%input%check_for_errors(st, error)
            call this%input%clear_exceptions()
            call this%err%equit(class//sname//error)
         end if         

         call this%err%werrfl2(class//sname//' ended')

      end subroutine json_file_get_object
!
      subroutine json_file_get_integer(this, path, val)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         integer,intent(out) :: val 
! local data
         character(len=38), save :: sname = 'json_file_get_integer:'
         character(len=:), allocatable :: error
         logical :: st

         call this%err%werrfl2(class//sname//' started')

         call this%input%get(path, val)
         if (this%input%failed()) then
            call this%input%check_for_errors(st, error)
            call this%input%clear_exceptions()
            call this%err%equit(class//sname//error)
         end if         

         call this%err%werrfl2(class//sname//' ended')

      end subroutine json_file_get_integer
!
      subroutine json_file_get_double(this, path, val)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         real,intent(out) :: val 
! local data
         character(len=38), save :: sname = 'json_file_get_double:'
         character(len=:), allocatable :: error
         logical :: st

         call this%err%werrfl2(class//sname//' started')

         call this%input%get(path, val)
         if (this%input%failed()) then
            call this%input%check_for_errors(st, error)
            call this%input%clear_exceptions()
            call this%err%equit(class//sname//error)
         end if         

         call this%err%werrfl2(class//sname//' ended')

      end subroutine json_file_get_double
!
      subroutine json_file_get_logical(this, path, val)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         logical,intent(out) :: val 
! local data
         character(len=38), save :: sname = 'json_file_get_logical:'
         character(len=:), allocatable :: error
         logical :: st

         call this%err%werrfl2(class//sname//' started')

         call this%input%get(path, val)
         if (this%input%failed()) then
            call this%input%check_for_errors(st, error)
            call this%input%clear_exceptions()
            call this%err%equit(class//sname//error)
         end if         

         call this%err%werrfl2(class//sname//' ended')

      end subroutine json_file_get_logical
!
      subroutine json_file_get_string(this, path, val)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         character(len=:),allocatable,intent(out) :: val 
! local data
         character(len=38), save :: sname = 'json_file_get_string:'
         character(len=:), allocatable :: error
         logical :: st

         call this%err%werrfl2(class//sname//' started')

         call this%input%get(path, val)
         if (this%input%failed()) then
            call this%input%check_for_errors(st, error)
            call this%input%clear_exceptions()
            call this%err%equit(class//sname//error)
         end if         

         call this%err%werrfl2(class//sname//' ended')

      end subroutine json_file_get_string
!
      subroutine json_file_get_integer_vec(this, path, vec)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         integer,dimension(:),allocatable,intent(out) :: vec 
! local data
         character(len=38), save :: sname = 'json_file_get_integer_vec:'
         character(len=:), allocatable :: error
         logical :: st

         call this%err%werrfl2(class//sname//' started')

         call this%input%get(path, vec)
         if (this%input%failed()) then
            call this%input%check_for_errors(st, error)
            call this%input%clear_exceptions()
            call this%err%equit(class//sname//error)
         end if         

         call this%err%werrfl2(class//sname//' ended')

      end subroutine json_file_get_integer_vec
!
      subroutine json_file_get_double_vec(this, path, vec)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         real,dimension(:),allocatable,intent(out) :: vec 
! local data
         character(len=38), save :: sname = 'json_file_get_double_vec:'
         character(len=:), allocatable :: error
         logical :: st

         call this%err%werrfl2(class//sname//' started')

         call this%input%get(path, vec)
         if (this%input%failed()) then
            call this%input%check_for_errors(st, error)
            call this%input%clear_exceptions()
            call this%err%equit(class//sname//error)
         end if         

         call this%err%werrfl2(class//sname//' ended')

      end subroutine json_file_get_double_vec
!
      subroutine json_file_get_logical_vec(this, path, vec)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         logical,dimension(:),allocatable,intent(out) :: vec 
! local data
         character(len=38), save :: sname = 'json_file_get_logical_vec:'
         character(len=:), allocatable :: error
         logical :: st

         call this%err%werrfl2(class//sname//' started')

         call this%input%get(path, vec)
         if (this%input%failed()) then
            call this%input%check_for_errors(st, error)
            call this%input%clear_exceptions()
            call this%err%equit(class//sname//error)
         end if         

         call this%err%werrfl2(class//sname//' ended')

      end subroutine json_file_get_logical_vec
!
      subroutine json_file_get_string_vec(this, path, vec)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         character(len=*),dimension(:),allocatable,intent(out) :: vec 
! local data
         character(len=38), save :: sname = 'json_file_get_string_vec:'
         character(len=:), allocatable :: error
         logical :: st

         call this%err%werrfl2(class//sname//' started')

         call this%input%get(path, vec)
         if (this%input%failed()) then
            call this%input%check_for_errors(st, error)
            call this%input%clear_exceptions()
            call this%err%equit(class//sname//error)
         end if         

         call this%err%werrfl2(class//sname//' ended')

      end subroutine json_file_get_string_vec
!
      subroutine json_file_get_alloc_string_vec(this, path, vec, ilen)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         character(len=:),dimension(:),allocatable,intent(out) :: vec
         integer,dimension(:),allocatable,intent(out) :: ilen 
! local data
         character(len=38), save :: sname = 'json_file_get_alloc_string_vec:'
         character(len=:), allocatable :: error
         logical :: st

         call this%err%werrfl2(class//sname//' started')

         call this%input%get(path, vec, ilen)
         if (this%input%failed()) then
            call this%input%check_for_errors(st, error)
            call this%input%clear_exceptions()
            call this%err%equit(class//sname//error)
         end if         

         call this%err%werrfl2(class//sname//' ended')

      end subroutine json_file_get_alloc_string_vec
!
      subroutine json_file_get_root(this,p)

         implicit none

         class(input_json),intent(inout) :: this
         type(json_value),pointer,intent(out) :: p
! local data
         character(len=38), save :: sname = 'json_file_get_root:'
         call this%err%werrfl2(class//sname//' started')

         call this%input%get(p)

         call this%err%werrfl2(class//sname//' ended')

      end subroutine json_file_get_root                            
!
      subroutine  json_file_variable_info(this,path,n_children)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         integer,intent(out) :: n_children

! local data
         character(len=38), save :: sname = 'json_file_variable_info:'
         character(len=:), allocatable :: error
         logical :: st
         integer :: var_type
         character(len=:),allocatable :: name

         call this%err%werrfl2(class//sname//' started')

         call this%input%info(path, var_type=var_type, n_children=n_children, name=name)

         if (this%input%failed()) then
            call this%input%check_for_errors(st, error)
            call this%input%clear_exceptions()
            call this%err%equit(class//sname//error)
         end if         

         call this%err%werrfl2(class//sname//' ended')

      end subroutine json_file_variable_info
!
      end module input_class