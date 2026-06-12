! Micro-calibration program for dropsonde_gdb.py against real gdb/gfortran.
! Mimics every structure the dumper depends on, with module names chosen so
! the constituent-capture code paths run verbatim:
!   - module constituents / cnst_name        (CAM-role capture)
!   - module cam_constituents / const_props  (SIMA-role capture)
!   - assumed-shape and explicit-shape dummies, non-square 2-D/3-D arrays
!   - a strided actual argument (array slice)
!   - a scheme called both from "toplevel" and nested inside another scheme
module constituents
  implicit none
  character(len=16) :: cnst_name(5)
end module constituents

module cam_constituents
  implicit none
  type props_t
    character(len=:), allocatable :: var_std_name
  end type props_t
  type prop_ptr_t
    type(props_t), pointer :: prop => null()
  end type prop_ptr_t
  type(prop_ptr_t), pointer :: const_props(:) => null()
  integer :: num_constituents = 0
end module cam_constituents

module cam_comp
  implicit none
contains
  subroutine cam_run1()
    ! timestep sentinel: dropsonde counts hits here
  end subroutine cam_run1
end module cam_comp

module phys
  implicit none
contains
  subroutine inner_run(ncol, pver, a, b, errflg)
    integer, intent(in) :: ncol, pver
    real(8), intent(in) :: a(:, :)
    real(8), intent(out) :: b(ncol, pver)
    integer, intent(out) :: errflg
    b(:, :) = a(:, :) * 2.0_8
    errflg = 0
  end subroutine inner_run

  subroutine outer_run(ncol, pver, q3, o, errflg)
    integer, intent(in) :: ncol, pver
    real(8), intent(in) :: q3(:, :, :)
    real(8), intent(out) :: o(:, :)
    integer, intent(out) :: errflg
    call inner_run(ncol, pver, q3(:, :, 2), o, errflg)
  end subroutine outer_run
end module phys

program cal
  use constituents
  use cam_constituents
  use cam_comp
  use phys
  implicit none
  integer :: i, errflg
  real(8) :: q3(4, 3, 2), o(4, 3), os(3, 3)
  cnst_name = [character(len=16) :: 'Q', 'CLDLIQ', 'CLDICE', 'NUMLIQ', &
               'NUMICE']
  allocate(const_props(3))
  num_constituents = 3
  do i = 1, 3
    allocate(const_props(i)%prop)
    const_props(i)%prop%var_std_name = 'std_name_'//achar(48 + i)
  end do
  q3 = reshape([(real(i, 8), i = 1, 24)], [4, 3, 2])
  call cam_run1()                              ! "timestep" 1
  call outer_run(4, 3, q3, o, errflg)          ! toplevel + nested inner
  call cam_run1()                              ! "timestep" 2
  call inner_run(3, 3, q3(1:3, :, 1), os, errflg)  ! strided actual arg
  print *, o(1, 1), os(1, 1)
end program cal
