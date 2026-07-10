! Micro-calibration for dropsonde's derived-type expansion (--capture)
! against real gdb + a Fortran compiler, shaped like the SE dycore
! (cal.f90 covers the CCPP-scheme paths; this file covers the pseudo-SDF /
! SE dycore paths -- keep them separate):
!   - element_t-like array-of-struct dummy with fixed-shape components
!     (CAM style), an allocatable component (CAM-SIMA style), nested
!     state/derived subtypes, an integer scalar and a character component
!   - scalar struct dummies: TimeLevel_t-like (integers, rotated between
!     timesteps) and hvcoord_t-like (real array + scalar)
!   - an unassociated pointer array-of-struct dummy (fvm(:) with CSLAM off)
!   - NO errmsg=''/errflg=0 idiom: entry capture must fall back to the
!     first executable statement
!   - a nested call (compute_like from inside prim_step_like) and two
!     "subcycles" per timestep, mirroring the dycore call pattern
module constituents
  implicit none
  character(len=16) :: cnst_name(3)
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
  integer :: nstep_count = 0
contains
  subroutine cam_run1()
    ! timestep sentinel: dropsonde counts hits here (side effect keeps
    ! the call alive under optimization)
    nstep_count = nstep_count + 1
  end subroutine cam_run1
end module cam_comp

module se_element_mod
  implicit none
  integer, parameter :: np = 2, nlev = 3, qsize = 2, timelevels = 2

  type state_t
    real(8) :: v(np, np, 2, nlev, timelevels)   ! fixed shape (CAM style)
    real(8) :: t(np, np, nlev, timelevels)
    real(8), allocatable :: qdp(:, :, :, :, :)  ! per-element heap (SIMA style)
  end type state_t

  type derived_t
    real(8) :: ft(np, np, nlev)
  end type derived_t

  type element_t
    integer :: localid
    type(state_t) :: state
    type(derived_t) :: derived
    real(8) :: spheremp(np, np)
    character(len=8) :: tag                     ! must be skipped gracefully
  end type element_t

  type timelevel_t
    integer :: nm1, n0, np1, nstep
  end type timelevel_t

  type hvcoord_t
    real(8) :: hyai(nlev + 1)
    real(8) :: ps0
  end type hvcoord_t
end module se_element_mod

module se_dycore
  use se_element_mod
  implicit none
contains
  subroutine compute_like(elem, nets, nete, np1, dt2, hvcoord)
    type(element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete, np1
    real(8), intent(in) :: dt2
    type(hvcoord_t), intent(in) :: hvcoord
    integer :: ie
    do ie = nets, nete
      elem(ie)%state%v(:, :, :, :, np1) = elem(ie)%state%v(:, :, :, :, np1) &
           + dt2 * hvcoord%ps0
      elem(ie)%derived%ft = elem(ie)%derived%ft + hvcoord%hyai(1)
    end do
  end subroutine compute_like

  subroutine prim_step_like(elem, fvm, nets, nete, dt, tl, hvcoord)
    type(element_t), intent(inout) :: elem(:)
    type(element_t), pointer :: fvm(:)
    integer, intent(in) :: nets, nete
    real(8), intent(in) :: dt
    type(timelevel_t), intent(in) :: tl
    type(hvcoord_t), intent(in) :: hvcoord
    integer :: ie
    do ie = nets, nete
      elem(ie)%state%t(:, :, :, tl%np1) = elem(ie)%state%t(:, :, :, tl%n0) + dt
      elem(ie)%state%qdp(:, :, :, :, 2) = elem(ie)%state%qdp(:, :, :, :, 1) &
           * 2.0_8
      elem(ie)%localid = elem(ie)%localid + 1
    end do
    call compute_like(elem, nets, nete, tl%np1, dt * 0.5_8, hvcoord)
  end subroutine prim_step_like
end module se_dycore

program cal_se
  use constituents
  use cam_constituents
  use cam_comp
  use se_element_mod
  use se_dycore
  implicit none
  integer, parameter :: nelem = 4
  type(element_t) :: elem(nelem)
  type(element_t), pointer :: fvm_null(:) => null()
  type(timelevel_t) :: tl
  type(hvcoord_t) :: hv
  integer :: ie, i, k

  cnst_name = [character(len=16) :: 'Q', 'CLDLIQ', 'RAINQM']
  allocate(const_props(2))
  num_constituents = 2
  do i = 1, 2
    allocate(const_props(i)%prop)
    const_props(i)%prop%var_std_name = 'std_name_' // achar(48 + i)
  end do

  do ie = 1, nelem
    allocate(elem(ie)%state%qdp(np, np, nlev, qsize, 2))
    elem(ie)%localid = ie
    elem(ie)%tag = 'elem' // achar(48 + ie)
    elem(ie)%spheremp = 100.0_8 * ie + 1.0_8
    elem(ie)%derived%ft = 100.0_8 * ie + 2.0_8
    do k = 1, timelevels
      elem(ie)%state%v(:, :, :, :, k) = 100.0_8 * ie + 10.0_8 * k
      elem(ie)%state%t(:, :, :, k) = 100.0_8 * ie + 20.0_8 * k
    end do
    elem(ie)%state%qdp(:, :, :, :, 1) = 100.0_8 * ie + 30.0_8
    elem(ie)%state%qdp(:, :, :, :, 2) = 0.0_8
    ! interior bumps so a stride/dim-order mistake shifts values instead
    ! of reproducing them
    elem(ie)%state%t(1, 2, 3, 1) = elem(ie)%state%t(1, 2, 3, 1) + 0.5_8
    elem(ie)%state%v(2, 1, 2, 3, 1) = elem(ie)%state%v(2, 1, 2, 3, 1) &
         + 0.25_8
    elem(ie)%state%qdp(2, 1, 3, 2, 1) = elem(ie)%state%qdp(2, 1, 3, 2, 1) &
         + 0.125_8
  end do

  tl%nm1 = 1; tl%n0 = 1; tl%np1 = 2; tl%nstep = 0
  hv%hyai = [(0.1_8 * i, i = 1, nlev + 1)]
  hv%ps0 = 1000.0_8

  call cam_run1()                                             ! timestep 1
  call prim_step_like(elem, fvm_null, 1, nelem, 0.5_8, tl, hv) ! subcycle 1
  call prim_step_like(elem, fvm_null, 1, nelem, 0.5_8, tl, hv) ! subcycle 2
  call cam_run1()                                             ! timestep 2
  tl%n0 = 2; tl%np1 = 1; tl%nstep = 1                          ! rotate levels
  call prim_step_like(elem, fvm_null, 1, nelem, 0.5_8, tl, hv)
  print *, elem(1)%state%t(1, 1, 1, 1), elem(nelem)%state%qdp(1, 1, 1, 1, 2)
end program cal_se
