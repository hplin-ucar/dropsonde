! cal_mpi.f90 -- synthetic MPI mini-model for the MPMD one-rank calibration
! (cal_mpi_run.sh drives it through the real ./dropsonde driver).
!
! Shape mirrors a real model run: a once-per-timestep sentinel
! (step_begin), a compared science subroutine with intent-declared dummies
! (flux_calc), rank-local state derived deterministically from
! (rank, step, index), and a per-step MPI_Allreduce so the ranks are
! genuinely coupled -- which makes "gdb killed rank 0 and the whole job
! exited" a real assertion rather than a vacuous one.
!
! A planted divergence is read from perturb.txt in the run directory
! ("<in|out> <rank> <step> <index>"); absent file = clean run. "out"
! perturbs b INSIDE flux_calc (a bug in the subroutine: inputs match,
! outputs differ); "in" perturbs a upstream of the call (inputs differ).
module phys_mpi
  implicit none
  integer :: cur_step = 0
  integer :: my_rank = -1
  character(len=3) :: pkind = ""
  integer :: prank = -1, pstep = -1, pidx = -1
contains
  subroutine step_begin(t)
    integer, intent(in) :: t
    cur_step = t
  end subroutine step_begin

  subroutine read_perturb()
    integer :: u, ios
    open(newunit=u, file="perturb.txt", status="old", action="read", &
         iostat=ios)
    if (ios /= 0) return
    read(u, *, iostat=ios) pkind, prank, pstep, pidx
    if (ios /= 0) pkind = ""
    close(u)
  end subroutine read_perturb

  subroutine flux_calc(a, b, scale, ncol)
    integer, intent(in) :: ncol
    real(8), intent(in) :: a(ncol)
    real(8), intent(out) :: b(ncol)
    real(8), intent(in) :: scale
    integer :: i
    do i = 1, ncol
      b(i) = 2.0d0 * a(i) + scale
    end do
    ! the planted "bug": an output-side divergence inside the subroutine
    if (pkind == "out" .and. my_rank == prank .and. cur_step == pstep) then
      b(pidx) = b(pidx) + 1.0d0
    end if
  end subroutine flux_calc
end module phys_mpi

program cal_mpi
  use mpi
  use phys_mpi
  implicit none
  integer, parameter :: ncol = 8, nsteps = 3
  integer :: ierr, nranks, t, i
  real(8) :: a(ncol), b(ncol), mysum, total
  call read_perturb()
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
  do t = 1, nsteps
    call step_begin(t)
    do i = 1, ncol
      a(i) = dble(my_rank * 1000 + t * 100 + i)
    end do
    ! input-side planted divergence: upstream of the compared subroutine
    if (pkind == "in" .and. my_rank == prank .and. t == pstep) then
      a(pidx) = a(pidx) + 1.0d0
    end if
    call flux_calc(a, b, dble(t), ncol)
    mysum = sum(b)
    call MPI_Allreduce(mysum, total, 1, MPI_REAL8, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    if (my_rank == 0) print '(a,i0,a,f14.2)', 'step ', t, ' total ', total
  end do
  call MPI_Finalize(ierr)
end program cal_mpi
