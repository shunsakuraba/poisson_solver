module run_test_poisson
contains
  subroutine run_poisson(ndim, nprocxyz, &
       gridlen, boundary_condition, &
       rhofunction, &
       label, expectedresult, transposed)
    use mpi
    use poisson_solver, poisson_init => init
    implicit none

    interface
       function coord_to_real(ix, iy, iz, rx, ry, rz) result(res)
         integer, intent(in) :: ix, iy, iz
         real(8), intent(in) :: rx, ry, rz
         real(8) :: res
       end function coord_to_real
    end interface

    integer, intent(in) :: ndim(3), nprocxyz(3)
    real(8), intent(in) :: gridlen(3)
    integer, intent(in) :: boundary_condition
    procedure(coord_to_real) :: rhofunction
    character(len=*), intent(in) :: label
    procedure(coord_to_real), optional :: expectedresult
    logical, optional, intent(in) :: transposed


    integer :: mpi_comm_cart, ierr

    integer :: rank, nworldsize
    integer :: nprocxyz_corder(3)
    integer :: coords_corder(3)
    integer :: myposition(3)
    integer :: ix, iy, iz, igx, igy, igz
    real(8) :: rx, ry, rz
    integer :: dim_eachproc(3)

    real(8), allocatable :: buf_for_error_calc(:, :, :), buf_for_recv(:, :, :)
    integer :: irank, other_coords(3)
    real(8) :: mse, msesum
    integer :: status(MPI_STATUS_SIZE)
    real(8) :: rhosum

    logical :: use_transposed_derived

    nprocxyz_corder = nprocxyz(3:1:-1)

    call MPI_Cart_create(MPI_COMM_WORLD, 3, &
         nprocxyz_corder, &
         [.false., .false., .false.], &
         .true., &
         mpi_comm_cart, ierr)
    if(ierr /= MPI_SUCCESS) then
       stop "mpi_cart_create"
    endif

    call MPI_Comm_rank(mpi_comm_cart, rank, ierr)
    if(ierr /= MPI_SUCCESS) then
       stop "mpi_comm_rank"
    endif
    call MPI_Comm_size(mpi_comm_cart, nworldsize, ierr)
    if(ierr /= MPI_SUCCESS) then
       stop "mpi_comm_size"
    endif
    if(rank == 0) print *, "Running Poisson solver test ", label

    call MPI_Cart_coords(mpi_comm_cart, rank, 3, coords_corder, ierr)
    if(ierr /= MPI_SUCCESS) then
       stop "cmpi_cart_coords"
    endif
    myposition = coords_corder(3:1:-1)

    dim_eachproc = ndim / nprocxyz

    use_transposed_derived = .false.
    if(present(transposed)) then
       use_transposed_derived = transposed
    end if
    call poisson_init(ndim, mpi_comm_cart, &
         gridlen, &
         myposition, &
         boundary_condition, &
         verbose = .true., &
         transposed = use_transposed_derived, &
         measure = .false.)

    do iz = 0, dim_eachproc(3) - 1
       igz = iz + dim_eachproc(3) * myposition(3)
       rz = igz * gridlen(3)
       do iy = 0, dim_eachproc(2) - 1
          igy = iy + dim_eachproc(2) * myposition(2)
          ry = igy * gridlen(2)
          do ix = 0, dim_eachproc(1) - 1
             igx = ix + dim_eachproc(1) * myposition(1)
             rx = igx * gridlen(1)
             phi_or_rho(ix + 1, iy + 1, iz + 1) = rhofunction(igx, igy, igz, rx, ry, rz)
          end do
       end do
    end do

    call run()

    mse = 0.0d0
    if(present(expectedresult)) then
       do iz = 0, dim_eachproc(3) - 1
          igz = iz + dim_eachproc(3) * myposition(3)
          rz = igz * gridlen(3)
          do iy = 0, dim_eachproc(2) - 1
             igy = iy + dim_eachproc(2) * myposition(2)
             ry = igy * gridlen(2)
             do ix = 0, dim_eachproc(1) - 1
                igx = ix + dim_eachproc(1) * myposition(1)
                rx = igx * gridlen(1)
                mse = mse + (phi_or_rho(ix + 1, iy + 1, iz + 1) - expectedresult(ix, iy, iz, rx, ry, rz)) ** 2
             end do
          end do
       end do

       msesum = 0.0d0
       call MPI_Reduce(mse, msesum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            0, mpi_comm_cart, ierr)
       if(ierr /= MPI_SUCCESS) stop "Reduce"
       mse = msesum
    else
       allocate(buf_for_error_calc(0:ndim(1)-1, 0:ndim(2)-1, 0:ndim(3)-1))
       allocate(buf_for_recv(dim_eachproc(1), dim_eachproc(2), dim_eachproc(3)))
       buf_for_error_calc(:, :, :) = 0.0d0
       do irank = 0, nworldsize - 1
          if(irank == 0 .and. rank == 0) then
             buf_for_error_calc(&
                  0:dim_eachproc(1) - 1,&
                  0:dim_eachproc(2) - 1,&
                  0:dim_eachproc(3) - 1) = &
                  phi_or_rho(1:dim_eachproc(1), &
                             1:dim_eachproc(2), &
                             1:dim_eachproc(3))
          elseif(irank == rank) then
             ! send
             call MPI_Send(phi_or_rho, dim_eachproc(1) * dim_eachproc(2) * dim_eachproc(3), MPI_DOUBLE_PRECISION, &
                  0, 0, mpi_comm_cart, ierr)
             if(ierr /= MPI_SUCCESS) stop "MPI_Send"
          elseif(rank == 0) then
             ! recv
             call MPI_Recv(buf_for_recv, dim_eachproc(1) * dim_eachproc(2) * dim_eachproc(3), MPI_DOUBLE_PRECISION, &
                  irank, 0, mpi_comm_cart, status, ierr)
             if(ierr /= MPI_SUCCESS) stop "MPI_Recv"
             call MPI_Cart_coords(mpi_comm_cart, irank, 3, other_coords, ierr)
             ! note the order, other_coords are in C order
             buf_for_error_calc(other_coords(3) * dim_eachproc(1) : (other_coords(3) + 1) * dim_eachproc(1) - 1, &
                                other_coords(2) * dim_eachproc(2) : (other_coords(2) + 1) * dim_eachproc(2) - 1, &
                                other_coords(1) * dim_eachproc(3) : (other_coords(1) + 1) * dim_eachproc(3) - 1) &
                                = buf_for_recv(:, :, :)
          end if
       end do

       if(rank == 0) then
          do igz = 0, ndim(3) - 1
             rz = igz * gridlen(3)
             do igy = 0, ndim(2) - 1
                ry = igy * gridlen(2)
                do igx = 0, ndim(1) - 1
                   rx = igx * gridlen(1)
                   mse = mse + &
                        (rhofunction(igx, igy, igz, rx, ry, rz) - &
                        (ddx(igx, igy, igz) + ddy(igx, igy, igz) + ddz(igx, igy, igz))) ** 2
                end do
             end do
          end do
       end if

       if(.true. .and. rank == 0) then
          do igz = 0, min(1, ndim(3) - 1)
          rz = igz * gridlen(3)
          do igy = 0, ndim(2) - 1
             ry = igy * gridlen(2)
             do igx = 0, ndim(1) - 1
                rx = igx * gridlen(1)
                print *, igx, igy, igz, rhofunction(igx, igy, igz, rx, ry, rz), &
                     ddx(igx, igy, igz) + ddy(igx, igy, igz) + ddz(igx, igy, igz)
             end do
          end do
       end do
       end if

       deallocate(buf_for_error_calc)
       deallocate(buf_for_recv)
    end if
    mse = sqrt(mse / (ndim(1) * ndim(2) * ndim(3)))
    if(rank == 0) print *, "Poisson ", label , " MSE = ", mse 

    call MPI_Comm_free(mpi_comm_cart, ierr)
    if(ierr /= MPI_SUCCESS) then
       stop "mpi_comm_free"
    end if

    call finalize()
  contains
    function is_on_border(igx, igy, igz, ndim) result(r)
      implicit none
      integer, intent(in) :: igx, igy, igz, ndim(3)
      logical :: r
      ! I already regretted this. Twice.
      if((igx == 0           .and. igy == 0) .or. &
           (igx == 0           .and. igy == ndim(2) - 1) .or. &
           (igx == ndim(1) - 1 .and. igy == 0) .or. &
           (igx == ndim(1) - 1 .and. igy == ndim(2) - 1) .or. &
           (igy == 0           .and. igz == 0) .or. &
           (igy == 0           .and. igz == ndim(3) - 1) .or. &
           (igy == ndim(2) - 1 .and. igz == 0) .or. &
           (igy == ndim(2) - 1 .and. igz == ndim(3) - 1) .or. &
           (igz == 0           .and. igx == 0) .or. &
           (igz == 0           .and. igx == ndim(1) - 1) .or. &
           (igz == ndim(3) - 1 .and. igx == 0) .or. &
           (igz == ndim(3) - 1 .and. igx == ndim(1) - 1)) then
         r = .true.
      else
         r = .false.
      end if
    end function is_on_border

    function ddx(igx, igy, igz) result(r)
      implicit none
      integer, intent(in) :: igx, igy, igz
      real(8) :: r
      r = ddn(1, [igx, igy, igz])
    end function ddx

    function ddy(igx, igy, igz) result(r)
      implicit none
      integer, intent(in) :: igx, igy, igz
      real(8) :: r
      r = ddn(2, [igx, igy, igz])
    end function ddy

    function ddz(igx, igy, igz) result(r)
      implicit none
      integer, intent(in) :: igx, igy, igz
      real(8) :: r
      r = ddn(3, [igx, igy, igz])
    end function ddz

    function ddn(ix, coords) result(r)
      integer, intent(in) :: ix, coords(3)
      real(8) :: r

      integer :: coordsm1(3), coordsp1(3)
      real(8) :: vm1, v, vp1, dx

      coordsm1(:) = coords(:)
      coordsm1(ix) = coordsm1(ix) - 1
      coordsp1(:) = coords(:)
      coordsp1(ix) = coordsp1(ix) + 1
      dx = gridlen(ix)

      if(coords(ix) == 0) then
         vm1 = 0.0d0
      else
         vm1 = buf_for_error_calc(coordsm1(1), coordsm1(2), coordsm1(3))
      end if
      v = buf_for_error_calc(coords(1), coords(2), coords(3))
      if(coords(ix) == ndim(ix) - 1) then
         vp1 = 0.0d0
      else
         vp1 = buf_for_error_calc(coordsp1(1), coordsp1(2), coordsp1(3))
      end if
      ! overwrite for boundary conditions
      if(coords(ix) == 0) then
         select case(boundary_condition)
         case(poisson_bc_dirichlet_zero)
            vm1 = 0.0d0
         case(poisson_bc_dirichlet_zero_staggered)
            vm1 = -v
         case(poisson_bc_neumann_zero)
            vm1 = vp1
         case(poisson_bc_neumann_zero_staggered)
            vm1 = v
         end select
      end if

      if(coords(ix) == ndim(ix) - 1) then
         select case(boundary_condition)
         case(poisson_bc_dirichlet_zero)
            vp1 = 0.0d0
         case(poisson_bc_dirichlet_zero_staggered)
            vp1 = -v
         case(poisson_bc_neumann_zero)
            vp1 = vm1
         case(poisson_bc_neumann_zero_staggered)
            vp1 = v
         end select
      end if

      ! Gotcha, typical three-point finite difference
      r = (vm1 * 1.0d0 - v * 2.0d0 + vp1 * 1.0d0) / (dx**2)
    end function ddn
  end subroutine run_poisson
end module run_test_poisson

function zero(ix, iy, iz, rx, ry, rz) result(res)
  implicit none
  integer, intent(in) :: ix, iy, iz
  real(8), intent(in) :: rx, ry, rz
  real(8) :: res
  
  res = 0.0d0
end function zero

! From PoisFFT
function sin_distr_dir(ix, iy, iz, rx, ry, rz) result(res)
  implicit none
  integer, intent(in) :: ix, iy, iz
  real(8), intent(in) :: rx, ry, rz
  real(8) :: res

  real(8), parameter :: d(3) = [1.0d0, 2.0d0, 3.0d0]
  integer, parameter :: world(3) = [10, 6, 5]
  real(8), parameter :: pi = 3.1415926535

  real(8) :: ls(3), x, y, z

  ls = d * (world + 1)
  x = (ix + 1) * d(1)
  y = (iy + 1) * d(2)
  z = (iz + 1) * d(3)
  
  res = &
       sin(3.0 * pi * (ix + 1) * d(1) / ls(1)) * &
       sin(0.5 * pi * (iy + 1) * d(2) / ls(2)) * & 
       sin(0.7 * pi * (iz + 1) * d(3) / ls(3))
end function sin_distr_dir

! From PoisFFT
function sin_simm_dir(ix, iy, iz, rx, ry, rz) result(res)
  implicit none
  integer, intent(in) :: ix, iy, iz
  real(8), intent(in) :: rx, ry, rz
  real(8) :: res

  real(8), parameter :: d(3) = [1.0d0, 1.0d0, 1.0d0]
  integer, parameter :: world(3) = [40, 40, 40]
  real(8), parameter :: pi = 3.1415926535

  real(8) :: ls(3), x, y, z

  ls = d * (world + 1)
  x = (ix + 1) * d(1)
  y = (iy + 1) * d(2)
  z = (iz + 1) * d(3)
  
  res = &
       sin(0.5 * pi * (ix + 1) * d(1) / ls(1)) * &
       sin(0.5 * pi * (iy + 1) * d(2) / ls(2)) * & 
       sin(0.5 * pi * (iz + 1) * d(3) / ls(3))
end function sin_simm_dir

function dir_results(ix, iy, iz, rx, ry, rz) result(res)
  implicit none
  integer, intent(in) :: ix, iy, iz
  real(8), intent(in) :: rx, ry, rz
  real(8) :: res

  real(8), parameter :: d(3) = [1.0d0, 2.0d0, 3.0d0]
  integer, parameter :: world(3) = [10, 6, 5]
  real(8), parameter :: pi = 3.1415926535

  real(8) :: ls(3)
  ls = d * (world + 1)


  res = -1.0d0 / (&
       (3 * pi) ** 2 * 2)
       
end function dir_results

program test_poisson_solver
  use mpi
  use poisson_solver
  use run_test_poisson
  implicit none
  integer :: ierr
  interface
     function zero(ix, iy, iz, rx, ry, rz) result(res)
       integer, intent(in) :: ix, iy, iz
       real(8), intent(in) :: rx, ry, rz
       real(8) :: res
     end function zero

     function sin_distr_dir(ix, iy, iz, rx, ry, rz) result(res)
       integer, intent(in) :: ix, iy, iz
       real(8), intent(in) :: rx, ry, rz
       real(8) :: res
     end function sin_distr_dir

     function sin_simm_dir(ix, iy, iz, rx, ry, rz) result(res)
       integer, intent(in) :: ix, iy, iz
       real(8), intent(in) :: rx, ry, rz
       real(8) :: res
     end function sin_simm_dir
  end interface

  call MPI_Init(ierr)
  ! 10/2 x 6/3 x 5/5 transform (in Fortran order)
  ! ==> 5/5 x 6/3 x 10/2 transform (in C order)
  ! ==> P2 = 2, Q0 = 1, Q1 = 2 (which is an only solution)
  ! ==> 5/(5*Q0) x 6/(3*Q1) x 10 transform in effect.
  ! ==> 6/5 x 10/6 x 5 transposed transform (remainder is rounded up, requiring 20 elements / process),
  !  or 5/5 x 6/3 x 10/2 untransposed (because FFT'ed input is re-transposed)
  !call run_poisson([10, 6, 5], [2, 3, 5], [1.0d0, 1.0d0, 1.0d0], &
  ! poisson_bc_dirichlet_zero, zero, "zero", expectedresult = zero)
  !call run_poisson([20, 1, 1], [1, 1, 1], [1.0d0, 1.0d0, 1.0d0], &
  !poisson_bc_dirichlet_zero, sin_simm_dir, "sin_simm_dir")

  !call run_poisson([4, 4, 4], [2, 2, 2], [1.0d0, 1.0d0, 1.0d0], &
  !poisson_bc_dirichlet_zero, sin_simm_dir, "sin_simm_dir")

  !call run_poisson([4, 4, 4], [2, 2, 2], [1.0d0, 1.0d0, 1.0d0], &
  !poisson_bc_dirichlet_zero, sin_distr_dir, "sin")

  call run_poisson([10, 6, 5], [2, 3, 5], [1.0d0, 1.0d0, 1.0d0], &
       poisson_bc_dirichlet_zero, sin_distr_dir, "sin")

  call run_poisson([10, 6, 5], [2, 3, 5], [1.0d0, 2.0d0, 3.0d0], &
       poisson_bc_dirichlet_zero, sin_distr_dir, "sin with unequal length")

  call run_poisson([10, 6, 5], [2, 3, 5], [1.0d0, 2.0d0, 3.0d0], &
       poisson_bc_dirichlet_zero, sin_distr_dir, "sin with unequal length transposed", transposed = .true.)

  call run_poisson([10, 6, 5], [2, 3, 5], [1.0d0, 2.0d0, 3.0d0], &
       poisson_bc_dirichlet_zero_staggered, sin_distr_dir, "sin with unequal length transposed, DS", transposed = .true.)

  call run_poisson([10, 6, 5], [2, 3, 5], [1.0d0, 2.0d0, 3.0d0], &
       poisson_bc_neumann_zero, sin_distr_dir, "sin with unequal length transposed, N", transposed = .true.)
  call run_poisson([10, 6, 5], [2, 3, 5], [1.0d0, 2.0d0, 3.0d0], &
       poisson_bc_neumann_zero_staggered, sin_distr_dir, "sin with unequal length transposed, NS/UT", transposed = .false.)
  call run_poisson([10, 6, 5], [2, 3, 5], [1.0d0, 2.0d0, 3.0d0], &
       poisson_bc_neumann_zero_staggered, sin_distr_dir, "sin with unequal length transposed, NS/T", transposed = .true.)
  if(.false.) then
  call run_poisson([5, 1, 1], [1, 1, 1], [1.0d0, 1.0d0, 1.0d0], &
       poisson_bc_neumann_zero_staggered, sin_distr_dir, "sin with unequal length transposed, NS/T", transposed = .true.)
end if

  call MPI_Finalize(ierr)
end program test_poisson_solver



  
