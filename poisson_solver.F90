! Copyright 2018 and 2019, Shun Sakuraba
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! This file is suffixed f90 but requires Fortran03 features such as iso_c_binding.
! (f90 just means free-form source.)

! FFTW3, packed into an module
module fftw3_mpi
  use iso_c_binding
  include 'fftw3-mpi.f03'  
end module fftw3_mpi

! PFFT, packed into an module
module pfft
  use iso_c_binding
  use fftw3_mpi
  include 'pfft.f03'
end module pfft

module poisson_solver
  use pfft
  use iso_c_binding
  implicit none

  ! Public variables.

  ! rho is stored to this array and output as phi.
  ! Note this array may be larger than the size you specified in subroutine init(),
  ! to achieve optimal performance.
  real(8), allocatable :: phi_or_rho(:, :, :)
!DIR$ ATTRIBUTES ALIGN: 64:: phi_or_rho

  ! Public parameters.

  ! Boundary conditions. For the best performance, it is advised to use staggered conditions.
  ! * poisson_bc_dirichlet_zero is one of the Dirichlet boundary, where zeros lies outside the boundary.
  !   For example, for phi(1:nx, 1:ny, 1:nz), points at phi(0, :, :) and phi(nx + 1, :, :) are all zero.
  !   Same applies to both Y and Z.
  !   (Note this zero-padding is performed only in a conceptual manner and you do not have to allocate memory
  !    to these 0 or nx+1 regions.)
  ! * poisson_bc_dirichlet_zero_staggered is similar to the above,
  !   but it zeros out half-shifted points.
  !   For phi(1:nx, 1:ny, 1:nz), points at phi(1/2, :, :) and phi(nx + 1/2, :, :) are all zero, in a conceptual manner.
  !   poisson_bc_dirichlet_zero_staggered is typically faster than poisson_bc_dirichlet_zero.
  ! * poisson_bc_neumann_zero is the Neumann boundary (first-order differntiation constraint)
  !   with phi' = 0 at the border.
  !   For example, for phi(1:nx, 1:ny, 1:nz), gradients at phi(0, :, :) and phi(nx + 1, :, :) are all zero.
  !   Same applies to both Y and Z.
  !   (Note this zero-padding is performed only in a conceptual manner and you do not have to allocate memory
  !    to these 0 or nx+1 regions.)
  ! * poisson_bc_neumann_zero_staggered is the Neumann boundary with phi' = 0 at the border.
  !   Unlike poisson_bc_neumann_zero, gradients at phi(1/2, :, :) and phi(nx + 1/2, :, :) are all zero.
  integer, parameter :: poisson_bc_dirichlet_zero           = 1
  integer, parameter :: poisson_bc_dirichlet_zero_staggered = 2
  integer, parameter :: poisson_bc_neumann_zero             = 3
  integer, parameter :: poisson_bc_neumann_zero_staggered   = 4


  ! Background scheme.
  ! Spectral approximation (Approximate the derivative from Fourier series)
  integer, parameter :: eigenvalue_spectrum = 1
  ! Finite difference approximation (Approximate the derivative from 3-point finite difference)
  integer, parameter :: eigenvalue_finite_difference_2 = 2

  ! Estimate, Measure or Exhaustive?
  integer, parameter :: initialize_estimate = 1
  integer, parameter :: initialize_measure = 2
  integer, parameter :: initialize_exhaustive = 3


  ! ----------------------------------------------------------------
  ! Private variables.

  ! TODO: move these from module variable (global) to a single type struct?
  type(c_ptr), private :: pfft_forward_plan = c_null_ptr
  type(c_ptr), private :: pfft_backward_plan = c_null_ptr

  real(8), allocatable, private :: transposed_buffer(:)
  integer(8), private :: local_no(3), local_o_start(3)
  integer(8), private :: global_no(3)
  integer(8), private :: local_ni(3), local_i_start(3)
  real(8), private :: gridlen(3)
  integer, private :: bc, eigenvalue_type_
  real(8), allocatable, private :: eigenvals_x(:), eigenvals_y(:), eigenvals_z(:)
  integer, private :: comm



  integer, parameter, private :: ind_forward = 1, ind_backward = 2
  integer, parameter, private :: plan_list(2, 4) = & ! of dimension (fwd/bwd, bc)
       reshape([ PFFT_RODFT00, PFFT_RODFT00, &
                 PFFT_RODFT10, PFFT_RODFT01, &
                 PFFT_REDFT00, PFFT_REDFT00, &
                 PFFT_REDFT10, PFFT_REDFT01 ], [ 2, 4 ])
  
  real(8), parameter, private :: pi = 3.14159265358979323846264338327950d0
  integer, parameter, private :: stderr = 6
  logical, parameter, private :: debug = .false.
  logical, private :: transposed_, verbosem
  real(8), parameter, private :: dbl_min = 2.2250738585072013D-308

  ! save all module variables
  save

  private :: multiply_factors, check_mpi_error
contains

  ! Initialize the 3D parallel poisson solver
  subroutine init(&
       globalxyz, comm_cart_3d, &
       gridlenxyz, &
       myposition, &
       boundary_condition, &
       eigenvalue_type, &
       transposed, &
       verbose, &
       measure_type, &
       initialize_fftw, &
       initialize_pfft)
    use mpi
    use pfft
    !$ use omp_lib

    implicit none

    integer, intent(in) :: globalxyz(3)
    integer, intent(in) :: comm_cart_3d
    real(8), intent(in) :: gridlenxyz(3)
    integer, intent(in) :: myposition(3)
    integer, intent(in) :: boundary_condition
    integer, optional, intent(in) :: eigenvalue_type
    logical, optional, intent(in) :: transposed                         ! Default: .true. (faster)
    logical, optional, intent(in) :: verbose                            ! Default: .false.
    integer, optional, intent(in) :: measure_type                       ! Default: initialize_measure (medium)
    logical, optional, intent(in) :: initialize_fftw, initialize_pfft   ! Default: .true. (initialize by this routine)

    integer :: cart_dims(3)
    logical :: is_periodic(3)
    integer :: coords(3)
    integer :: rank, wrank
    integer :: ierr
    integer(c_int) :: forward_transpose_type, backward_transpose_type
    integer(c_int) :: forward_dft_type, backward_dft_type
    integer(8) :: output_size
    real(8) :: wtime0, wtime1
    integer(c_int) :: measure_flag
    integer :: i, j, cartcomm_size
    integer(8) :: li
    integer(8) :: global_no_with_pad(3), dftsize(3)
    real(8) :: worldsize(3), offsets(3), normalize_factor

    integer(8), allocatable :: local_i_start_print(:, :), local_o_start_print(:, :)
    integer(8), allocatable :: local_ni_print(:, :), local_no_print(:, :)

    logical :: initialize_fftw_, initialize_pfft_

    initialize_fftw_ = .true.
    if(present(initialize_fftw)) initialize_fftw_ = initialize_fftw
    initialize_pfft_ = .true.
    if(present(initialize_pfft)) initialize_pfft_ = initialize_pfft

#ifndef NO_USE_FFTW_THREADS
    if(initialize_fftw_) then
       !$ ierr = fftw_init_threads()
       !$ if(ierr == 0) stop "Hybrid execution requrested but failed to init fftw with threads"
    endif
#endif

    global_no(:) = globalxyz(3:1:-1) ! in an order of z/y/x, to match the C natural order
    gridlen(:) = gridlenxyz(3:1:-1)  ! ibid.
    bc = boundary_condition
    if(present(eigenvalue_type)) then
       eigenvalue_type_ = eigenvalue_type
    else
       eigenvalue_type_ = eigenvalue_spectrum
    end if
    transposed_ = .true.
    if(present(transposed)) transposed_ = transposed

    call MPI_Comm_rank(comm_cart_3d, rank, ierr)
    verbosem = .false.
    if(present(verbose)) then
       verbosem = verbose .and. rank == 0
    end if
    call check_mpi_error(ierr, "poisson_solver%init(MPI_Comm_rank:comm_cart_3d)")
    call MPI_Comm_rank(MPI_COMM_WORLD, wrank, ierr)
    call check_mpi_error(ierr, "poisson_solver%init(MPI_Comm_rank:world)")

    call MPI_Cart_Get(comm_cart_3d, 3, &
         cart_dims, is_periodic, coords, ierr)
    call check_mpi_error(ierr, "poisson_solver%init(MPI_Cart_Get)")
    if(verbosem) print *, "Psolver: Cartesian dims are: ", cart_dims
    comm = comm_cart_3d

    if(initialize_pfft_) then
       if(verbosem) print *, "Psolver: Initializing PFFT"
       call pfft_init() ! XXX: other routines may use pfft in the future
       if(verbosem) print *, "Psolver: Finished initializing PFFT"
    endif

#ifndef NO_USE_FFTW_THREADS
    if(initialize_fftw_) then
       !$ call fftw_plan_with_nthreads(omp_get_max_threads())
    endif
#endif

    ! yay, assertions always save you a day
    if(.not. all(is_periodic .eqv. .false.)) then
       write(stderr, *) "poisson_solver%init: comm_cart_3d is periodic"
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    if(.not. all(myposition(3:1:-1) == coords)) then
       write(stderr, *) "Rank ", rank, ": poisson_solver%init: comm_cart_3d and myposition mismatches" // &
            "(comm_cart 3d(z/y/x): ", coords(3:1:-1), " myposition: ", myposition, ")"
       write(stderr, *) "Note: this typically means that one of dimension is not divisable, " // &
            "or 3D to 2D redistribution may be difficult " // &
            "(for the distribution over NX/PX * NY/PY * NZ/PZ, there must be factors QY, QZ s.t. " // &
            "QYQZ = PX and NY/(PY QY) and NZ/(PZ QZ) is divisable)."
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    if(transposed_) then
       if(verbosem) print *, "Psolver: using transposed mode for PFFT (faster!)"
       forward_transpose_type  = PFFT_TRANSPOSED_OUT
       backward_transpose_type = PFFT_TRANSPOSED_IN
    else
       if(verbosem) print *, "Psolver: using non-transposed mode for PFFT (slower, for debug!)"
       forward_transpose_type  = PFFT_TRANSPOSED_NONE
       backward_transpose_type = PFFT_TRANSPOSED_NONE
    end if

    ! Initialize pfft
    if(verbosem) print *, "Psolver: Computing memory allocation per rank (this may take a few sec)"
    wtime0 = MPI_Wtime()
    output_size = pfft_local_size_r2r_3d( &
         global_no, comm_cart_3d, &
         forward_transpose_type + PFFT_DESTROY_INPUT, &
         local_ni, local_i_start, &
         local_no, local_o_start)
    wtime1 = MPI_Wtime()
    if(verbosem) print *, "Psolver: Rank / grid allocation took", wtime1 - wtime0, " sec"
    if(verbosem) print *, "Psolver: allocating FFT buffer, size=", output_size
    ! output_size MAY exceed product(local_ni)
    allocate(phi_or_rho(local_ni(3), local_ni(2),&
         (output_size + local_ni(2) * local_ni(3) - 1) / (local_ni(2) * local_ni(3))))
    if(verbosem) print *, "Psolver: allocated phi_or_rho with size", shape(phi_or_rho)
    allocate(transposed_buffer(0:(output_size - 1)))
    if(verbosem) print *, "Psolver: local input size (Z/Y/X)", int(local_ni, 4)
    if(verbosem) print *, "Psolver: local output size (may be transposed)", int(local_no, 4)

    if(present(verbose)) then
       if(verbose .and. debug) then
          call MPI_Comm_size(comm_cart_3d, cartcomm_size, ierr)
          call check_mpi_error(ierr, "poisson_solver%init(MPI_Comm_size)")
          allocate(local_i_start_print(3, cartcomm_size))
          allocate(local_o_start_print(3, cartcomm_size))
          allocate(local_ni_print(3, cartcomm_size))
          allocate(local_no_print(3, cartcomm_size))
          call MPI_Gather(local_i_start(1), 3, MPI_LONG_LONG, &
               local_i_start_print(1, 1), 3, MPI_LONG_LONG, &
               0, comm_cart_3d, ierr)
          call check_mpi_error(ierr, "poisson_solver%init(MPI_Gather)")
          call MPI_Gather(local_o_start(1), 3, MPI_LONG_LONG, &
               local_o_start_print(1, 1), 3, MPI_LONG_LONG, &
               0, comm_cart_3d, ierr)
          call check_mpi_error(ierr, "poisson_solver%init(MPI_Gather)")
          call MPI_Gather(local_ni(1), 3, MPI_LONG_LONG, &
               local_ni_print(1, 1), 3, MPI_LONG_LONG, &
               0, comm_cart_3d, ierr)
          call check_mpi_error(ierr, "poisson_solver%init(MPI_Gather)")
          call MPI_Gather(local_no(1), 3, MPI_LONG_LONG, &
               local_no_print(1, 1), 3, MPI_LONG_LONG, &
               0, comm_cart_3d, ierr)
          call check_mpi_error(ierr, "poisson_solver%init(MPI_Gather)")

          if(verbosem) then
             do i = 0, cartcomm_size - 1
                print "(A,I3,A,3I3,A,3I3)", "Psolver: In (", i, "):", &
                     (local_i_start_print(j, i + 1), j = 1,3), " / ", &
                     (local_ni_print(j, i + 1), j = 1,3)
             end do
             do i = 0, cartcomm_size - 1
                print "(A,I3,A,3I3,A,3I3)", "Psolver: Out(", i, "):", &
                     (local_o_start_print(j, i + 1), j = 1,3), " / ", &
                     (local_no_print(j, i + 1), j = 1,3)
             end do
          end if
          deallocate(local_i_start_print)
          deallocate(local_o_start_print)
          deallocate(local_ni_print)
          deallocate(local_no_print)
       end if
    end if

    if(.not. all(global_no / cart_dims == local_ni)) then
       write(stderr, *) "Global no / cart_dims do not match with local_ni: " // &
            "Former:", global_no / cart_dims, &
            "Latter:", local_ni, &
            "cart3d rank = ", rank, &
            "world rank = ", wrank
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    if(.not. all(local_i_start(3:1:-1) == local_ni(3:1:-1) * myposition)) then
       write(stderr, *) "local_i_start and local_ni * myposition mismatch_ni: " // &
            "Former:", local_i_start(3:1:-1), &
            "Latter:", local_ni(3:1:-1) * myposition, &
            "cart3d rank = ", rank, &
            "world rank = ", wrank
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    if(verbosem) print *, "Psolver: creating & optimizing PFFT plan"
    wtime0 = MPI_Wtime()
    measure_flag = PFFT_MEASURE
    if(present(measure_type)) then
       select case(measure_type)
       case(initialize_estimate)
          if(verbosem) print *, "Psolver: (Measure diabled)"
          measure_flag = PFFT_ESTIMATE
       case(initialize_measure)
          if(verbosem) print *, "Psolver: (Measure enabled)"
          measure_flag = PFFT_MEASURE
       case(initialize_exhaustive)
          if(verbosem) print *, "Psolver: (Measure enabled with exhaustive search)"
          measure_flag = PFFT_EXHAUSTIVE + PFFT_TUNE
       end select
    end if
    forward_dft_type = plan_list(ind_forward, bc)
    backward_dft_type = plan_list(ind_backward, bc)
    pfft_forward_plan =  pfft_plan_r2r_3d( &
         global_no, &
         phi_or_rho, transposed_buffer, &
         comm_cart_3d, &
         [forward_dft_type, forward_dft_type, forward_dft_type], &
         forward_transpose_type + PFFT_DESTROY_INPUT + measure_flag)
    pfft_backward_plan =  pfft_plan_r2r_3d( &
         global_no, &
         transposed_buffer, phi_or_rho, &
         comm_cart_3d, &
         [backward_dft_type, backward_dft_type, backward_dft_type], &
         backward_transpose_type + PFFT_DESTROY_INPUT + measure_flag)
    wtime1 =  MPI_Wtime()
    if(verbosem) print *, "Psolver: finished PFFT plan (", wtime1 - wtime0, "sec )"

    select case(boundary_condition)
    case(poisson_bc_dirichlet_zero)
       global_no_with_pad(:) = global_no(:) + 1
       offsets(:) = 1
    case(poisson_bc_dirichlet_zero_staggered)
       global_no_with_pad(:) = global_no(:)
       offsets(:) = 1
    case(poisson_bc_neumann_zero)
       global_no_with_pad(:) = global_no(:) - 1
       offsets(:) = 0
    case(poisson_bc_neumann_zero_staggered)
       global_no_with_pad(:) = global_no(:)
       offsets(:) = 0
    end select
    dftsize(:) = global_no_with_pad(:) * 2
    worldsize(:) = global_no_with_pad(:) * gridlen(:)
    normalize_factor = product(dble(dftsize))

    if(verbosem) print *, "Psolver: preparing eigenvalues"
    
    allocate(eigenvals_z(0:(local_no(1) - 1)))
    allocate(eigenvals_y(0:(local_no(2) - 1)))
    allocate(eigenvals_x(0:(local_no(3) - 1)))
    
    ! Uses spectral or finite difference approximation.
    select case(eigenvalue_type_)
    case(eigenvalue_spectrum)
       do li = 0, local_no(1) - 1
          eigenvals_z(li) = - (2 / gridlen(1) * (pi * (local_o_start(1) + li + offsets(1)) / dftsize(1))) ** 2
       end do
       do li = 0, local_no(2) - 1
          eigenvals_y(li) = - (2 / gridlen(2) * (pi * (local_o_start(2) + li + offsets(2)) / dftsize(2))) ** 2
       end do
       do li = 0, local_no(3) - 1
          eigenvals_x(li) = - (2 / gridlen(3) * (pi * (local_o_start(3) + li + offsets(3)) / dftsize(3))) ** 2
       end do
    case(eigenvalue_finite_difference_2)
       do li = 0, local_no(1) - 1
          eigenvals_z(li) = - (2 / gridlen(1) * sin(pi * (local_o_start(1) + li + offsets(1)) / dftsize(1))) ** 2
       end do
       do li = 0, local_no(2) - 1
          eigenvals_y(li) = - (2 / gridlen(2) * sin(pi * (local_o_start(2) + li + offsets(2)) / dftsize(2))) ** 2
       end do
       do li = 0, local_no(3) - 1
          eigenvals_x(li) = - (2 / gridlen(3) * sin(pi * (local_o_start(3) + li + offsets(3)) / dftsize(3))) ** 2
       end do
    end select
    eigenvals_z(:) = eigenvals_z(:) * normalize_factor
    eigenvals_y(:) = eigenvals_y(:) * normalize_factor
    eigenvals_x(:) = eigenvals_x(:) * normalize_factor
  end subroutine init

  subroutine finalize()
    use pfft
    implicit none

    deallocate(eigenvals_x)
    deallocate(eigenvals_y)
    deallocate(eigenvals_z)
    deallocate(transposed_buffer)
    deallocate(phi_or_rho)

    call pfft_destroy_plan(pfft_forward_plan)
    pfft_forward_plan = c_null_ptr

    call pfft_destroy_plan(pfft_backward_plan)
    pfft_backward_plan = c_null_ptr
  end subroutine finalize

  ! Run the actual Poisson solver. This routine is typically used to solve \nabla^2 P = -4 pi rho.
  ! Before calling this routine, phi_or_rho should be set to -4 pi rho. 
  subroutine run()
    use mpi
    implicit none
    ! unused actually
    integer :: ierr
    real(8) :: wtime1, wtime2

    if((.not. c_associated(pfft_forward_plan)) .or. (.not. c_associated(pfft_backward_plan))) then
       write(stderr, *) "poisson_solver%run: uninitialized call"
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if


    ! Perform 3D FFT
    if(verbosem) print *, "Psolver: performing forward transform"
    wtime1 = MPI_Wtime()
    call pfft_execute(pfft_forward_plan)
    wtime2 = MPI_Wtime()
    if(verbosem) print *, "Psolver: forward took", wtime2 - wtime1, " s"
    
    ! Calculate and multiply (optionally transposed) factors
    wtime1 = MPI_Wtime()
    call multiply_factors()
    wtime2 = MPI_Wtime()
    if(verbosem) print *, "Psolver: scale took", wtime2 - wtime1, " s"

    ! Perform Inverse 3D FFT
    if(verbosem) print *, "Psolver: performing backward transform"
    wtime1 = MPI_Wtime()
    call pfft_execute(pfft_backward_plan)
    wtime2 = MPI_Wtime()
    if(verbosem) print *, "Psolver: backward took", wtime2 - wtime1, " s"

    
  end subroutine run

  ! User must provide values at half grid from x, y, z ends (e.g., (xmin - 1/2, y, z) or (xmax + 1/2, y, z))
  subroutine run_with_dirichlet_boundary(xmplane, xpplane, ymplane, ypplane, zmplane, zpplane)
    use MPI
    implicit none

    ! X- / X+ plane, Y- / Y+ plane, Z- / Z+ plane.
    real(8), intent(in) :: xmplane(1:local_ni(2), 1:local_ni(1)), xpplane(1:local_ni(2), 1:local_ni(1))
    real(8), intent(in) :: ymplane(1:local_ni(3), 1:local_ni(1)), ypplane(1:local_ni(3), 1:local_ni(1))
    real(8), intent(in) :: zmplane(1:local_ni(3), 1:local_ni(2)), zpplane(1:local_ni(3), 1:local_ni(2))
    
    integer :: cart_dims(3), coords(3)
    logical :: is_periodic(3)
    integer :: ierr

    if(bc /= poisson_bc_dirichlet_zero_staggered) then
       write(stderr, *) "poisson_solver%run_with_dirichlet_boundary: must be initialized with poisson_bc_dirichlet_zero_staggered"
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    if(eigenvalue_type_ /= eigenvalue_finite_difference_2) then
       write(stderr, *) "poisson_solver%run_with_dirichlet_boundary: must be initialized with eigenvalue_finite_difference_2"
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    ! Boundary erasure: see for example Ricker, Astrophys. J. Supp. Ser., 176:293-300, 2008.
    call MPI_Cart_Get(comm, 3, &
         cart_dims, is_periodic, coords, ierr)
    call check_mpi_error(ierr, "poisson_solver%run_with_dirichlet_boundary(MPI_Cart_Get)")

    ! note that the coords and gridlen has the order of z/y/x
    ! X- / X+
    if(coords(3) == 0) then
       phi_or_rho(1, 1:local_ni(2), 1:local_ni(1)) = &
            phi_or_rho(1, 1:local_ni(2), 1:local_ni(1)) - &
            (2.0d0 / gridlen(3) ** 2) * xmplane(1:local_ni(2), 1:local_ni(1))
    endif
    if(coords(3) == cart_dims(3) - 1) then
       phi_or_rho(local_ni(3), 1:local_ni(2), 1:local_ni(1)) = &
            phi_or_rho(local_ni(3), 1:local_ni(2), 1:local_ni(1)) - &
            (2.0d0 / gridlen(3) ** 2) * xpplane(1:local_ni(2), 1:local_ni(1))
    endif
    ! Y- / Y+
    if(coords(2) == 0) then
       phi_or_rho(1:local_ni(3), 1, 1:local_ni(1)) = &
            phi_or_rho(1:local_ni(3), 1, 1:local_ni(1)) - &
            (2.0d0 / gridlen(2) ** 2) * ymplane(1:local_ni(3), 1:local_ni(1))
    endif
    if(coords(2) == cart_dims(2) - 1) then
       phi_or_rho(1:local_ni(3), local_ni(2), 1:local_ni(1)) = &
            phi_or_rho(1:local_ni(3), local_ni(2), 1:local_ni(1)) - &
            (2.0d0 / gridlen(2) ** 2) * ypplane(1:local_ni(3), 1:local_ni(1))
    endif
    ! Z- / Z+
    if(coords(1) == 0) then
       phi_or_rho(1:local_ni(3), 1:local_ni(2), 1) = &
            phi_or_rho(1:local_ni(3), 1:local_ni(2), 1) - &
            (2.0d0 / gridlen(1) ** 2) * zmplane(1:local_ni(3), 1:local_ni(2))
    endif
    if(coords(1) == cart_dims(1) - 1) then
       phi_or_rho(1:local_ni(3), 1:local_ni(2), local_ni(1)) = &
            phi_or_rho(1:local_ni(3), 1:local_ni(2), local_ni(1)) - &
            (2.0d0 / gridlen(1) ** 2) * zpplane(1:local_ni(3), 1:local_ni(2))
    endif

    call run()
  end subroutine run_with_dirichlet_boundary
  
  subroutine multiply_factors()
    implicit none 

    integer(8) :: i, j, k
    integer(8) :: pos
    real(8) :: lambda_1, lambda_2, lambda_3
    real(8) :: lambda
    integer, parameter :: stderr = 6

    if((bc == poisson_bc_neumann_zero .or. &
        bc == poisson_bc_neumann_zero_staggered) .and. &
       all(local_o_start(:) == 0) .and. &
       all(local_no(:) >= 1)) then
       if(abs(transposed_buffer(0)) > 1.D-6 * product(dble(global_no(:)))) then
          write (stderr, *) "Warning: Neumann boundary is specified but the average of rho is non-zero"
       end if
       transposed_buffer(0) = 0.0d0
    end if

    ! do this in C-order.
    if(transposed_) then
       ! transposed version
       !$omp parallel do collapse(2) default(none) &
       !$omp private(i, j, k, lambda_1, lambda_2, lambda_3, lambda, pos) &
       !$omp shared(local_no, eigenvals_x, eigenvals_y, eigenvals_z, transposed_buffer)
       do i = 0, local_no(2) - 1
          do j = 0, local_no(3) - 1
             lambda_2 = eigenvals_y(i)
             lambda_3 = eigenvals_x(j)
             ! prevent nan at (0, 0, 0)
             lambda = min(lambda_2 + lambda_3, -DBL_MIN) ! lambdas should be negative
             do k = 0, local_no(1) - 1
                lambda_1 = eigenvals_z(k)
                pos = i * local_no(3) * local_no(1) + j * local_no(1) + k
                transposed_buffer(pos) = &
                     transposed_buffer(pos) / (lambda_1 + lambda)
             end do
          end do
       end do
    else
       !$omp parallel do collapse(2) default(none) &
       !$omp private(i, j, k, lambda_1, lambda_2, lambda_3, lambda, pos) &
       !$omp shared(local_no, eigenvals_x, eigenvals_y, eigenvals_z, transposed_buffer)
       do i = 0, local_no(1) - 1
          do j = 0, local_no(2) - 1
             lambda_1 = eigenvals_z(i)
             lambda_2 = eigenvals_y(j)
             ! prevent nan at (0, 0, 0)
             lambda = min(lambda_1 + lambda_2, -DBL_MIN) ! lambdas should be negative
             do k = 0, local_no(3) - 1
                lambda_3 = eigenvals_x(k)
                pos = i * local_no(2) * local_no(3) + j * local_no(3) + k
                transposed_buffer(pos) = &
                     transposed_buffer(pos) / (lambda_3 + lambda)
             end do
          end do
       end do
    end if
  end subroutine multiply_factors

  subroutine check_mpi_error(ierr, wherestring)
    use mpi
    implicit none
    integer, intent(in) :: ierr
    character(len=*) :: wherestring
    integer :: resultlen
    character(len=MPI_MAX_ERROR_STRING) :: errstr
    integer :: wtf, wtf2
    integer :: rank

    if(ierr == MPI_SUCCESS) return
   
    call MPI_Error_string(ierr, errstr, resultlen, wtf)
    if(wtf /= MPI_SUCCESS) goto 999
    
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, wtf)
    if(wtf /= MPI_SUCCESS) goto 999

    write(stderr, *) "Rank(in world): ", rank, " error: ", errstr(1:resultlen), " at ", wherestring
    flush(stderr)
    call MPI_Abort(MPI_COMM_WORLD, ierr, wtf)
    stop

999 continue
    ! wtf
    write(stderr, *) "Failed to translate errcode, MPI_Error_string / MPI_Comm_rank returned ", wtf
    write(stderr, *) "Original error code = ", ierr, " position ", wherestring
    call MPI_Abort(MPI_COMM_WORLD, wtf, wtf2)
    stop
  end subroutine check_mpi_error
  
end module poisson_solver
  
