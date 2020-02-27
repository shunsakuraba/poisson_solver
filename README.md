# Parallel Poisson solver

This library provides [a parallel Poission solver][VmolFFT] implemented in [Vmol][Vmol]. A notable feature is that the non-zero Dirichlet boundary condition can be used (solved using [Ricker's boundary erasure][Ricker]).

## Prerequisites

This library requires [FFTW](https://www.fftw.org/) and [PFFT](https://github.com/mpip/pfft).

## Usage

First, use `MPI_Cart_create` to generate the 3-dimensional mapping of MPI ranks to the 3-D cartesian coordinate.
```F90
   use MPI
   use poission_solver
   ...

   call MPI_Cart_create(MPI_COMM_WORLD, 3, &      ! 3 dimensions
         nprocxyz_corder, &                       ! Desired # of processors at each dimension, in a reversed order. Typically (PZ, PY, PX). PX*PY*PZ must be equal to total ranks.
         [.false., .false., .false.], &           ! Periodic = .false.
         .true., &                                ! Allow reordering (.true. for better performance, .false. may very slightly slow down but MPI_COMM_WORLD's node order matches the mappping)
         mpi_comm_cart, ierr)                     ! output communicator
	call MPI_Comm_rank(mpi_comm_cart, rank, ierr) ! When reorder = .true., this rank might not be equal to MPI_COMM_WORLD's.
    call MPI_Cart_coords(mpi_comm_cart, rank, 3, coords_corder, ierr)
    myposition = coords_corder(3:1:-1)            ! (px, py, pz) of my rank, 0 <= px < PX and such like.
```

Then, initialize the library by calling poisson_init.

```F90
    call poisson_init(ndim, mpi_comm_cart, &      ! ndim = total system size (NX, NY, NZ). (NX, NY, NZ) must be a multiple of (PX, PY, PZ).
         gridlen, &                               ! grid length (LX, LY, LZ)
         myposition, &                            ! my rank, only used to do sanity check
		 eigenvalue_type, &                       ! Use "eigenvalue_finite_difference_2" to use non-zero Dirichlet BC
         boundary_condition, &                    ! Use "poisson_bc_dirichlet_zero_staggered" to use non-zero Dirichlet BC
         verbose = .true., &                      ! Be loud and noisy
         measure = initialize_measure)            ! Use initial benchmark to choose faster routine (in PFFT and FFTW). Other choices are "initialize_estimate" (fast init, slow runtime) and "initialize_exhaustive" (slow init, fast runtime).
```

After the initialization is complete, the buffer `phi_or_rho` within poisson_solver library is allocated. This variable is a 3-dimensional matrix of `real(8)`. Set this variable to the RHS of the Poisson equation, typically -4 * pi * rho (see also notes below). Then, if you want to use non-zero Dirichlet boundary condition, run the solver by

```F90
    call run_with_dirichlet_boundary(xmplane, xpplane, ymplane, ypplane, zmplane, zpplane) ! xmplane is 2-dimensional (y,z) plane at (xmin - 1/2, y, z), of the callee's rank.
```

Or, if you want the boundary condition with zeros (Dirichlet, Neumann etc.), just do `call run()`. 

After executing `run_with_dirichlet_boundary` or `run`, `phi_or_rho` is set to the electrostatic potential. If you want to recalculate with different charge densities, just set `phi_or_rho` and call `run_with_dirichlet_boundary` or `run` again.

If you want to release all the resources allocated by this library, call `finalize`.

See also an example in `test_poisson_solver.f90`.

## Compile & link
The typical compilation commandline is:
```sh
mpif90 -O -c -I/path/to/FFTW/include -I/path/to/pfft/include -std=f2003 poission_solver.F90 -o poisson_solver.o
```

The typical linking order is:
```sh
mpif90 -o executable_file the_object_files_using_solver.o poisson_solver.o -L/path/to/pfft/lib -lpfft -L/path/to/fftw/lib -lfftw3_mpi -lfftw3
```

And please add as many optimization options you desire.

## Notes

### Size of `phi_or_rho`

The size of `phi_or_rho` might be larger than the size of (NX/PX, NY/PY, NZ/PZ). This is a feature of PFFT, required to achieve the better performance. Thus, for example,
```F90
    phi_or_rho(:, :, :) = -4 * pi * your_rho_variable(:, :, :)
```
*may fail due to the size difference.* You need to explicitly specify the size you want to save:
```F90
    CX = NX / PX
	CY = NY / PY
	CZ = NZ / PZ
    phi_or_rho(1:CX, 1:CY, 1:CZ) = -4 * pi * your_rho_variable(:, :, :)
```

### Fortran 90 or 2003?

The file name is suffixed with `.F90` but we use Fortran 2003. This is because Intel Fortran does not recognize `.F03` nor `.f03` as a suffix.

### Asking to cite our works
If you use this library in your academic work, we humbly ask you to cite [our work][VmolFFT]:
```
H. Takahashi, S. Sakuraba & A. Morita.
Large-Scale Parallel Implementation of Hartree-Fock Exchange Energy on the Real-Space Grids Using 3D-Parallel FFT.
Journal of Chemical Information and Modeling, in press. DOI: https://doi.org/10.1021/acs.jcim.9b01063 .
```

### See also

* Kudos to [PoisFFT](https://github.com/LadaF/PoisFFT), a library for solving the Poisson equation in parallel.
* [BigDFT](http://bigdft.org/Wiki/index.php?title=BigDFT_website) contains a Poisson solver library `Psolver`.

### TODOs
* Move ugly module global variables into a single type.

[Vmol]:  https://doi.org/10.1002/jcc.1082 "Takahashi, H.; Hori, T.; Hashimoto, H.; Nitta, T. A Hybrid QM/MM Method Employing Real Space Grids for QM Water in the TIP4P Water Solvents. J. Comput. Chem. 2001, 22, 1252--1261.; Takahashi, H.; Suzuoka, D.; Sakuraba, S.; Morita, A. Role of the Photosystem II as an Environment in the Oxidation Free Energy of the Mn Cluster from S1 to S2. J. Phys. Chem. B 2019, 123, 7081-7091."
[VmolFFT]: https://doi.org/10.1021/acs.jcim.9b01063 "Takahashi, H.; Sakuraba S.; Morita A. Large-Scale Parallel Implementation of Hartree-Fock Exchange Energy on the Real-Space Grids Using 3D-Parallel FFT. J. Chem. Inform. Model., in press. (DOI: https://doi.org/10.1021/acs.jcim.9b01063)."
[Ricker]: https://doi.org/10.1086/526425 "Ricker, P.M. A Direct Multigrid Poisson Solver for Oct-tree Adaptive Meshes. Astrophys. J. Suppl. Ser. 2008, 176, 293-300."

