! This program demonstrates hybrid programming.

      program madhybrid
      use omp_lib
      implicit none

      include 'mpif.h'

      integer :: rank, nproc, ierr, i, id
      integer, dimension(MPI_STATUS_SIZE) :: status

! Initialise MPI.

      call MPI_INIT(ierr)

! Get the topology of the communicator.
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

!$OMP PARALLEL PRIVATE(id)
      id = OMP_GET_THREAD_NUM()

      print *, 'I am process ', rank, ', thread ', id

!$OMP END PARALLEL

! Shut down MPI.
      call MPI_FINALIZE(ierr)

      end program madhybrid
