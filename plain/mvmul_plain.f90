program mvmul_plain
        use, intrinsic :: iso_fortran_env
        use :: hdf5
        implicit none
        integer, parameter :: dp = REAL64 ! double precision
        integer, parameter :: i32 = INT32 ! 32-bit integers

        ! matrix MxN, vector NxO=1 (easily extendable to matmul)
        integer(i32), parameter :: M=10000_i32
        integer(i32), parameter :: N=10000_i32
        integer(i32), parameter :: O=1_i32

        ! matrix and vectors
        real(dp), dimension(:,:), allocatable :: mat
        real(dp), dimension(:,:), allocatable :: vec
        real(dp), dimension(:,:), allocatable :: prod

        ! support for indexes
        integer(i32) :: i, j, k

        ! support for timing
        real(dp) :: startT, endT


        ! ALLOCATE 
        allocate(mat(M,N))
        allocate(vec(N,O))
        allocate(prod(M,O))


        ! FILL UP ARRAYS
        prod = 0.0_dp
        call random_seed()
        call random_number(mat)
        call random_number(vec)


        ! MAT-VEC MULTIPLICATION
        call cpu_time(startT)
        do k = 1, N
                do j = 1, O
                        do i = 1, M
                                prod(i,j) = prod(i,j) + mat(i,k)*vec(k,j)
                        enddo
                enddo
        enddo
        call cpu_time(endT)

        write(*,*) 'MAT-VEC multiplication took: ', (endT-startT), ' s'

        deallocate(mat)
        deallocate(vec)
        deallocate(prod)


end program mvmul_plain
