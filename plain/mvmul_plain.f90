program mvmul_plain
        use, intrinsic :: iso_fortran_env
        use :: hdf5
        implicit none
        integer, parameter :: dp = REAL64 ! double precision
        integer, parameter :: i32 = INT32 ! 32-bit integers

        ! matrix MxN, vector NxO=1 (easily extendable to matmul)
        integer(i32), parameter :: M=100_i32
        integer(i32), parameter :: N=100_i32
        integer(i32), parameter :: O=1_i32

        ! matrix and vectors
        real(dp), dimension(:,:), allocatable :: mat
        real(dp), dimension(:,:), allocatable :: vec
        real(dp), dimension(:,:), allocatable :: prod

        ! support for indexes
        integer(i32) :: i, j, k

        ! support for timing
        real(dp) :: startT, endT

        ! HDF5 parameters
        integer(hid_t) :: file_id       ! File identifier
        integer(hid_t) :: dset_id       ! Dataset identifier
        integer(hid_t) :: dspace_id     ! Dataspace idenfitier
        integer(hid_t) :: mspace_id     ! Memory dataspace identifier
        integer(hid_t) :: crp_list      ! Dataset creation property identifier
        integer(i32) :: error           ! Error checks
        integer(i32) :: space_rank      ! Number of dimensions in the data space
        integer(hsize_t) :: data_dims(1) ! array containing the length of each dimension
        character(len=20) :: filename
        filename = "mvmul.h5"


        ! ALLOCATE ARRAYS
        allocate(mat(M,N))
        allocate(vec(N,O))
        allocate(prod(M,O))


        ! FILL UP ARRAYS
        prod = 0.0_dp
        call random_seed()
        call random_number(mat)
        call random_number(vec)


        ! open HDF5 interface
        call h5open_f(error)
        
        ! open the HDF5 file
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

        space_rank = 1
        data_dims(1) = N

        ! open data space
        call h5screate_simple_f(space_rank, data_dims, dspace_id, error)

        ! create dataset
        call h5dcreate_f(file_id, "vec", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

        ! write data set
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vec(1:N,1), data_dims, error)

        call h5dclose_f(dset_id,error)   
        call h5sclose_f(dspace_id, error)
        call h5fclose_f(file_id, error)
        call h5close_f(error)


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
