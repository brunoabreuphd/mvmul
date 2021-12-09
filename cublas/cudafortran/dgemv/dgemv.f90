program dgemv
        !! https://docs.nvidia.com/hpc-sdk/compilers/fortran-cuda-interfaces/index.html#dp-dgemv
        use cudafor
        use cublas
        implicit none

        ! matrix dimensions
        integer, parameter :: ORD1 = 32
        integer, parameter :: ORD2 = 32
        integer, parameter :: ORD3 = 1

        ! cublas dgemv variables
        ! host
        integer :: t                    ! transposition
        integer :: lda, incx, incy      ! leading dimension of matrix, increments
        integer :: m, n                 ! indexes of matrix
        real(8) :: alpha, beta             ! scalars for y -> alpha*A*x + beta*y
        real(8), allocatable :: A(:,:)     ! matrix
        real(8), allocatable :: x(:), y(:) ! vectors
        ! device
        integer, device :: t_d
        integer, device :: lda_d, incx_d, incy_d
        integer, device :: m_d, n_d
        real(8), device :: alpha_d, beta_d
        real(8), device, allocatable :: A_d(:,:)
        real(8), device, allocatable :: x_d(:), y_d(:)

        ! helpers
        integer :: istat
        ! cublas handle
        type(cublasHandle) :: handle 


        ! allocate
        allocate(A(ORD1,ORD2), x(ORD2), y(ORD1))
        allocate(A_d(ORD1,ORD2), x_d(ORD2), y_d(ORD1))

        ! create cublas handle
        istat = cublasCreate(handle)
        if (istat /= 0) then
                write(*,*) '** cublas handle creation failed **'
                stop
        endif

        ! initiate host arrays
        call random_number(A)
        call random_number(x)
        y = 0.0
        ! initiate host variables
        t = 0
        m = ORD1
        n = ORD2
        lda = ORD1
        incx = 1
        incy = 1
        alpha = 1.0
        beta = 0.0

        ! transfer everything to device
        A_d = A
        x_d = x
        y_d = y
        alpha_d = alpha
        beta_d = beta
        
        ! call cublas
        istat = cublasDgemv_v2(handle, t, m, n, alpha_d, A_d, lda, x_d, incx, beta_d, y_d, incy)
        if (istat /= 0) then
                write(*,*) '** cublas call failed **'
                stop
        endif

        ! transfer result to host
        y = y_d
        
        ! print out
        write(*,*) 'y(1:10) = '
        write(*,*) y(1:10)


        ! clean up
        deallocate(A, x, y)
        deallocate(A_d, x_d, y_d)
        istat = cublasDestroy(handle)

end program dgemv
