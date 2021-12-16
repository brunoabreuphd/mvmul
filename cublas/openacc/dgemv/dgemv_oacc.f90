program dgemv_oacc
        use openacc
        use cublas
        implicit none

        ! matrix, vectors dimensions
        integer, parameter :: ORD1 = 32
        integer, parameter :: ORD2 = 32
        integer, parameter :: ORD3 = 1

        ! cublas dgemv variables
        integer :: t                            ! transposition
        integer :: lda, incx, incy              ! leading dimension of matrix, increments
        integer :: m, n                         ! matrix dimensions
        real(8) :: alpha, beta                  ! y -> alpha*A*x + beta*y
        real(8), allocatable :: A(:,:)          ! matrix
        real(8), allocatable :: x(:), y(:)      ! vectors

        ! cublas handle
        type(cublasHandle) :: handle
        
        ! helpers
        integer :: istat


        ! allocate arrays
        allocate(A(ORD1, ORD2))
        allocate(x(ORD2))
        allocate(y(ORD1))

        ! create cublas handle
        istat = cublasCreate(handle)
        if (istat /= 0 ) then
                write(*,*) ' ** cublas handle creation failed ** '
                stop
        endif
        ! force OpenACC and cublas to use same stream
        istat = cublasSetStream(handle, acc_get_cuda_stream(acc_async_sync))

        ! instantiate arrays with random numbers at host
        call random_number(A)
        call random_number(x)
        y = 0.0
        ! instantiate cublas dgemv variables
        t = 0
        m = ORD1
        n = ORD2
        lda = ORD1
        incx = 1
        incy = 1
        alpha = 1.0
        beta = 0.0

        !$acc data copyin(A,x) copy(y)
        istat = cublasDgemv_v2(handle, t, m, n, alpha, A, lda, x, incx, beta, y, incy)
        !$acc end data

        write(*,*) '    y(1:5):'
        write(*,*) y(1:5)

        ! clean up
        istat = cublasDestroy(handle)
        deallocate(A,x,y)

end program dgemv_oacc
