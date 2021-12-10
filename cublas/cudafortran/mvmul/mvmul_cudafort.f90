program mvmul_cudafort
        use cudafor
        use cublas
        use mvmul_h5
        implicit none

        ! cuda and cublas dgemv variables
        ! HOST
        integer :: t                    ! transposition
        integer :: lda, incx, incy      ! leading dimension of matrix, increments
        integer :: m, n                 ! indexes of matrix
        real(8) :: alpha, beta             ! scalars for y -> alpha*A*x + beta*y
        ! DEVICE
        real(8), device :: alpha_d, beta_d
        ! helpers
        type(cudaEvent) :: startEvent, stopEvent
        integer :: istat
        real :: time
        ! cublas handle
        type(cublasHandle) :: handle 


        ! variables from previous mvmul code
        real(8) :: acc_diff, abs_diff
        real(8), dimension(:,:), allocatable :: mat, vec, prod  ! arrays from h5 file
        real(8), dimension(:,:), allocatable :: prod_cublas        ! product from cublas
        integer, dimension(2) :: shp_mat, shp_vec, shp_prod
        integer :: D, E, F
        integer :: i
        integer :: flag

        ! device arrays
        real(8), device, dimension(:,:), allocatable :: mat_d, vec_d, prod_d


        ! read data from h5 file and allocate host arrays
        call read_mvmul_data(mat,vec,prod)
        shp_mat = SHAPE(mat)
        shp_vec = SHAPE(vec)
        shp_prod = SHAPE(prod)
        allocate(prod_cublas(shp_prod(1), shp_prod(2)))
        D = shp_mat(1)
        E = shp_mat(2)
        F = shp_vec(2)

        ! making the translation into cublasDgemv variables
        prod_cublas = 0.0
        m = D
        n = E
        lda = D
        incx = 1
        incy = 1
        alpha = 1.0
        beta = 0.0
        t = 0
        

        ! allocate device arrays
        allocate(mat_d(D,E))
        allocate(vec_d(E,F))
        allocate(prod_d(D,F))

        ! create cublas handle
        istat = cublasCreate(handle)
        if (istat /= 0) then
                write(*,*) '** cublas handle creation failed **'
                stop
        endif
        ! create cuda events
        istat = cudaEventCreate(startEvent)
        istat = cudaEventCreate(stopEvent)


        ! transfer everything to device
        mat_d = mat
        vec_d = vec
        prod_d = prod
        alpha_d = alpha
        beta_d = beta
        
        ! call cublas
        istat = cudaEventRecord(startEvent,0)
        istat = cublasDgemv_v2(handle, t, m, n, alpha_d, mat_d, lda, vec_d(:,1), incx, beta_d, prod_d(:,1), incy)
        istat = cudaEventRecord(stopEvent, 0)
        istat = cudaEventSynchronize(stopEvent)
        istat = cudaEventElapsedTime(time, startEvent, stopEvent)
        if (istat /= 0) then
                write(*,*) '** cublas call failed **'
                stop
        endif
        write(*,*) 'CUBLAS DGEMV took: ', time, ' ms'

        ! transfer result to host
        prod_cublas = prod_d
        

        ! compare results
        flag = 0
        acc_diff = 0.0
        do i = 1, E
                abs_diff = abs(prod(i,1) - prod_cublas(i,1))
                if (abs_diff /= 0.0) then
                        acc_diff = acc_diff + abs_diff
                        flag = 1
                endif
        enddo
        if (flag == 0) then
                write(*,*) "Products are equal!"
        else
                write(*,*) "Products are not equal. Accumulated difference (absloute values) is: ", acc_diff
        endif


        ! clean up
        deallocate(mat, vec, prod, prod_cublas)
        deallocate(mat_d, vec_d, prod_d)
        istat = cublasDestroy(handle)

end program mvmul_cudafort
