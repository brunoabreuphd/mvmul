program mvmul_oacc
        use openacc
        use cublas
        use mvmul_h5
        implicit none

        ! cublas dgemv variables
        integer :: t                            ! transposition
        integer :: lda, incx, incy              ! lead dim of matrix, increments
        integer :: m, n                         ! dims of matrix
        real(8) :: alpha, beta                  ! scalar multipliers
        ! cublas handle
        type(cublasHandle) :: handle
        ! helper
        integer :: istat

        ! variables from previous mvmul code
        real(8) :: acc_diff, abs_diff
        real(8), dimension(:,:), allocatable :: mat, vec, prod  ! arrays from h5 file
        real(8), dimension(:,:), allocatable :: prod_cublas     ! product from cublas
        integer, dimension(2) :: shp_mat, shp_vec, shp_prod     ! shapes read from h5 file
        ! helpers
        integer :: D, E, F
        integer :: i
        integer :: flag


        ! read data from h5 file, allocate array
        call read_mvmul_data(mat,vec,prod)
        shp_mat = SHAPE(mat)
        shp_vec = SHAPE(vec)
        shp_prod = SHAPE(prod)
        allocate(prod_cublas(shp_prod(1), shp_prod(2)))
        D = shp_mat(1)
        E = shp_mat(2)
        F = shp_vec(2)

        ! instantiate prod_cublas
        prod_cublas = 0.0
        ! set cublas dgemv variables
        t = 0
        m = D
        n = E
        lda = D
        incx = 1
        incy = 1
        alpha = 1.0
        beta = 0.0
 
        ! create cublas handle
        istat = cublasCreate(handle)
        if (istat /= 0) then
                write(*,*) ' ** cublas handle creation failed ** '
                stop
        endif
        ! force OpenACC and cublas to use same stream
        istat = cublasSetStream(handle, acc_get_cuda_stream(acc_async_sync))

        !$acc data copyin(mat,vec) copy(prod_cublas) copyout(istat)
        istat = cublasDgemv_v2(handle, t, m, n, alpha, mat, lda, &
                vec(:,1), incx, beta, prod_cublas(:,1), incy)
        !$acc end data
        if (istat /= 0) then
                write(*,*) ' ** cublas call failed ** '
                stop
        endif


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
                write(*,*) "Products are not equal. Accumulated difference (absolute values) is: ", acc_diff
        endif

        
        ! clean up
        istat = cublasDestroy(handle)
        deallocate(mat, vec, prod, prod_cublas)


end program mvmul_oacc
