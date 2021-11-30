program mvmul_intelmkl
        use mvmul_h5
        use mkl_service
        use, intrinsic :: iso_fortran_env
        implicit none
        include "mkl_lapack.fi"

        integer, parameter :: dp = REAL64 ! double precision
        integer, parameter :: i32 = INT32 ! 32-bit integer

        real(dp), dimension(:,:), allocatable :: mat, vec, prod  ! arrays from h5 file
        real(dp), dimension(:,:), allocatable :: prod_mkl        ! product from mkl
        integer, dimension(2) :: shp_mat, shp_vec, shp_prod
        integer(i32) :: M, N, O
        integer :: i
        integer :: flag
        real(dp) :: startT, endT

        ! read vectors from h5 file and allocate for new vec
        call read_mvmul_data(mat,vec,prod)
        shp_mat = SHAPE(mat)
        shp_vec = SHAPE(vec)
        shp_prod = SHAPE(prod)
        allocate(prod_mkl(shp_prod(1), shp_prod(2)))
        M = shp_mat(1)
        N = shp_mat(2)
        O = shp_vec(2)

        ! call mkl
        call cpu_time(startT)
        call dgemv('N', M, N, 1.0_dp, mat, M, vec, 1, 0.0_dp, prod_mkl, 1)
        call cpu_time(endT)
        write(*,*) 'MKL DGEMV took: ', (endT-startT), ' s'

        flag = 0
        do i = 1, N
                if (prod(i,1) /=  prod_mkl(i,1)) then
                        write(*,*) "Products are not the same!", prod(i,1), prod_mkl(i,1)
                        flag = 1
                endif
        enddo
        if (flag == 0) then
                write(*,*) "Products are equal!"
        endif

        deallocate(mat)
        deallocate(vec)
        deallocate(prod)
        deallocate(prod_mkl)

end program mvmul_intelmkl
