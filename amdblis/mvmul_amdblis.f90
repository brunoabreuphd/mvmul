program mvmul_amdblis
        use mvmul_h5
        use, intrinsic :: iso_fortran_env
        implicit none

        integer, parameter :: dp = REAL64 ! double precision
        integer, parameter :: i32 = INT32 ! 32-bit integer

        real(dp) :: acc_diff, abs_diff
        real(dp), dimension(:,:), allocatable :: mat, vec, prod  ! arrays from h5 file
        real(dp), dimension(:,:), allocatable :: prod_blis        ! product from blis
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
        allocate(prod_blis(shp_prod(1), shp_prod(2)))
        M = shp_mat(1)
        N = shp_mat(2)
        O = shp_vec(2)

        ! call blis
        call cpu_time(startT)
        call dgemv('N', M, N, 1.0_dp, mat, M, vec, 1, 0.0_dp, prod_blis, 1)
        call cpu_time(endT)
        write(*,*) 'BLIS DGEMV took: ', (endT-startT), ' s'

        flag = 0
        acc_diff = 0.0_dp
        do i = 1, N
                abs_diff = abs(prod(i,1) - prod_blis(i,1))
                if (abs_diff /=  0.0_dp) then
                        acc_diff = acc_diff + abs_diff
                        flag = 1
                endif
        enddo
        if (flag == 0) then
                write(*,*) "Products are equal!"
        else
                write(*,*) "Products are not equal. Acc diff is: ", acc_diff
        endif

        deallocate(mat)
        deallocate(vec)
        deallocate(prod)
        deallocate(prod_blis)

end program mvmul_amdblis
