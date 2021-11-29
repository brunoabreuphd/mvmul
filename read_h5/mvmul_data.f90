program mvmul_data
        use mvmul_h5
        implicit none

        double precision, dimension(:,:), allocatable :: mat, vec, prod
        integer :: i

        call read_mvmul_data(mat,vec,prod)

        do i = 1, 10
                write(*,*) vec(i,1), prod(i,1), mat(i,1)
        enddo

        deallocate(mat)
        deallocate(vec)
        deallocate(prod)

end program mvmul_data
