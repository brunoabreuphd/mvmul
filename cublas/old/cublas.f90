module cublas
!
! This module defines the interface to NVIDIA's C code cublasDgemv
!
        interface cuda_mvmul
        !
        ! Parameters to this cublas function can be found here:
        ! https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-gemv
        !
                subroutine cuda_dgemv(cta, m, n, alpha, A, lda, x, incx, &
                                beta, y, incy) bind(C, name='cublasDgemv')
                        use iso_c_binding
                        implicit none
                        character(1,c_char), value :: cta
                        integer(c_int), value :: m, n, lda, incx, incy
                        real(c_double), value :: alpha, beta
                        real(c_double), device, dimension(lda,*) :: A
                        real(c_double), device, dimension(*) :: x, y
                
                end subroutine cuda_dgemv
        
        end interface

end module cublas
