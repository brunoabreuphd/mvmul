program test
        use cudafor
        use cublas
        implicit none

        integer, parameter :: N = 1024
        real, device, allocatable :: x_d(:)
        real, allocatable :: x(:)
        integer, device :: k_d
        integer :: k

        integer :: istat
        type(cublasHandle) :: handle

        ! allocate
        allocate(x(N))
        allocate(x_d(N))

        ! create cublas handle
        istat = cublasCreate(handle)
        if (istat /= 0) then
                write(*,*) '** cublas handle creation failed **'
                stop
        endif

        ! initiate host vector
        call random_number(x)
        ! transfer it to device
        x_d = x
        
        ! call cublas
        istat = cublasIsamax_v2(handle, N, x_d, 1, k_d)
        if (istat /= 0) then
                write(*,*) '** cublas call failed **'
                stop
        endif

        ! transfer result to host
        k = k_d

        ! print out result
        write(*,*) 
        write(*,*) 'Index of max value in x is: ', k
        

        ! clean up
        deallocate(x)
        deallocate(x_d)
        istat = cublasDestroy(handle)

end program test
