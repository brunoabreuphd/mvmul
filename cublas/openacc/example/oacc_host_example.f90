!! This code is from the following NVIDIA HPC SDK Documentation:
!! https://docs.nvidia.com/hpc-sdk/compilers/fortran-cuda-interfaces/#cflib-blas-oacc-host
!!
program testcublas
        use openacc
        use cublas
        implicit none

        integer, parameter :: n = 1000
        ! host arrays
        integer :: a(n), b(n)
        ! cublas handle
        type(cublasHandle) :: handle
        ! helpers
        integer :: istat

        ! create handle
        istat = cublasCreate(handle)

        ! force OpenACC kernels and cuBLAS to use the OpenACC stream
        istat = cublasSetStream(handle, acc_get_cuda_stream(acc_async_sync))
        
        !$acc data copyout(a, b)
        !$acc kernels
        a = 1
        b = 2
        !$acc end kernels
        call sswap(n, a, 1, b, 1)
        call cublasSswap(n, a, 1, b, 1)
        !$acc end data

        ! compare
        if (all(a == 1) .and. all(b == 2)) then
                write(*,*) ' Test PASSED '
        else
                write(*,*) ' Test FAILED '
        endif

end program testcublas
