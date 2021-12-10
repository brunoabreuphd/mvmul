module mvmul_h5
use hdf5
implicit none

contains
        subroutine read_mvmul_data(mat, vec, prod)
                double precision, dimension(:,:), allocatable, &
                        intent(out) :: mat, vec, prod
                character(len=20) :: filename
                integer :: error
                integer(hid_t) :: file_id, root_id, dset_id, dspace_id
                character(len=5) :: dset_name
                integer(hid_t), dimension(2) :: dims, maxdims

                filename="mvmul.h5"

                call h5open_f(error)
                call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
                call h5gopen_f(file_id, "/", root_id, error)

                dset_name = "vec"
                call h5dopen_f(root_id, dset_name, dset_id, error)
                call h5dget_space_f(dset_id, dspace_id,error)
                call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)
                allocate(vec(dims(1), dims(2)))
                call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vec, dims, error)
                call h5dclose_f(dset_id, error)

                dset_name = "mat"
                call h5dopen_f(root_id, dset_name, dset_id, error)
                call h5dget_space_f(dset_id, dspace_id,error)
                call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)
                allocate(mat(dims(1), dims(2)))
                call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mat, dims, error)
                call h5dclose_f(dset_id, error)

                dset_name = "prod"
                call h5dopen_f(root_id, dset_name, dset_id, error)
                call h5dget_space_f(dset_id, dspace_id,error)
                call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)
                allocate(prod(dims(1), dims(2)))
                call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, prod, dims, error)

                call h5dclose_f(dset_id, error)
                call h5gclose_f(root_id, error)
                call h5fclose_f(file_id, error)
                call h5close_f(error)

        end subroutine read_mvmul_data

end module mvmul_h5
