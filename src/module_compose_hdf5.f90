module module_compose_hdf5
  use hdf5
  implicit none
  private

  public :: write_compose_hdf5

contains

  subroutine write_compose_hdf5(filename, nb, t, yq, mn, mp, &
       q1, q2, q3, q4, q5, q6, q7, cs2, &
       ye, yn, yp, yh2, yh3, yhe3, yhe4, ynuc, &
       anuc, znuc, abar)

    character(*), intent(in) :: filename
    real(8), intent(in) :: nb(:), t(:), yq(:)
    real(8), intent(in) :: mn, mp

    ! Internally all 3-D arrays are stored as (t, yq, nb).
    ! The HDF5 Fortran interface reverses dimensions in the file,
    ! so h5py/h5ls see (nb, yq, t), matching PyCompOSE DD2.h5.
    real(8), intent(in) :: q1(:,:,:), q2(:,:,:), q3(:,:,:), q4(:,:,:)
    real(8), intent(in) :: q5(:,:,:), q6(:,:,:), q7(:,:,:), cs2(:,:,:)
    real(8), intent(in) :: ye(:,:,:), yn(:,:,:), yp(:,:,:)
    real(8), intent(in) :: yh2(:,:,:), yh3(:,:,:), yhe3(:,:,:), yhe4(:,:,:)
    real(8), intent(in) :: ynuc(:,:,:), anuc(:,:,:), znuc(:,:,:), abar(:,:,:)

    integer(hid_t) :: file_id
    integer :: hdferr

    call check_shapes(nb, t, yq, q1, q2, q3, q4, q5, q6, q7, cs2, &
         ye, yn, yp, yh2, yh3, yhe3, yhe4, ynuc, anuc, znuc, abar)

    call h5open_f(hdferr)
    call check_hdf(hdferr, "h5open_f")

    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdferr)
    call check_hdf(hdferr, "creating "//trim(filename))

    call write_real_1d(file_id, "nb", nb, "baryon number density [fm^-3]")
    call write_real_1d(file_id, "t",  t,  "temperature [MeV]")
    call write_real_1d(file_id, "yq", yq, "charge fraction")

    call write_real_scalar(file_id, "mn", mn, "neutron mass [MeV]")
    call write_real_scalar(file_id, "mp", mp, "proton mass [MeV]")

    call write_real_3d(file_id, "Q1",  q1, "pressure over number density: p/nb [MeV]")
    call write_real_3d(file_id, "Q2",  q2, "entropy per baryon [kb]")
    call write_real_3d(file_id, "Q3",  q3, "scaled and shifted baryon chemical potential: mu_b/m_n - 1")
    call write_real_3d(file_id, "Q4",  q4, "scaled charge chemical potential: mu_q/m_n")
    call write_real_3d(file_id, "Q5",  q5, "scaled electron-lepton chemical potential: mu_l/m_n")
    call write_real_3d(file_id, "Q6",  q6, "scaled free energy per baryon: f/(nb*m_n) - 1")
    call write_real_3d(file_id, "Q7",  q7, "scaled internal energy per baryon: e/(nb*m_n) - 1")
    call write_real_3d(file_id, "cs2", cs2, "sound speed squared [c^2]")

    call write_real_3d(file_id, "Y[e]",   ye,   "electron")
    call write_real_3d(file_id, "Y[n]",   yn,   "neutron")
    call write_real_3d(file_id, "Y[p]",   yp,   "proton")
    call write_real_3d(file_id, "Y[H2]",  yh2,  "deuteron")
    call write_real_3d(file_id, "Y[H3]",  yh3,  "tritium")
    call write_real_3d(file_id, "Y[He3]", yhe3, "helium 3")
    call write_real_3d(file_id, "Y[He4]", yhe4, "alpha particle")
    call write_real_3d(file_id, "Y[N]",   ynuc, "average nucleus")

    call write_real_3d(file_id, "A[N]", anuc, "average mass number of representative heavy nucleus")
    call write_real_3d(file_id, "Z[N]", znuc, "average charge number of representative heavy nucleus")
    call write_real_3d(file_id, "Abar", abar, "average mass number")

    call h5fclose_f(file_id, hdferr)
    call check_hdf(hdferr, "h5fclose_f")

    call h5close_f(hdferr)
    call check_hdf(hdferr, "h5close_f")

  end subroutine write_compose_hdf5


  subroutine write_real_scalar(file_id, name, value, desc)
    integer(hid_t), intent(in) :: file_id
    character(*), intent(in) :: name, desc
    real(8), intent(in) :: value

    integer(hid_t) :: space_id, dset_id
    integer(hsize_t) :: dims(1)
    integer :: hdferr

    dims(1) = 1_hsize_t

    call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
    call check_hdf(hdferr, "creating scalar dataspace for "//trim(name))

    call h5dcreate_f(file_id, trim(name), H5T_NATIVE_DOUBLE, &
         space_id, dset_id, hdferr)
    call check_hdf(hdferr, "creating dataset "//trim(name))

    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, value, dims, hdferr)
    call check_hdf(hdferr, "writing dataset "//trim(name))
    call write_string_attribute(dset_id, "desc", desc)

    call h5dclose_f(dset_id, hdferr)
    call check_hdf(hdferr, "closing dataset "//trim(name))
    call h5sclose_f(space_id, hdferr)
    call check_hdf(hdferr, "closing dataspace "//trim(name))

  end subroutine write_real_scalar


  subroutine write_real_1d(file_id, name, a, desc)
    integer(hid_t), intent(in) :: file_id
    character(*), intent(in) :: name, desc
    real(8), intent(in) :: a(:)

    integer(hid_t) :: space_id, dset_id
    integer(hsize_t) :: dims(1)
    integer :: hdferr

    dims(1) = int(size(a, 1), hsize_t)

    call h5screate_simple_f(1, dims, space_id, hdferr)
    call check_hdf(hdferr, "creating dataspace for "//trim(name))

    call h5dcreate_f(file_id, trim(name), H5T_NATIVE_DOUBLE, &
         space_id, dset_id, hdferr)
    call check_hdf(hdferr, "creating dataset "//trim(name))

    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, a, dims, hdferr)
    call check_hdf(hdferr, "writing dataset "//trim(name))
    call write_string_attribute(dset_id, "desc", desc)

    call h5dclose_f(dset_id, hdferr)
    call check_hdf(hdferr, "closing dataset "//trim(name))
    call h5sclose_f(space_id, hdferr)
    call check_hdf(hdferr, "closing dataspace "//trim(name))

  end subroutine write_real_1d


  subroutine write_real_3d(file_id, name, a, desc)
    integer(hid_t), intent(in) :: file_id
    character(*), intent(in) :: name, desc
    real(8), intent(in) :: a(:,:,:)

    integer(hid_t) :: space_id, dset_id
    integer(hsize_t) :: dims(3)
    integer :: hdferr

    dims = [int(size(a,1),hsize_t), int(size(a,2),hsize_t), &
         int(size(a,3),hsize_t)]

    call h5screate_simple_f(3, dims, space_id, hdferr)
    call check_hdf(hdferr, "creating dataspace for "//trim(name))

    call h5dcreate_f(file_id, trim(name), H5T_NATIVE_DOUBLE, &
         space_id, dset_id, hdferr)
    call check_hdf(hdferr, "creating dataset "//trim(name))

    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, a, dims, hdferr)
    call check_hdf(hdferr, "writing dataset "//trim(name))
    call write_string_attribute(dset_id, "desc", desc)

    call h5dclose_f(dset_id, hdferr)
    call check_hdf(hdferr, "closing dataset "//trim(name))
    call h5sclose_f(space_id, hdferr)
    call check_hdf(hdferr, "closing dataspace "//trim(name))

  end subroutine write_real_3d


  subroutine write_string_attribute(obj_id, name, value)
    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name, value

    integer(hid_t) :: space_id, type_id, attr_id
    integer(hsize_t) :: dims(1)
    integer(size_t) :: nchar
    integer :: hdferr

    if (len_trim(value) == 0) return

    call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
    call check_hdf(hdferr, "creating attribute dataspace "//trim(name))

    call h5tcopy_f(H5T_FORTRAN_S1, type_id, hdferr)
    call check_hdf(hdferr, "copying string type for attribute "//trim(name))

    nchar = int(len_trim(value), size_t)
    call h5tset_size_f(type_id, nchar, hdferr)
    call check_hdf(hdferr, "setting string size for attribute "//trim(name))

    call h5acreate_f(obj_id, trim(name), type_id, space_id, attr_id, hdferr)
    call check_hdf(hdferr, "creating attribute "//trim(name))

    dims(1) = 1_hsize_t
    call h5awrite_f(attr_id, type_id, trim(value), dims, hdferr)
    call check_hdf(hdferr, "writing attribute "//trim(name))

    call h5aclose_f(attr_id, hdferr)
    call check_hdf(hdferr, "closing attribute "//trim(name))
    call h5tclose_f(type_id, hdferr)
    call check_hdf(hdferr, "closing attribute type "//trim(name))
    call h5sclose_f(space_id, hdferr)
    call check_hdf(hdferr, "closing attribute dataspace "//trim(name))
  end subroutine write_string_attribute


  subroutine check_shapes(nb, t, yq, q1, q2, q3, q4, q5, q6, q7, cs2, &
       ye, yn, yp, yh2, yh3, yhe3, yhe4, ynuc, anuc, znuc, abar)

    real(8), intent(in) :: nb(:), t(:), yq(:)
    real(8), intent(in) :: q1(:,:,:), q2(:,:,:), q3(:,:,:), q4(:,:,:)
    real(8), intent(in) :: q5(:,:,:), q6(:,:,:), q7(:,:,:), cs2(:,:,:)
    real(8), intent(in) :: ye(:,:,:), yn(:,:,:), yp(:,:,:)
    real(8), intent(in) :: yh2(:,:,:), yh3(:,:,:), yhe3(:,:,:), yhe4(:,:,:)
    real(8), intent(in) :: ynuc(:,:,:), anuc(:,:,:), znuc(:,:,:), abar(:,:,:)

    integer :: nt, ny, nn

    nt = size(t)
    ny = size(yq)
    nn = size(nb)

    call assert_shape_3d("Q1",    q1,   nt, ny, nn)
    call assert_shape_3d("Q2",    q2,   nt, ny, nn)
    call assert_shape_3d("Q3",    q3,   nt, ny, nn)
    call assert_shape_3d("Q4",    q4,   nt, ny, nn)
    call assert_shape_3d("Q5",    q5,   nt, ny, nn)
    call assert_shape_3d("Q6",    q6,   nt, ny, nn)
    call assert_shape_3d("Q7",    q7,   nt, ny, nn)
    call assert_shape_3d("cs2",   cs2,  nt, ny, nn)
    call assert_shape_3d("Y[e]",  ye,    nt, ny, nn)
    call assert_shape_3d("Y[n]",  yn,    nt, ny, nn)
    call assert_shape_3d("Y[p]",  yp,    nt, ny, nn)
    call assert_shape_3d("Y[H2]", yh2,   nt, ny, nn)
    call assert_shape_3d("Y[H3]", yh3,   nt, ny, nn)
    call assert_shape_3d("Y[He3]",yhe3,  nt, ny, nn)
    call assert_shape_3d("Y[He4]",yhe4,  nt, ny, nn)
    call assert_shape_3d("Y[N]",  ynuc,  nt, ny, nn)
    call assert_shape_3d("A[N]",  anuc,  nt, ny, nn)
    call assert_shape_3d("Z[N]",  znuc,  nt, ny, nn)
    call assert_shape_3d("Abar",  abar,  nt, ny, nn)

  end subroutine check_shapes


  subroutine assert_shape_3d(name, a, n1, n2, n3)
    character(*), intent(in) :: name
    real(8), intent(in) :: a(:,:,:)
    integer, intent(in) :: n1, n2, n3

    if (size(a,1) /= n1 .or. size(a,2) /= n2 .or. size(a,3) /= n3) then
       write(*,'(a,1x,a,1x,3(i0,1x),a,1x,3(i0,1x))') &
            "ERROR: wrong shape for", trim(name), &
            size(a,1), size(a,2), size(a,3), &
            "expected", n1, n2, n3
       error stop
    endif
  end subroutine assert_shape_3d


  subroutine check_hdf(hdferr, where)
    integer, intent(in) :: hdferr
    character(*), intent(in) :: where

    if (hdferr < 0) then
       write(*,'(a,1x,a,1x,i0)') "HDF5 error in", trim(where), hdferr
       error stop
    endif
  end subroutine check_hdf

end module module_compose_hdf5
