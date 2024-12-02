import sys, os
from .Package import Package


class HDF5(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.4/src/hdf5-1.14.4-3.tar.gz',
        }
        defaults.update(kwargs)
        super(HDF5, self).__init__(**defaults)
        self.parallel = kwargs.get('parallel', True)
        self.libs=['hdf5']
        self.extra_libs=[
            [], ['z'], ['m'], ['z', 'm'],
        ]
        if self.parallel:
            self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <mpi.h>

#define FILE "testfile.h5"

int main(int argc, char* argv[]) {
  hid_t plist_id, file_id, mem_space, file_space, dset_id;
  hsize_t dims[1];
  dims[0] = 10;
  int data[10];
  int rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  file_space = H5Screate_simple(1, dims, NULL);
  dset_id = H5Dcreate(file_id, "test", H5T_NATIVE_INT, file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, H5T_NATIVE_INT, H5P_DEFAULT, file_space, plist_id, data);

  H5Pclose(plist_id);
  H5Dclose(dset_id);
  H5Sclose(file_space);
  H5Fclose(file_id);
  MPI_Finalize();

  remove(FILE);
  return EXIT_SUCCESS;
}
'''
        else:
            self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
int main(int argc, char* argv[]) {
   hid_t plist_id, file_id, space;
   hsize_t dims[1];
   dims[0] = 10;
   space = H5Screate_simple(1, dims, NULL);
   return EXIT_SUCCESS;
}
'''

        # Setup the build handler. I'm going to assume this will work for all architectures.
        self.set_build_handler([
            "mkdir -p ${PREFIX}",
            "cd ${SOURCE_DIR} && \
            mkdir -pv build && cd build && \
            cmake -DCMAKE_BUILD_TYPE=RELEASE \
                  -DHDF5_BUILD_CPP_LIB=OFF \
                  -DHDF5_ENABLE_PARALLEL=ON \
                  -DBUILD_SHARED_LIBS=ON \
                  -DCMAKE_INSTALL_PREFIX=${PREFIX} \
                  CC=mpicc \
                  .. && \
            make && make install",
        ])

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for HDF5 ...          ')
        self.check_options(env)

        res = super(HDF5, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
