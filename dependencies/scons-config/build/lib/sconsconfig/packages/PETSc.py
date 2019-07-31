import sys, os
from distutils import sysconfig
from Package import Package

petsc_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <petsc.h>
int main(int argc, char* argv[]) {
   PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
   printf("%d\n", MPI_VERSION);
   printf("%d\n", MPI_SUBVERSION);
   PetscFinalize();
   return EXIT_SUCCESS;
}
'''


def parse_conf(ctx, conf_path, lib_dirs, libs):
    vars = {}
    sysconfig.parse_makefile(conf_path, vars)
    flag_dict = ctx.env.ParseFlags(vars['PACKAGES_LIBS'])
    lib_dirs.extend(flag_dict['LIBPATH'])
    for ii in range(len(libs)):
        libs[ii].extend(flag_dict['LIBS'])

def find_conf(ctx, base, inc_dirs, lib_dirs, libs, extra_libs):
    # PETSc 3.1
    conf_path = os.path.join(base, 'conf', 'petscvariables')
    if os.path.exists(conf_path):
        parse_conf(ctx, conf_path, lib_dirs, libs)

    # PETSC 2.3.3
    conf_path = os.path.join(base, 'bmake', 'petscconf')
    if os.path.exists(conf_path):
        vars = {}
        sysconfig.parse_makefile(conf_path, vars)
        if 'PETSC_ARCH' in vars:
            arch = vars['PETSC_ARCH']
            inc_dirs.extend([os.path.join(base, 'bmake', arch)])
            lib_dirs.extend([os.path.join(base, 'lib', arch)])
            conf_path = os.path.join(base, 'bmake', arch, 'petscconf')
            parse_conf(ctx, conf_path, lib_dirs, libs)

class PETSc(Package):
  
    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.6.tar.gz',
        }
        defaults.update(kwargs)
        super(PETSc, self).__init__(**defaults)
        #self.ext = '.c'
        #self.sub_dirs = [
        #    ('include/mysql', 'lib'),
        #    ('include/mysql', 'lib64'),
        #]
        #self.headers = ['mysql.h']
        self.libs = [['petsc'], ['petscksp', 'petscvec', 'petsc']]
        self.extra_libs = ['blas']        # the system tries to include one of them after other, if linking else fails
        self.check_text = petsc_text
        self.static = False
        
        # Setup the build handler.
        self.set_build_handler([
            './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=no \
            --with-lapack-lib=${LAPACK_DIR}/lib/liblapack.so\
            --with-blas-lib=${LAPACK_DIR}/lib/libblas.so',
            'make all',     # do not add -j option, because it is not supported by Makefile of PETSc
            'make install',
            'make test',
        ])

        self.number_output_lines = 1850
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for PETSc ... ')
        self.check_options(env)

        res = super(PETSc, self).check(ctx, loc_callback=find_conf)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]