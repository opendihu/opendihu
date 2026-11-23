import sys, os, multiprocessing, subprocess, socket
from .Package import Package

#
# Python requires bzip2 (for matplotlib to work)
#
class Python(Package):

  def __init__(self, **kwargs):
    defaults = {
      'download_url': 'https://www.python.org/ftp/python/3.12.16/Python-3.12.16.tgz'
    }

    defaults.update(kwargs)
    super(Python, self).__init__(**defaults)
    self.ext = '.cpp'
    self.sub_dirs = [
        ('include', 'lib'),
    #    ('include/mysql', 'lib64'),
    ]
    #self.headers = ['mysql.h']
    #self.libs = ['mysqlclient']
    #self.extra_libs = ['lapack', 'blas']
    self.check_text = r'''
      #include <stdlib.h>
      #include <stdio.h>
      #include <Python.h>
      int main(int argc, char* argv[])
      {
        Py_Initialize();
        return EXIT_SUCCESS;
      }
'''

  def use_existing_python(self, ctx):
    """
    Configure an already installed python (e.g. the system python) instead of downloading one.
    The interpreter is "python3" from the PATH, or PYTHON_DIR/bin/python3 if PYTHON_DIR is given.
    Its include directory, library and prefix are queried with sysconfig, and a python_home.h
    with the prefix of this python is generated, which is needed by the opendihu core.
    Returns the python version as string, e.g. "3.12", or None on failure.
    """
    env = ctx.env
    python_dir = self.get_option(env, 'PYTHON_DIR')
    executable = os.path.join(python_dir, 'bin', 'python3') if python_dir else 'python3'

    query = 'import sys, sysconfig; print(sysconfig.get_config_var("VERSION")); print(sys.base_prefix); ' \
      'print(sysconfig.get_paths()["include"]); print(sysconfig.get_config_var("LIBDIR"))'
    try:
      output = subprocess.check_output([executable, '-I', '-c', query], universal_newlines=True)
      version, prefix, include_dir, lib_dir = output.strip().split('\n')
    except Exception as e:
      ctx.Log('Could not query python interpreter "{}": {}\n'.format(executable, e))
      return None
    ctx.Log('Using existing python "{}": version {}, prefix {}, include dir {}, lib dir {}\n'.format(
      executable, version, prefix, include_dir, lib_dir))

    # write python_home.h to a directory inside opendihu, because the python include directory is usually not writable
    python_home_include_dir = os.path.join(self.determine_dependencies_dir(ctx), 'python', 'existing_python', 'include')
    os.makedirs(python_home_include_dir, exist_ok=True)
    with open(os.path.join(python_home_include_dir, 'python_home.h'), 'w') as f:
      f.write('#define PYTHON_HOME_DIRECTORY "{}"\n'.format(prefix))

    # let the package search use exactly these directories (PYTHON_INC_DIR/PYTHON_LIB_DIR replace PYTHON_DIR)
    if 'PYTHON_DIR' in env:
      del env['PYTHON_DIR']
    env['PYTHON_INC_DIR'] = include_dir + ';' + python_home_include_dir
    env['PYTHON_LIB_DIR'] = lib_dir
    self.libs = ['python' + version]
    self.existing_python_configured = True
    return version

  def check(self, ctx):
    env = ctx.env
    self.headers = ["Python.h"]

    # download and build python, unless PYTHON_DOWNLOAD = False, then use the existing (system) python
    download = self.have_any_options(env, 'PYTHON_DOWNLOAD', 'DOWNLOAD_ALL', 'PYTHON_REBUILD')
    user_given_paths = not getattr(self, 'existing_python_configured', False) \
      and self.have_any_options(env, 'PYTHON_INC_DIR', 'PYTHON_LIB_DIR', 'PYTHON_LIBS')
    if download or user_given_paths:
      ctx.Message('Checking for Python 3.12 ...   ')
      self.libs = ["python3.12"]
    else:
      version = self.use_existing_python(ctx)
      ctx.Message('Checking for Python {} ...'.format(version or '').ljust(31))
      if version is None:
        self.check_required(False, ctx)
        ctx.Result(False)
        return False

    self.set_build_handler([
      'mkdir -p ${PREFIX}',
      'cd ${SOURCE_DIR} && chmod +x ./configure && ./configure --enable-shared --prefix=${PREFIX} \
        LDFLAGS="-Wl,--rpath=${PREFIX}/lib -L${DEPENDENCIES_DIR}/bzip2/install/lib -L${DEPENDENCIES_DIR}/zlib/install/lib" \
        CPPFLAGS="-I${DEPENDENCIES_DIR}/bzip2/install/include -I${DEPENDENCIES_DIR}/zlib/install/include" CXX=g++ CC=gcc \
        && make && make install',
      '$export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PREFIX}/lib',
      'cd ${PREFIX}/include && echo "#define PYTHON_HOME_DIRECTORY \\"${PREFIX}\\"\n" > python_home.h',
    ])
    self.number_output_lines = 11540*1.45

    self.check_options(env)
    res = super(Python, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
