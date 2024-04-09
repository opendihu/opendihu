import sys, os, multiprocessing
from .Package import Package
import subprocess

class precice(Package):

  def __init__(self, **kwargs):
    defaults = {
        'download_url': 'https://github.com/precice/precice/archive/refs/tags/v3.1.0.zip',
    }
    defaults.update(kwargs)
    super(precice, self).__init__(**defaults)
    self.ext = '.cpp'
    self.set_rpath = True
    self.check_text = r'''
    #include <iostream>
    #include <cstdlib>
    #include <fstream>
    #include <precice/precice.hpp>
    
    int main()
    {
      std::ofstream file("install_precice-config.xml");
      file << R"(<?xml version="1.0"?>
<precice-configuration>
    
    <data:scalar name="Data"/>
    <mesh name="Mesh1" dimensions="3">
      <use-data name="Data"/>
    </mesh>
    <mesh name="Mesh2" dimensions="3">
      <use-data name="Data"/>
    </mesh>
    
    <participant name="Participant1">
      <provide-mesh name="Mesh1"/>
      <write-data name="Data" mesh="Mesh1"/>    
    </participant>

    <participant name="Participant2">
      <receive-mesh name="Mesh1" from="Participant1"/>
      <provide-mesh name="Mesh2"/>
      <read-data name="Data" mesh="Mesh2"/>
      <mapping:nearest-neighbor
        direction="read"
        from="Mesh1"
        to="Mesh2"
        constraint="consistent" />
    </participant>
    
    <m2n:sockets acceptor="Participant1" connector="Participant2" network="lo" />
    <coupling-scheme:serial-explicit>
      <participants first="Participant1" second="Participant2"/>
      <time-window-size value="0.01"/>
      <max-time value="0.05"/>
      <exchange data="Data" mesh="Mesh1" from="Participant1" to="Participant2"/>
    </coupling-scheme:serial-explicit>

</precice-configuration>
)";
      file.close();
    
      precice::Participant participant("Participant1","install_precice-config.xml",0,1);
      //participant.initialize();
      //participant.finalize();
      
      int ret = system("rm -f install_precice-config.xml");
    
      return EXIT_SUCCESS;
    }
'''

    self.number_output_lines = 15536*1.5
      
    self.libs = ['precice']
    self.extra_libs = [[],
                       ['boost_filesystem', 'boost_log_setup', 'boost_log', 'boost_program_options', 'boost_system', 'boost_thread', 'boost_unit_test_framework', 'dl', 'boost_regex'],
    									 ['boost_atomic', 'boost_chrono', 'boost_filesystem', 'boost_log_setup', 'boost_log', 'boost_prg_exec_monitor', 'boost_program_options', 'boost_system', 'boost_test_exec_monitor', 'boost_thread', 'boost_unit_test_framework', 'dl', 'boost_regex']]    
    self.headers = ["precice/precice.hpp"]

  def check(self, ctx):
    env = ctx.env
    ctx.Message('Checking for preCICE ...       ')
    self.check_options(env)

    if True:
      # first try is using a system boost
      self.set_build_handler([
        'mkdir -p ${PREFIX}/include',

        # precice
        'cd ${SOURCE_DIR} && mkdir -p build && cd build && '+ctx.env["cmake"]+' -DCMAKE_INSTALL_PREFIX=${PREFIX} \
          -DCMAKE_BUILD_TYPE=Release \
          -DPRECICE_FEATURE_PYTHON_ACTIONS=OFF \
          ..',
        'cd ${SOURCE_DIR}/build && make precice install'
      ])
      
      res = super(precice, self).check(ctx)
    
    
    self.check_options(env)
    res = super(precice, self).check(ctx)
  

    # try again downloading libxml2 -> this is necessary for the ipvs epyc cluster
    # to run in the ipvs_epyc cluster, you also need to add 
    # export PKG_CONFIG_PATH="/path/to/opendihu/dependencies/petsc/install/lib/pkgconfig:${PKG_CONFIG_PATH}"



    if not res[0]:
      ctx.Log('Retry (1) with manually building libxml2\n')
      self.set_build_handler([
        'mkdir -p ${PREFIX}/include',

        # libxml2
        'cd ${SOURCE_DIR} && if [ ! -f ${SOURCE_DIR}/libxml2-2.9.9.tar.gz ]; then wget ftp://xmlsoft.org/libxml2/libxml2-2.9.9.tar.gz; fi; \
         tar xf libxml2-2.9.9.tar.gz; cd libxml2-2.9.9; ./configure --prefix=${PREFIX} --without-python && make install -j 16',

        # precice
        'cd ${SOURCE_DIR} && mkdir -p build && cd build && '+ctx.env["cmake"]+' -DCMAKE_INSTALL_PREFIX=${PREFIX} \
          -DCMAKE_BUILD_TYPE=Release \
          -DPRECICE_FEATURE_PYTHON_ACTIONS=OFF \
          -DLIBXML2_INCLUDE_DIR=${PREFIX}/include/libxml2 -DLIBXML2_LIBRARY=${PREFIX}/lib/libxml2.so \
          ..',
        'cd ${SOURCE_DIR}/build && make precice install'
      ])
      
      self.check_options(env)
      res = super(precice, self).check(ctx)

    if not res[0]:
      ctx.Log('\n\nInstallation of preCICE failed. Rebuild with\n  make clean; scons PRECICE_REBUILD=True\n\n')
      
    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
