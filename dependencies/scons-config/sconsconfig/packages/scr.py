from .Package import Package

scr_text = r"""
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>

#include "mpi.h"
#include "scr.h"

int checkpoint(int size_mb) {
  int rank;
  char tmp[256];
  char file[SCR_MAX_FILENAME];
  char dir[SCR_MAX_FILENAME];
  char dname;
  int rc;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  SCR_Start_checkpoint();

  sprintf(tmp, "rank_%d", rank);
  std::cout << "In: " << tmp << "\n";

  SCR_Route_file(tmp, file);

  std::cout << "Out: " << file << "\n";

  SCR_Complete_checkpoint(1);

  return 0;
}

int main(int argc, char **argv) {
  int size_mb = 1;

  MPI_Init(NULL, NULL);
  SCR_Configf("SCR_CHECKPOINT_INTERVAL=%d", 1);
  SCR_Configf("SCR_USER_NAME=%s", "scons");
  SCR_Init();

  int flag = 0;
  SCR_Need_checkpoint(&flag);
  if (flag) {
    checkpoint(size_mb);
  }

  SCR_Finalize();
  MPI_Finalize();

  return 0;
}
"""


class SCR(Package):
    def __init__(self, **kwargs):
        defaults = {
            "download_url": "https://github.com/LLNL/scr/archive/refs/tags/v3.0.1.tar.gz"
        }
        defaults.update(kwargs)
        super(SCR, self).__init__(**defaults)
        self.sub_dirs = [
            # could be either lib or lib64
            ("include", "lib"),
            ("include", "lib64")
        ]
        self.libs = ["scr"]
        self.headers = ["scr.h"]

        self.ext = '.cpp'
        self.check_text = scr_text
        self.static = False
        self.number_output_lines = 11540 * 1.45

    def check(self, ctx):
        env = ctx.env
        ctx.Message("Checking for SCR ...           ")

        self.set_build_handler(
            [
                "mkdir -p ${PREFIX}",
                "sed -i 's/#!\\/bin\\/bash/#!\\/usr\\/bin\\/env bash/g' ${SOURCE_DIR}/bootstrap.sh",
                "sed -i '/INSTALL_DIR=/d' ${SOURCE_DIR}/bootstrap.sh",
                "cd ${SOURCE_DIR} && \
                 INSTALL_DIR=${PREFIX} ./bootstrap.sh --tag && \
                 mkdir -pv build && cd build && \
                 cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=${PREFIX} .. && \
                 make && make install",
            ]
        )

        self.check_options(env)
        res = super(SCR, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
