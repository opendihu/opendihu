#include "utility/path.h"

#include <errno.h>
#include <libgen.h>
#include <string.h>
#include <sys/stat.h>

namespace Path {
int mkpath(const char *dir, mode_t mode) {
  struct stat sb;

  if (!dir) {
    errno = EINVAL;
    return -1;
  }

  if (!stat(dir, &sb))
    return 0;

  mkpath(dirname(strdupa(dir)), mode);

  return mkdir(dir, mode);
}

bool fileExists(const char *file) {
  struct stat buffer;
  return (stat(file, &buffer) == 0);
}
} // namespace Path
