#pragma once

#include <sys/stat.h>
#include <sys/types.h>

namespace Path {
int mkpath(const char *dir, mode_t mode = S_IRWXU);

bool fileExists(const char *file);
} // namespace Path
