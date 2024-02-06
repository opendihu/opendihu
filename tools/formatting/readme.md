# Formatting the code

In OpenDiHu we stick to the [LLVM style guide](https://llvm.org/docs/CodingStandards.html) for `.cpp`, `.hpp` and `.tpp` files. To make sure that the files comply with the style, we use the [clang-format](https://clang.llvm.org/docs/ClangFormat.html) tool to format the code.

### How to use clang-format?
1. Install clang-format: `sudo apt install clang-format`.
2. Run clang-format in the parent directory `. tools/formatting/format-all`. You may also run clang-format for a single file, e.g, `clang-format -i file.cpp`.
3. Check if a file is already formatted: `clang-format -dry-run -Werror file.cpp`.