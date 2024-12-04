FROM ubuntu:22.04

# Install prerequisites
RUN apt-get update && \
    apt-get install -y \
    git \
    wget unzip \
    build-essential petsc-dev bison flex libeigen3-dev libxml2-dev libboost-all-dev \
    cmake \
    libffi-dev \
    vim \
    python-is-python3

# git -> Dockerfile
# wget unzip -> Download dependencies
# build-essential petsc-dev bison flex libeigen3-dev libxml2-dev libboost-all-dev -> PETSc
# cmake -> base64
# libffi-dev -> Solves error at run time: ModuleNotFoundError: No module named '_ctypes’
# vim -> optional
# python-is-python3 -> make `python` point to python3 (needed for some scripts that call python instead of python3)

WORKDIR /workspace

# Clone opendihu
RUN git clone --branch update-docker --depth=1 https://github.com/opendihu/opendihu.git  
WORKDIR opendihu

# Build opendihu, debug and release target for core, build release unit tests and run tests, build all examples
RUN make release_without_tests; echo "Build and run tests";

