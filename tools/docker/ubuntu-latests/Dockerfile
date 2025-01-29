FROM ubuntu:latest

# Install prerequisites
RUN apt-get update && \
    apt-get install -y \
    clang-format git

# git -> Dockerfile
# wget unzip -> Download dependencies
# build-essential petsc-dev bison flex libeigen3-dev libxml2-dev libboost-all-dev -> PETSc
# cmake -> base64
# libffi-dev -> Solves error at run time: ModuleNotFoundError: No module named '_ctypes’
# vim -> optional
# python-is-python3 -> make `python` point to python3 (needed for some scripts that call python instead of python3)

