## Compiling CL2QCD

CL2QCD uses [CMake](https://cmake.org/) to produce makefiles and then to be compiled.
To avoid polluting the original repository, you are encouraged to create a `build` folder, for example inside the cloned `cl2qcd` repository.
This is a desired behaviour as it makes accidently adding build or execution results to the version control system less likely.
In principle, the following usual steps should be enough to create all the executables.

```bash
mkdir build
cd build
cmake ..
make
```

### External dependencies

However, before building CL2QCD as described above, you should make sure that all the following required libraries are installed and accessible on your system.
* [OpenCL](http://www.khronos.org/opencl)
* [LIME](http://usqcd.jlab.org/usqcd-docs/c-lime/)
* [libxml2](http://xmlsoft.org)
* [Boost](http://www.boost.org/)
* [GMP](http://gmplib.org/)
* [MPFR](http://www.mpfr.org/)
* [Nettle](http://www.lysator.liu.se/~nisse/nettle/)


The `C++` compiler must be capable of basic `C++11` features.
Using one of the following compilers is recommended.
* [GCC](https://gcc.gnu.org/) &ge; `4.7`
* [Clang](https://clang.llvm.org/) &ge; `3.2`


The following software is not required, but used to build documentation and by auxiliary scripts.
* [doxygen](http://www.stack.nl/~dimitri/doxygen/)
* [Python 2](http://python.org)


### Helping `CMake`

The build configuration uses some non-standard `cmake` scripts.

To find the OpenCL library a script called [`FindOpenCL.cmake`](https://gitlab.com/Marix/FindOpenCL/raw/master/FindOpenCL.cmake) is used.
Make sure it can be found by CMake.
One option is to symlink it into some path searched for by CMake, e.g.

```bash
mkdir -p ${HOME}/.cmake/modules
cd ${HOME}/.cmake/modules
wget https://gitlab.com/Marix/FindOpenCL/raw/master/FindOpenCL.cmake
```

and adding the following to your `.bashrc` file,

```bash
alias cmake='cmake -DCMAKE_MODULE_PATH=~/.cmake/modules'
```

If you have the LIME-Library installed in some non-standard location you should add it to your `LIBRARY_PATH`,
```bash
export LIBRARY_PATH=/path/to/lime/lib:${LIBRARY_PATH}
```

#### Different approach

Of course, you can also simply pass the proper paths to `cmake` when you run it to set up the `build` folder.
The variables that could be used, among others, are
* `CMAKE_PREFIX_PATH`;
* `CMAKE_LIBRARY_PATH`;
* `CMAKE_INCLUDE_PATH`;
* `CMAKE_PROGRAM_PATH`.

Refer to the [CMake documentation](https://cmake.org/documentation/) to know how to use which, an example could read like the following.

```bash
cmake -DCMAKE_MODULE_PATH=${HOME}/.cmake/modules
      -DCMAKE_LIBRARY_PATH='/opt/AMDAPPSDK-2.9-1/lib/x86_64;/opt/lib64/'
      -DCMAKE_PREFIX_PATH='/opt2/' ..
```
