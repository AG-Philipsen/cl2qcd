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
You might consider to use `ccmake ..` as third step, in order to explore possible different configurations of the software.
Type, then, `h` in order to get help either about the selected variable or on the general usage.

### External dependencies

However, before building CL2QCD as described above, you should make sure that all the following required libraries are installed and accessible on your system.
* [OpenCL](http://www.khronos.org/opencl)
* [LIME](http://usqcd.jlab.org/usqcd-docs/c-lime/)
* [libxml2](http://xmlsoft.org)
* [Boost](http://www.boost.org/) &ge; `1.59.0`
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
One option is to download and place it into some path searched for by CMake, e.g.

```bash
mkdir -p ${HOME}/.cmake/modules
cd ${HOME}/.cmake/modules
wget https://gitlab.com/Marix/FindOpenCL/raw/master/FindOpenCL.cmake
```

and running, then, `cmake` always with the `-DCMAKE_MODULE_PATH=~/.cmake/modules` option.
You might add

```bash
alias cmake='cmake -DCMAKE_MODULE_PATH=~/.cmake/modules'
alias ccmake='ccmake -DCMAKE_MODULE_PATH=~/.cmake/modules'
```
to your shell config file, or simply remember to specify the option when configuring the code.

If you have the LIME or Nettle Libraries installed in some non-standard location you should add them to your `LIBRARY_PATH`,
```bash
export LIBRARY_PATH=/path/to/lime/lib:/path/to/nettle/lib:${LIBRARY_PATH}
```
Actually this is true in general and you should add **CMake** to find the correct files if you did install them in standard positions.
Keep in mind that **CMake** uses the find_package functionality to locate files and you might check the [online modules documentation](https://cmake.org/cmake/help/latest/manual/cmake-modules.7.html) in case you wish to know how to ask to look for files in a non standard location.
For instance, it might happen that you have an older Boost version installed in the standard system location, but you would like to use a more recent one that you installed locally in the `${HOME}/Programs` folder.
According to the [FindBoost module documentation](https://cmake.org/cmake/help/latest/module/FindBoost.html), it should be sufficient to specify the options `-D BOOST_ROOT:PATH=${HOME}/Programs` and probably `-D Boost_NO_SYSTEM_PATHS=ON` when running `cmake` or `ccmake`.

#### Different approach

Of course, you can also *simply* pass the proper paths to `cmake` or to `ccmake` when you run it to set up the `build` folder (without setting up any alias or modifying the environment variables).
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

If you decide to use such an approach, once figured out the correct paths to configure `CL2QCD`, then you might consider to add a function to your shell config file, so that you can always retrieve your command.
Using `bash` for demonstration purposes, such a function could read
```bash
function GetCommandToConfigureCL2QCD()
{
    echo 'cmake -DCMAKE_MODULE_PATH=${HOME}/.cmake/modules ..'
}
```
where additional `-D` options might have to be added.
