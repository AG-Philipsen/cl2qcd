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
* [Boost](http://www.boost.org/) &ge; `1.60.0`
* [GMP](http://gmplib.org/)
* [MPFR](http://www.mpfr.org/)
* [Nettle](http://www.lysator.liu.se/~nisse/nettle/)


The `C++` compiler must be capable of basic `C++17` features.
Using one of the following compilers is recommended.
* [GCC](https://gcc.gnu.org/) &ge; `9.3`
* [Clang](https://clang.llvm.org/) &ge; `5.0`


The following software is not required, but used to build documentation and by auxiliary scripts.
* [doxygen](http://www.stack.nl/~dimitri/doxygen/)
* [Python 2](http://python.org)


### Helping `CMake`

Running `cmake` will automatically look for all dependency and libraries.

If you have any required library installed in some non-standard location you should add them to your `LIBRARY_PATH`, e.g.

```bash
export LIBRARY_PATH=/path/to/lime/lib:/path/to/nettle/lib:${LIBRARY_PATH}
```

to help **CMake** to find the correct files.
Keep in mind that **CMake** uses the `find_package` functionality to locate files and you might check the [online modules documentation](https://cmake.org/cmake/help/latest/manual/cmake-modules.7.html) in case you wish to know how to ask to look for files in a non standard location.
For instance, it might happen that you have an older Boost version installed in the standard system location, but you would like to use a more recent one that you installed locally in the `${HOME}/Programs` folder.
According to the [FindBoost module documentation](https://cmake.org/cmake/help/latest/module/FindBoost.html), it should be sufficient to specify the options `-DBOOST_ROOT:PATH=${HOME}/Programs` and probably `-DBoost_NO_SYSTEM_PATHS=ON` when running `cmake` or `ccmake`.

In principle, `cmake` should be able to locate OpenCL as well on modern architectures.
However, it might not be the case, e.g. if you are using the AMD SDK development kit.
If you get a failure during the `cmake` configuration, you might ask `cmake` to use an alternative procedure, by running `cmake -DCMAKE_MODULE_PATH=../cmake/optional ..` (from the `build` directory).
In this way `cmake` will not use the standard procedure, but a customized one shipped within the repository.
If the failure still persists, you need to resort to the different approach explained here below.


#### Different approach

Of course, you can also *simply* pass the proper paths to `cmake` or to `ccmake` when you run it to set up the `build` folder (without setting up any alias or modifying the environment variables).
The variables that could be used, among others, are
* `CMAKE_PREFIX_PATH`;
* `CMAKE_LIBRARY_PATH`;
* `CMAKE_INCLUDE_PATH`;
* `CMAKE_PROGRAM_PATH`.

Refer to the [CMake documentation](https://cmake.org/documentation/) to know how to use which, an example could read like the following.

```bash
cmake -DCMAKE_LIBRARY_PATH='/opt/AMDAPPSDK-2.9-1/lib/x86_64;/opt/lib64/'
      -DCMAKE_PREFIX_PATH='/opt2/' ..
```

If you decide to use such an approach, once figured out the correct paths to configure `CL2QCD`, then you might consider to add a function to your shell config file, so that you can always retrieve your command.
Using `bash` for demonstration purposes, such a function could read

```bash
function GetCommandToConfigureCL2QCD()
{
    printf "cmake -DCMAKE_PREFIX_PATH='/opt/rocm/opencl' ..\n"
}
```
where additional `-D` options might have to be added.
