# raytrace_cpu
A CPU port of the CUDAKerr ray tracing code (Wilkins et al. 2012, 2016). While CUDAKerr was written to run on NVIDIA GPUs, this version provides the same functionality running on the CPU, and can be compiled to run on a range of platforms.

CUDAKerr and raytrace_cpu predict light paths in the curved spacetime around black holes by calculating the null geodesics in the Kerr spacetime. While CUDAKerr is a general ray tracing code, it specialises in creating models of X-ray reflection and reverberation from the accretion discs around black holes.

Ray tracing simulations are built from three types of component:
1. The core Raytracer class (src/raytrace/raytracer.h and raytracer.cpp). You do not need to access this library directly, since it is implemented by the ray source classes.
2. Ray source classes that define where rays are started (e.g. PointSource, ImagePlane, also in src/raytracer)
3. Simulation applications that implement one of the ray source classes, run the ray trace from it, then run a series of computations on the result, to produce the final scientific output.

Typically, a user will only create or interact with number (3).

All of the source code is organised in the src/ directory, and arranged into subdirectories that group together code for specific scientific applications, for example simulating emissivity profiles of accretion discs (src/emissivity), or creating ray traced images of the accretion disc (imageplane/). The full GPU version of CUDAKerr contains the full suite of simlations, only a subset are provided in this CPU version.

## Prerequesites
- A C++ compiler (GCC or the Apple-provided clang)
- cmake
- cfitsio
- GSL (only required for some of the tools, can comment out requirement in CMakeLists.txt)

This code should, for the most part, be platform-independent. Linux and Mac are supported, and there is no reason you cannot compile it on Windows or any other platform.

All of these packages are most easily installed from your package manager on Linux, or via MacPorts on Mac.

## Building rayrace_cpu
Compilation of the raytrace_cpu code is handled by the cmake build system, which you will need to have installed. The cmake configuration is contained within the file CMakeLists.txt. There is such a file in the top level directory which contains the global build options, then each of the subdirectories within src/ have their own CMakeLists.txt file which defines the different simulation applications for each directory

You can create a file in the top level directory named `CMakeLocal.txt`. This file allows you to set specific cmake options for each system the code will be compiled on, and therefore is not synced by git. Use this, for example, if you get library not found errors and need to explicitly pont cmake to some of the required libraries.

The following steps assume your terminal is in the top level `raytrce_cpu` directory.

The first time you want to compile on your system, configure cmake for raytrace_cpu by running:

`cmake .`

You will only need to run this command the one time, unless you change any of the options in any of the CMakeLists.txt files. If you get any errors, check that you have the prerequisites installed. If you do and you still get errors, you may need to create a CMakeLocal.txt file (see above) to explicitly point cmake to the libraries on your system.

Once cmake is set up, you can then use it to build any of the simulation applications. For example, to build `emissivity`, run:

`$ cmake --build . --target emissivity`

The executable files, once compiled, are placed in the `bin/` directory, where usually you should execute them from.

## Parameter files
Each ray tracing simulation will have a number of variable parameters that need to be run before you execute the code. Usually these are placed in parameter files, within the `par/` directory (in the top level raytrace_cpu directory). This directory is not synced by git, so you will probably need to create it the first time around, however you can find some example parameter files in the `par_example` directory. 

The parameter files are plain text files with one parameter per line, of the format `name = value`, for example `spin = 0.998`.

By default, each simulation will look for a file within `par/` that corresponds to its own name with the extension `.par`, for example `emissivity` will look for `par/emissivity.par`, however you can manually specify a parameter file each time you run the application using, for example `./emissivity --parfile my_parameters.par`.

There are some parameters that you need to explicitly define, otherwise you will get an error when you run the code. For others, if the parameter is missing from the parameter file, a default value will be assumed. You can see the parameters defined for each simulation near the top of the `main()` function in its source code. Look for statements along the lines of `par_file.get_parameter<double>("spin")`. This means search for a parameter named `spin` and interpret its value as a double-precision floating point number.

There are also a few parameters that can be overridden on the command line (look for statements in the code along the lines of `par_args.get_parameter<double>("--spin")`). For these parameters, you can optionally pass them to the application via the command line, which will override the value in the parameter file, for example `./emissivity --outfile spin0.dat --spin 0`. This can be useful if you want to run several similar simulations for which you vary just a few of the parameters.

## Citation and acknowledgement
This code is free to use for any research project, however please bear in mind that creation of software such as this involves a significant effort on the part of the authors. Therefore please acknowledge the use of this code in any publications that may result from its use, and cite the relevant publications:
- Wilkins and Fabian 2012, MNRAS 424, 1284-1296
- Wilkins et al. 2016, MNRAS 458, 200
- Wilkins et al. 2020, MNRAS 493, 5532
- 