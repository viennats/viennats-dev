ViennaTS
--------------------------

Developer repository for ViennaTS, a C++, OpenMP-parallelized Topography simulator.
Releases are tagged on the master branch and available [in the releases section](https://github.com/viennats/viennats-dev/releases).

For more information and help to getting started, visit [viennats.github.io](https://viennats.github.io/).

System requirements
--------------------------

* C++ compiler
* OpenMP
* Boost C++ Libraries ([boost.org](https://www.boost.org/))
* VTK C++ Libraries ([vtk.org](https://vtk.org/))
* HDF5 (optional - required for TDR file support)

Currently supported operating systems
--------------------------
* GNU/Linux (32/64Bit)

Building instructions
--------------------------

To build ViennaTS, clone the repository and issue the following suggested commands:

<pre>
$> cd viennats-dev    # the checked-out GIT folder
$> mkdir build        # the build folder
</pre>

Configure the build, default build type is the 'optimized/release' mode:
<pre>
$> cd build/
$> cmake ..
</pre>
Watch for Warning/Error messages during the configuration stage.

Now build the 'viennats' simulation executable
<pre>
$> make
</pre>

CMake Options
--------------------------

<pre>
CMAKE_BUILD_TYPE   = Debug, (Release) # Turn off/on optimizations (default: Release, i.e., optimized mode)
VIENNATS_STATIC_BUILD = ON/OFF # build ViennaTS without dynamically linked libraries (VTK must be built statically when this is ON)
</pre>

Authors and Contact
------------------------

Current contributors: Lado Filipovic, Paul Manstetten, Xaver Klemenschits and Josef Weinbub

Contact us via: viennats@iue.tuwien.ac.at

Founder and initial developer was Otmar Ertl; not active anymore.

ViennaTS was developed under the aegis of the 'Institute for Microelectronics' at the 'TU Wien'.
http://www.iue.tuwien.ac.at/

License
--------------------------
See file LICENSE in the base directory.
