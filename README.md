ViennaTS
--------------------------

Developer repository for ViennaTS, a C++, OpenMP-parallelized Topography simulator.
ViennaTS is currently in a prototype state.

System requirements
--------------------------

* C++ compiler
* OpenMP
* Boost C++ Libraries
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
CMAKE_BUILD_TYPE   = debug, (release) # Turn off/on optimizations (default: release, i.e., optimized mode)
</pre>

Authors and Contact
------------------------

Current contributors: Lado Filipovic, Paul Manstetten, and Josef Weinbub

Contact us via: viennats AT iue DOT tuwien DOT ac DOT at

Founder and initial developer was Otmar Ertl; not active anymore.

ViennaTS was developed under the aegis of the 'Institute for Microelectronics' at the 'TU Wien'.
http://www.iue.tuwien.ac.at/

License
--------------------------
See file LICENSE in the base directory.
