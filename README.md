ViennaTS
--------------------------

Developer repository for ViennaTS, a C++, OpenMP-parallelized Topography simulator.
ViennaTS is currently in a prototype state.

System requirements
--------------------------

* C++ compiler
* OpenMP
* Boost C++ Libraries
* SPRNG 

Currently supported operating systems
--------------------------
* GNU/Linux

Building instructions
--------------------------

To build ViennaTS, clone the repository and issue the following suggested commands:

<pre>
$> cd viennats-dev  # the checked-out GIT folder
$> mkdir build        # the build folder
</pre>

Configure the build, default build type is the 'optimized/release' mode:
<pre>
$> cd build/
$> cmake ..
</pre>

Now build the 'viennats' executable 
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

Founder and initial author was Otmar Ertl; not active anymore.

ViennaTS was developed under the aegis of the 'Institute for Microelectronics' at the 'TU Wien'.


License
--------------------------
See file LICENSE in the base directory.
