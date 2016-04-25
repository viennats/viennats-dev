#!/bin/bash

echo ""
echo ""
echo "#==================================================#"
echo "#                                                  #"
echo "#     This test is meant to enable ViennaTS to     #"
echo "#     read a Dr.Litho mask and perform etch        #"
echo "#     simulations with the generated mask          #"
echo "#                                                  #"
echo "#==================================================#"
echo ""
echo ""

viennaTS="../../build/viennats"

#-- Generate par_plasma.txt file
echo "Generate par_plasma.txt file"
rm -f par_plasma.txt
cp par.txt par_plasma.txt
echo "geometry_file  =   \"wafer_volume.vtk\";" >> par_plasma.txt
echo "output_path = \"./output_plasma/\";" >> par_plasma.txt
echo "" >> par_plasma.txt
echo "surface_geometry = false;" >> par_plasma.txt
echo "report_import_errors = false;" >> par_plasma.txt
echo "" >> par_plasma.txt
echo "" >> par_plasma.txt
echo "cfl_condition = 0.1;" >> par_plasma.txt
echo "" >> par_plasma.txt
echo "input_scale = 1e-4;" >> par_plasma.txt
echo "grid_delta  = 2.5e-7;" >> par_plasma.txt
echo "" >> par_plasma.txt
echo "remove_bottom = true;" >> par_plasma.txt
echo "" >> par_plasma.txt
echo "num_dimensions = 3;" >> par_plasma.txt
echo "omp_threads=3;" >> par_plasma.txt
echo "" >> par_plasma.txt
echo "boundary_conditions = {" >> par_plasma.txt
echo "	{PERIODIC,PERIODIC}," >> par_plasma.txt
echo "	{PERIODIC,PERIODIC}," >> par_plasma.txt
echo "	{INFINITE,INFINITE}" >> par_plasma.txt
echo "};" >> par_plasma.txt
echo "" >> par_plasma.txt
echo "open_boundary=\"+z\";" >> par_plasma.txt
echo "default_disk_orientation={0,0,0};" >> par_plasma.txt
echo "" >> par_plasma.txt
echo "processes = {" >> par_plasma.txt
IFS=''
while read -r line
do
   echo $line >> par_plasma.txt
done < "generate_LS"
echo "   {" >> par_plasma.txt
echo "	print_coverages=true;" >> par_plasma.txt
echo "	print_rates=true;" >> par_plasma.txt
echo "	print_velocities=true;" >> par_plasma.txt
echo "	print_materials=true;" >> par_plasma.txt
echo "" >> par_plasma.txt
echo "	output_times = {0};" >> par_plasma.txt
echo "	output_times_period_length = 0.05;" >> par_plasma.txt
echo "	output_times_periodicity=20;" >> par_plasma.txt
echo "	final_output=true;" >> par_plasma.txt
echo "" >> par_plasma.txt
echo "	process_time= 1;" >> par_plasma.txt
echo "	model_name=\"SF6_O2PlasmaEtching\";" >> par_plasma.txt
echo "	add_layer=0;" >> par_plasma.txt
echo "	parameters= {" >> par_plasma.txt
echo "		direction={0,0,-1};" >> par_plasma.txt
echo "		statistical_accuracy=500;" >> par_plasma.txt
echo "		min_ion_energy=100;" >> par_plasma.txt
echo "		delta_ion_energy=40;" >> par_plasma.txt
echo "		flux_ion=1e16;" >> par_plasma.txt
echo "		flux_oxygen=3e17;" >> par_plasma.txt
echo "		flux_fluorine=5.5e18;" >> par_plasma.txt
echo "		a_oxygen=1.;" >> par_plasma.txt
echo "	};" >> par_plasma.txt
echo "	start_iteration_cycles=10;" >> par_plasma.txt
echo "	iteration_cycles=0;" >> par_plasma.txt
echo "   }" >> par_plasma.txt
echo "};" >> par_plasma.txt

#-- Run a sample SF6/O2 plasma etch:
echo "Run a sample SF6/O2 plasma etch"
rm -rf output_plasma
mkdir output_plasma
$viennaTS par_plasma.txt

#-- Output is in the output_plasma folder
echo "Output is in the output_plasma folder"
echo ""
