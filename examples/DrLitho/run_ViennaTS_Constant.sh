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

#-- Generate par_constant.txt file
echo "Generate par_constant.txt file"
rm -f par_constant.txt
cp par.txt par_constant.txt
echo "geometry_file  =   \"wafer_volume.vtk\";" >> par_constant.txt
echo "output_path = \"./output_constant/\";" >> par_constant.txt
echo "" >> par_constant.txt
echo "surface_geometry = false;" >> par_constant.txt
echo "report_import_errors = false;" >> par_constant.txt
echo "" >> par_constant.txt
echo "" >> par_constant.txt
echo "cfl_condition = 0.1;" >> par_constant.txt
echo "grid_delta = 0.001;" >> par_constant.txt
echo "" >> par_constant.txt
echo "remove_bottom = true;" >> par_constant.txt
echo "" >> par_constant.txt
echo "num_dimensions = 3;" >> par_constant.txt
echo "omp_threads=3;" >> par_constant.txt
echo "" >> par_constant.txt
echo "boundary_conditions = {" >> par_constant.txt
echo "	{PERIODIC,PERIODIC}," >> par_constant.txt
echo "	{PERIODIC,PERIODIC}," >> par_constant.txt
echo "	{INFINITE,INFINITE}" >> par_constant.txt
echo "};" >> par_constant.txt
echo "" >> par_constant.txt
echo "open_boundary=\"+z\";" >> par_constant.txt
echo "default_disk_orientation={0,0,0};" >> par_constant.txt
echo "" >> par_constant.txt
echo "processes = {" >> par_constant.txt
while read -r line
do
   echo $line >> par_constant.txt
done < "generate_LS"
echo "   {" >> par_constant.txt
echo "	print_velocities=true;" >> par_constant.txt
echo "	print_materials=true;" >> par_constant.txt
echo "	process_time= 10;" >> par_constant.txt
echo "	model_name=\"ConstantRates\";" >> par_constant.txt
echo "	parameters= {" >> par_constant.txt
echo "		constant_rates={-0.005,0.};" >> par_constant.txt
echo "	};" >> par_constant.txt
echo "	output_times={0,1,2,3,4,5,6,7,8,9,10};" >> par_constant.txt
echo "   }" >> par_constant.txt
echo "};" >> par_constant.txt

#-- Run a sample constant etch with a mask:
echo "Run a sample constant etch with a mask"
rm -rf output_constant
mkdir output_constant
${viennaTS} par_constant.txt

#-- Output is in the output_constant folder
echo "Output is in the output_constant folder"
echo ""
