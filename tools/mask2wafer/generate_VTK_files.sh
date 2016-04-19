#!/bin/bash

viennaTS="../build/viennats"

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

#-- Generate executable
echo "Generate executable: mask2wafer"
g++ -o mask2wafer litho.cpp

#-- Generate mask_volume.vtk and wafer_volume.vtk (for ViennaTS import)
echo "Generate mask_boundary.vtk (volume)"
./mask2wafer ./resist_profile_filled.vtu

#-- Use boolean to create the mask and wafer boundaries
echo "Generate wafer_boundary.vtk and mask_boundary.vtk"
rm -rf output_wafer
mkdir output_wafer
../../build/viennats generate_wafer.txt
mv output_wafer/Interface_0_0.vtk ../../examples/DrLitho/mask_boundary.vtk
mv output_wafer/Interface_0_1.vtk ../../examples/DrLitho/wafer_boundary.vtk
rm -rf output_wafer
