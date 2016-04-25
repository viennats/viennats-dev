#!/bin/bash

mask=${1}

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
echo ""
g++ -o mask2wafer litho.cpp

#-- Generate wafer_volume.vtk (for ViennaTS import)
echo "Generate wafer_volume.vtk (volume)"
echo ""
./mask2wafer ${mask}

#-- Generate ViennaTS parameters file
echo "Generate ViennaTS parameters file (generate_LS)"
echo ""
rm -f generate_LS
echo "   {" >> generate_LS
echo "	model_name=\"BooleanOperation\";" >> generate_LS
echo "	parameters= {" >> generate_LS
echo "		geometry_file=\"${mask}\";" >> generate_LS
echo "		surface_geometry=true;" >> generate_LS
echo "		remove_bottom=false;" >> generate_LS
echo "		level=-1;" >> generate_LS
echo "	};" >> generate_LS
echo "   }," >> generate_LS
echo "   {" >> generate_LS
echo "	model_name=\"BooleanOperation\";" >> generate_LS
echo "	parameters= {" >> generate_LS
echo "		geometry_file=\"${mask}\";" >> generate_LS
echo "		surface_geometry=true;" >> generate_LS
echo "		remove_bottom=true;" >> generate_LS
echo "		level=+1;" >> generate_LS
echo "	};" >> generate_LS
echo "   }," >> generate_LS
echo "   {" >> generate_LS
echo "	model_name=\"Mask\";" >> generate_LS
echo "	parameters= {" >> generate_LS
echo "		mask_file=\"resist_profile.vtk\";" >> generate_LS
echo "		surface_geometry=true;" >> generate_LS
echo "		remove_bottom=false;" >> generate_LS
echo "	};" >> generate_LS
echo "   }," >> generate_LS

#-- Move the wafer_volume.vtk, ${mask}.vtk, and generate_LS
echo "Move the wafer_volume.vtk, ${mask}, and generate_LS to examples"
echo ""
cp wafer_volume.vtk ../../examples/DrLitho/wafer_volume.vtk
cp ${mask} ../../examples/DrLitho/${mask}
cp generate_LS ../../examples/DrLitho/generate_LS
