#!/bin/bash

viennaTS="../../build/viennats"

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

#-- Run a sample Constant Etch with a mask:
rm -rf output_constant
mkdir output_constant
$viennaTS par_constant.txt

#-- output is in the output_constant folder
