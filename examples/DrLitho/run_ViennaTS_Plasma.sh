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

#-- Run a sample SF6/O2 Plasma Etch:
rm -rf output_plasma
mkdir output_plasma
$viennaTS par_plasma.txt

#-- output is in the output_plasma folder
