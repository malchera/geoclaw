#!/bin/bash 

file="rpn2_geoclaw_soa.f90"
phase="vec,loop"
#file=$1                      
#phase=$2

./generate_opt_report.sh <<EOF
    $file   
    5
    stdout
    $phase 
EOF
