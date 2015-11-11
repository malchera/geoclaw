#!/bin/bash 

#file="flux2fw.f"            
file=$1                      
phase=$2

./generate_opt_report.sh <<EOF
    $file   
    5
    stdout
    $phase 
     
EOF
