#!/bin/bash

#exa -lhFD --icons
#echo ""
#read -p "Enter input folder: " r1
#read -p "Enter output folder: " r2

r1="derivatives"
r2="data_packages"

find $r1/* -type f \( -name "maMc*sub-15*" \) -exec zip -9 $r2/derivatives_dcm_ready_without_gainmat_single_subject.zip {} +
find $r1/* -type f \( -name "maMc*sub-15*" -or -name "SPMgain*sub-01*" \) -exec zip -9 $r2/derivatives_dcm_ready_with_gainmat_single_subject.zip {} +
find $r1/* -type f \( -name "maMc*sub-15*" -or -name "SPMgain*sub-01*" \) -exec zip -9 -j $r2/flattened_dcm_ready_with_gainmat_single_subject.zip {} +

find $r1/* -type f \( -name "maMc*" \) -exec zip -9 $r2/derivatives_dcm_ready_without_gainmat_all_subjects.zip {} +
find $r1/* -type f \( -name "maMc*" -or -name "SPMgain*" \) -exec zip -9 $r2/derivatives_dcm_ready_with_gainmat_all_subjects.zip {} +
find $r1/* -type f \( -name "maMc*" -or -name "SPMgain*" \) -exec zip -9 -j $r2/flattened_dcm_ready_with_gainmat_all_subjects.zip {} +

find $r1/* -type f \( -name "*aMc*" -or -name "SPMgain*" \) -exec tar -rvf $r2/derivatives_processed_data_with_gainmat_all_subjects.tar {} \;
xz -v -9 -z --threads=0 $r2/derivatives_processed_data_with_gainmat_all_subjects.tar
