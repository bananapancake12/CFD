#!/bin/bash

# Check submissions from all candidates for plagiarism
# Written by James Taylor 
# November 2023

# Loop over all subdirectories and uncompress archives
find . -name '*.tar.gz' -execdir tar -xzf '{}' \;

# Check the archive has been uncompressed successfully
for dir in ./*
do
    path=$dir/SaveSrc	
    if [ -d "$path" ]; then
        echo "$path exists."
    else
        echo "$path does not exist."
    fi
done

# Check for presence of example code in the makefile
makefiles=( ./*/SaveSrc/makefile )
for ((i = 0; i < ${#makefiles[@]} ; i++))
do
    if grep -q _jvt "${makefiles[i]}"; then
        echo "${makefiles[i]} contains example subroutines"
    fi
done
        
# Check all candidates subroutines against each other
codefiles=( "apply_bconds.f90" "calc_areas.f90" "check_mesh.f90" "euler_iteration.f90" "flow_guess.f90" "generate_mesh.f90" "plot_contours.py" "routines.py" "read_settings.f90" "set_secondary.f90" "set_timestep.f90" "flux_stencil.f90" "smooth_stencil.f90") 
for codefile in ${codefiles[@]}
do
    # Create list of all candidates paths to the same file from codefiles
    submitfiles=( ./*/SaveSrc/$codefile )
    echo Checking $codefile

    # Loop over all candidates in pairs
    for ((i = 0; i < ${#submitfiles[@]} ; i++))
    do
        for ((j = i+1; j < ${#submitfiles[@]}; j++))
        do

            # Check for matches, no point checking a candidate against themselves
            if [ $i -ne $j ]; then
                if diff -q "${submitfiles[i]}" "${submitfiles[j]}" > /dev/null
#                if diff -q "${submitfiles[i]}" "${submitfiles[j]}"
                then
                    echo "Files ${submitfiles[i]#*/} and ${submitfiles[j]#*/} are identical"
                fi
            fi
        done
    done
done

