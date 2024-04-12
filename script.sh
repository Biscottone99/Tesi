#!/bin/bash
ifort basis.f90 -o basis.e
./basis.e
ifort input.f90 -o input.e
./input.e

ifort module.f90 elix.f90 -o geometria.e  -qopenmp -traceback
ifort module.f90 flash.f90 -o flash.e -qmkl  -qopenmp
for i in $(seq 1 201); do
    cp "input_${i}.dat" "input.dat"
    ./geometria.e
    ./flash.e
    cp "results.dat" "results_${i}.dat"
    cp "fort.100" "GS_${i}.dat"
    cp "fort.102" "S1_${i}.dat"
    cp "post.dat" "post_${i}.dat"
    cp "post1.dat" "post1_${i}.dat"
    cp "mcrot_r.dat" "mcrot_r_${i}.dat"
    cp "mcrot_i.dat" "mcrot_i_${i}.dat"
    cp "bcrot_r.dat" "bcrot_r_${i}.dat"
    cp "bcrot_i.dat" "bcrot_i_${i}.dat"
  
done
ifort post.f90 -o post.e
./post.e

ifort post2.f90 -o post2.e
./post2.e

ifort axis.f90 -o axis.e
./axis.e
python3 axis.py
echo "Script completato."
