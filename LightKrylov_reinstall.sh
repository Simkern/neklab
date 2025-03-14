#!/bin/bash

cd LightKrylov

############################################
######                                 #####
######     INSTALLING LIGHT KRYLOV     #####
######                                 #####
############################################

# Sets the compilation options depending on the detected compiler.
if command -v mpiifort >/dev/null 2>&1; then
    FPM_FFLAGS="-Ofast -xHost -g -traceback -DMPI"
    FPM_FC="mpiifort"
else
    FPM_FFLAGS="-march=native -O3 -funroll-loops -DMPI"
    FPM_FC="mpifort"
fi

export FPM_CC
export FPM_FFLAGS
export FPM_FC

# Install LightKrylov.
fpm install

# Run the unit tests.
read -p "Do you want to run the tests for LightKrylov? (y/n): " confirm
if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
    fpm test
fi

cd -

# Finalize.
echo "LightKrylov setup complete."
