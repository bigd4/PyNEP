#!/bin/sh
PYNEP_PATH=$1
LAMMPS_PATH=$2
curdir=$(pwd)
cp $PYNEP_PATH/nep_cpu/src/nep.*  $PYNEP_PATH/inference/lammps/USER-NEP/
cp -r $PYNEP_PATH/inference/lammps/USER-NEP/ $LAMMPS_PATH/src
cd $LAMMPS_PATH/src
make clean-serial
make no-user-nep
make yes-user-nep
make serial
mv lmp_serial $curdir
