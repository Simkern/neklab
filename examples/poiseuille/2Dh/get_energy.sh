#! /usr/bin/bash
# to be run after nekb poiseuile (or mpi) so that a new logfile is written


if [[ -f energy_2Dh.txt ]]; then echo 'energy_2Dh.txt exists!' && exit ; fi

grep 'energy  ' logfile > energy_2Dh.txt


