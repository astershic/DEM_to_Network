#!/bin/bash

find post/* -name "*generat*" -print0 | xargs -0 rm
rm log.liggghts
rm post/*.png

source ~/.bash_profile
Liggghts
np=10
#mpirun -np $np lmp_openmpi -e both -v r 0.0045 -v insertSteps 120000 -in in.generate  #coarse
#mpirun -np $np lmp_openmpi -e both -v r 0.003 -v insertSteps 180000 -in in.generate   #mid
#mpirun -np $np lmp_openmpi -e both -v r 0.0015 -v insertSteps 600000 -in in.generate  #fine
mpirun -np $np lmp_openmpi -e both -v r 0.003 -v insertSteps 200000 -in in.generate-distribution   #distribution
cd post
lpp dump*generation
cd ..
