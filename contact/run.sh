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


#ffmpeg -r 1 -i data/input-%4d.png -pix_fmt yuv420p -r 30 data/output.mp4

#contact status (generation*.dump): id p1x p1y p1z p2x p2y p2z p1id p2id periodic fx fy fz contactarea delta

#######
#TO DO:
#create imagery for srdjan
#(PY) contact boundary Z override 
#(PY) voronoi boundary Z override 
#(PY) connection to fabric tensor calculator
#(PY) truss network to <voltage,thermal,mechanical> moose simulation 
#(PY) voronoi network to <voltage,thermal,mechanical> moose simulation 
#python:
#  remove particles above certain height
#  run until kinetic under limit
#  compression simulations
#ryan's insitu placement
#######
#DONE:
#(PY) contact boundary conditions (physical groups)
#(PY) voronoi boundary conditions (physical groups)
#(PY) voronoi connection to liggghts data (voronoi.py)
#(PY) voronoi connection to ethz data (voronoi.py)
#(PY) contact dump to truss network (findContact.py)
#(PY) element dataa (delta, contact area, contact force)
