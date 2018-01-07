#!/bin/bash

time python findContactETH.py --clustering ../ethz/90wt_0bar.txt
time python findContactETH.py --clustering ../ethz/90wt_2000bar.txt
#time python findContact.py ../DEM/post/dump104000.generation ../DEM/post/generation.104000.contact_dump 
#time python findContact.py ../DEM/post/dump300000.generation ../DEM/post/generation.300000.contact_dump

#time python2.7 clustering.py --style=scipy      90wt_0bar_contact.msh &> clustering_90wt_0bar_scipy.txt     #should see 963 and 1183
#time python2.7 clustering.py --style=networkx   90wt_0bar_contact.msh &> clustering_90wt_0bar_networkx.txt     #should see 963 and 1183
#time python2.7 clustering.py --style=eigen      90wt_0bar_contact.msh &> clustering_90wt_0bar_eigen.txt     #should see 963 and 1183
#time python2.7 clustering.py --style=bruteforce 90wt_0bar_contact.msh &> clustering_90wt_0bar_bruteforce.txt     #should see 963 and 1183
#time python2.7 clustering.py 90wt_2000bar_contact.msh

#time python2.7 clustering.py --style=scipy
#time python2.7 clustering.py --style=networkx
#time python2.7 clustering.py --style=eigen
#time python2.7 clustering.py --style=bruteforce

#eigval is wrong with test problem, BF wrong with real problem
#like 0 should not be on list, I'm pretty sure


#test the new styles on the real problem 





#add clustering/disconnect as option to findContact*.py
#check msh file okay
#see if this solves truss problem
#only print contact to FT contact list when connected to main body
