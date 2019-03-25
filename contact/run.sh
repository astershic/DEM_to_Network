#!/bin/bash

time python findContactETH.py --clustering --writeOnlyConnected ../ethz/90wt_0bar.txt
time python findContactETH.py --clustering --writeOnlyConnected ../ethz/90wt_2000bar.txt
#time python findContactDEM.py --clustering --writeOnlyConnected ../DEM/post/dump104000.generation ../DEM/post/generation.104000.contact_dump 
#time python findContactDEM.py --clustering --writeOnlyConnected ../DEM/post/dump300000.generation ../DEM/post/generation.300000.contact_dump

#time python2.7 clustering.py --style=scipy      90wt_0bar_contact.msh &> clustering_90wt_0bar_scipy.txt     #should see 963 and 1183
#time python2.7 clustering.py --style=networkx   90wt_0bar_contact.msh &> clustering_90wt_0bar_networkx.txt     #should see 963 and 1183
#time python2.7 clustering.py --style=eigen      90wt_0bar_contact.msh &> clustering_90wt_0bar_eigen.txt     #should see 963 and 1183  #this takes by far longest
#time python2.7 clustering.py --style=bruteforce 90wt_0bar_contact.msh &> clustering_90wt_0bar_bruteforce.txt     #should see 963 and 1183
#time python2.7 clustering.py 90wt_2000bar_contact.msh

#time python2.7 clustering.py --style=scipy
#time python2.7 clustering.py --style=networkx
#time python2.7 clustering.py --style=eigen
#time python2.7 clustering.py --style=bruteforce

#first check planes, make sure enough bcs being made
#then check clustering; only include (multiply????-single may not contribute) connected elements
#check msh file okay
#see if this solves truss problem
