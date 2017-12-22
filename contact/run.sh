#!/bin/bash

#time python findContactETH.py ../ethz/90wt_0bar.txt
#time python findContactETH.py ../ethz/90wt_2000bar.txt
#time python findContact.py ../DEM/post/dump104000.generation ../DEM/post/generation.104000.contact_dump 
#time python findContact.py ../DEM/post/dump300000.generation ../DEM/post/generation.300000.contact_dump

time python clustering.py 90wt_0bar_contact.msh     #should see 963 and 1183
#time python clustering.py 90wt_2000bar_contact.msh

#eigval is wrong with test problem, BF wrong with real problem
#like 0 should not be on list, I'm pretty sure 
