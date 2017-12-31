#!/sw/bin/python

#trying to find main contact cluster in mesh, so that boundary conditions are sufficient - no free bodies

#python clustering.py 90wt_0bar_contact.msh
import getopt
import sys
import math
import numpy as np
from scipy.sparse import *
from scipy.sparse.linalg import *

def build_adjacency_matrix(elements):

    #get list of nodes
    nodes = set()
    for elem in elements:
        for node in elem:
           nodes.add(node)

#    print nodes
    nodes = list(nodes) 
#    print nodes

    #nodes-index mapping: key = node#, value = index
    nodeToIndex = {}
    for i in range(0, len(nodes)):
        node = nodes[i]
        nodeToIndex[node] = i 

    #build sparse adjacency matrix
    n = len(nodes)
    A = lil_matrix( (n,n), dtype=np.int8 )
    for elem in elements:
        n1 = elem[0]
        n2 = elem[1]
        i1 = nodeToIndex[n1]
        i2 = nodeToIndex[n2]
        A[i1, i2] = 1
        A[i2, i1] = 1

    #check that every row has an entry
    Asum = A.sum() 
    assert( np.product(Asum) != 0.0 )

    #convert to CSR matrix to make faster
    A = csr_matrix(A)
#    print A
#    print A.toarray()   

    return A, nodes, nodeToIndex

def calc_clusters(elements):


    #adjacency matrix
    print " #Building Adjacency Matrix"
    A, nodes, nodeToIndex = build_adjacency_matrix(elements)
    print "  *Matrix Shape: ", A.shape
    print "  *Non-zero: ", A.nnz

    #diagonal matrix
    print " #Building Diagonal Matrix"
    n = len(nodes)
    Adiag = sum(A).toarray()[0]
#    print A.toarray()
    D = diags( Adiag, 0, format="csr", dtype=np.int8)
#    print D.toarray()

    #laplacian
    print " #Building Laplacian Matrix"
    L = D - A
    L = L.astype(float)
#    print "L = ",L.toarray()

    #normalized symmetric laplacian
    print " #Building Symmetric Laplacian Matrix"
    sqrtD = np.sqrt(Adiag)
    sqrtD = diags( sqrtD, 0, format="csc")
    invSqrtD = inv(sqrtD)
    Lnorm = eye(n) - invSqrtD*A*invSqrtD    
#    print "Lnorm = ",Lnorm.toarray()

    #compute eigenvalues 
    print " #Computing Eigenvalues"
#    print "n matrix size = ", n
    laplacian = L 
#    laplacian = Lnorm 
    vals,vecs = eigsh( laplacian.asfptype() , k = n/2, sigma=-1)   #n/2 values:   n/2 vecs, length n
#    print "vals = ",vals
#    print "vecs = ",vecs
    tol = 1e-3
    vals[np.abs(vals) < tol] = 0
#    print "vals = ",vals

    #normalize (Lnorm)
#    vecs[np.abs(vecs) < tol] = 0
#    vecs[np.abs(vecs) >= tol] = 1

    #normalize (L)
    for i in range(0,  len(vecs[0,:]) ):
        vec = vecs[:,i]
        vec = vec/np.max(np.abs(vec))
        vec[np.abs(vec) < 0.9] = 0
        vec[np.abs(vec) >= 0.9] = 1
        vecs[:,i] = vec
#    print "normalized vectors = ", vecs

    #find connected groups
    print " #Finding Connected Groups"
    nGroups = list(vals).count(0)
    groups = []
    groupSize = []
    #iterate over columns
    for i in range(0, nGroups):
        vec = vecs[:,i]
        nonZero = vec.nonzero()[0]   #convert from indices (only nodes that are in elements) to real node numbers
        nonZeroNodes = [ nodes[x] for x in nonZero ]
        #group = indices of non-zero entries
        groups.append( nonZeroNodes ) 
        groupSize.append( len(nonZeroNodes) )
    print groups
    print groupSize
    assert( len(groups) == len(vals[np.abs(vals) < tol]) )
    return nodes,groups,groupSize,nodeToIndex


def readMshFile(filename):
    elements = []
    nElem = 0
    with open(filename, 'r') as inFile:
        elemSection = 0
        for line in inFile:
            if (elemSection == 0 and '$Elements'.lower() in line.lower()):
                print 'found elements'
                elemSection = 1
            elif (elemSection == 2 and '$EndElements'.lower() in line.lower()):
                print 'found end elements'
                elemSection = 0
            elif (elemSection == 1):
                nElem = int(line.strip())
                print "# elements = ",nElem
                elemSection = 2
            elif (elemSection == 2):
                ls = line.strip()
                ls = ls.split()

                #only interested in "real" elements
                if ( int(ls[1]) == 1):
                    n1 = int(ls[3])
                    n2 = int(ls[4])
                    elements.append([n1, n2])

#(real element) element number, type: 1 = 2-node line, number of physical group tags = 0, nodeID1, nodeID2 
#(physical tag) "element" number, type 15=point, number of physical IDs, physical ID, nodeID
            else:
                pass
            
    return elements

def identifyDisconnected(nodes, groups, groupSize, nodeToIndex):
   
    #find largest group 
    groupIdx = 0
    maxGroupSize = groupSize[0]
    for i in range(0, len(groupSize)):
        gs = groupSize[i]
        if (gs > maxGroupSize):
            maxGroupSize = gs
            groupIdx = i
#    print "Largest group = ",groupIdx, "with size",maxGroupSize,"nodes members",groups[groupIdx]

    #make set of nodes not in largest group
    disconnected = set()
    for i in range(0, len(groups)):

        if (i == groupIdx):
            continue

        group = groups[i]
        for node in group:
#            print node, type(node)
            disconnected.add( node ) 
    print "disconnected: ", disconnected

    return disconnected

def calc_clusters_bruteForce(elements):

    A, nodes, nodeToIndex = build_adjacency_matrix(elements)

    includedInGroup = [-1] * len(nodes)
    groupIdx = -1
    nodeIdx = -1
    groups = []
    #while not all included in groups
    while any( x == -1 for x in includedInGroup ):
        nodeIdx += 1
        node = nodes[nodeIdx]
        index = nodeToIndex[node]
#        print "group ",groupIdx+1
#        print "checking node ", nodeIdx
#        print "in group? ",includedInGroup[nodeIdx]
        if (includedInGroup[index] > -1):
            continue
        #this is an unincluded node, add to group and all its neighbors 
        #make new group
        thisGroup = set()
        groupIdx += 1
        includedInGroup[index] = groupIdx
        thisGroup.add(nodeIdx)
        #search through all elements, for all nodes connected to this group...
        #careful with node# vs index #
        changed = 1
        while (changed == 1):
#            print "looking at elements"
            changed = 0
            for elem in elements:
                n1 = elem[0]
                n2 = elem[1]

                idx1 = nodeToIndex[n1]
                idx2 = nodeToIndex[n2]

                if (n1 in thisGroup and includedInGroup[idx2] == -1):
                    thisGroup.add(n2)
                    includedInGroup[idx2] = groupIdx
                    changed = 1
#                    print "added node",n2
                if (n2 in thisGroup and includedInGroup[idx1] == -1):
                    thisGroup.add(n1)
                    includedInGroup[idx1] = groupIdx
                    changed = 1
#                    print "added node",n1
#            print "changed = ",changed
#            print "   *this group = ",thisGroup

        groups.append(thisGroup)

    groupSize = []
    for group in groups:
        groupSize.append(len(group))

    print groups
#    print includedInGroup
    print groupSize
    return nodes, groups, groupSize, nodeToIndex

def cluster_main(filename,style):

    if (filename == None):
        elements = [[0,1], [1,2], [2, 4], [4, 0], [4, 1], [5, 6], [3,5], [3,6], [3,7] , [8,9] , [0,10], [11,11], [0,15], [13, 15] ]
    else:
        elements = readMshFile(filename)
 
    print "Filename = ",filename

    if   (style.lower() == 'eigen'):
        nodes,groups,groupSize,nodeToIndex = calc_clusters_eigen(elements)
    elif (style.lower() == 'networkx'):
        nodes,groups,groupSize,nodeToIndex = calc_clusters_networkx(elements)
    elif (style.lower() == 'scipy'):
        nodes,groups,groupSize,nodeToIndex = calc_clusters_scipy(elements)
    elif (style.lower() == 'bruteforce'):
        nodes,groups,groupSize,nodeToIndex = calc_clusters_bruteForce(elements)
    else:
        print "Style <{}> not implemented".format(style)
        assert(1 == 0)

    disconnected = identifyDisconnected(nodes, groups, groupSize, nodeToIndex)

if __name__ == "__main__":
    optlist,args = getopt.getopt(sys.argv[1:],'',longopts=['style='])
    style = 'scipy' #default
    for item in optlist:
        if (item[0] == '--style'):
            style = item[1] 
    print "style = ",style
    if (len(args) > 0):
        for filename in args:
            cluster_main(filename,style)
    else:
        cluster_main(None,style)
