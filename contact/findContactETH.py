#!/sw/bin/python

#general:
#ITEM: ATOMS id type type x y z ix iy iz vx vy vz fx fy fz omegax omegay omegaz radius c_contactnum

#histogram:   --> cellDatatoPointData --> Histogram
#line thickness: ---> cellDatatoPointData ---> ExtractSurface --> Tube
import os
import math
import sys
from clustering import calc_clusters,identifyDisconnected

#run: python sphereContact.py *.txt

def calcDistance(cI, cJ):
    xI = cI[0]
    yI = cI[1]
    zI = cI[2]
    xJ = cJ[0]
    yJ = cJ[1]
    zJ = cJ[2]

    dist = math.sqrt( (xI-xJ)*(xI-xJ) + (yI-yJ)*(yI-yJ) + (zI-zJ)*(zI-zJ) )

    return dist

def findPhysical(ids, rads, centroids, limits):
    physical = {}

    xmin = limits[0][0]
    xmax = limits[0][1]
    ymin = limits[1][0]
    ymax = limits[1][1]
    zmin = limits[2][0]
    zmax = limits[2][1]

    #might consider a percentage-based BC
    ratio = 0.90
    xrange = xmax - xmin
    xmean = 0.5*(xmin + xmax)
    xmin = xmean - 0.5*ratio*xrange
    xmax = xmean + 0.5*ratio*xrange
    yrange = ymax - ymin
    ymean = 0.5*(ymin + ymax)
    ymin = ymean - 0.5*ratio*yrange
    ymax = ymean + 0.5*ratio*yrange
    zrange = zmax - zmin
    zmean = 0.5*(zmin + zmax)
    zmin = zmean - 0.5*ratio*zrange
    zmax = zmean + 0.5*ratio*zrange

    for i in range(0, len(ids)):
        idx = ids[i]

        xc = centroids[i][0]
        yc = centroids[i][1]
        zc = centroids[i][2]
        rad = rads[i]

        bc = []
        if (xc - rad <= xmin):
            bc.append(1) 
        if (xc + rad >= xmax):
            bc.append(2) 
        if (yc - rad <= ymin):
            bc.append(3) 
        if (yc + rad >= ymax):
            bc.append(4) 
        if (zc - rad <= zmin):
            bc.append(5) 
        if (zc + rad >= zmax):
            bc.append(6)
        physical[idx] = bc 

    return physical

def getPhysicalLength(physical):
    physLength = 0

    for item in physical:
        physLength += len(item)

    return physLength

def findContacts(ids, rads, centroids):
    idpairs = []
    unitvecs = []
    areas = []
    deltas = []
    radpairs = []
    nodesInElements = set()

    for i in range(0, len(ids)-1):
        idxI = ids[i]
        rI = rads[i]

        for j in range(i+1, len(ids)):
            idxJ = ids[j]
            rJ = rads[j]
             
            dist = calcDistance(centroids[i], centroids[j])
    
            if (rI + rJ >= dist):
                xI = centroids[i][0]
                yI = centroids[i][1]
                zI = centroids[i][2]
                xJ = centroids[j][0]
                yJ = centroids[j][1]
                zJ = centroids[j][2]
 
                vec = [xI-xJ, yI-yJ, zI-zJ]
                norm = math.sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] )
                unitvec = [vec[0]/norm, vec[1]/norm, vec[2]/norm]
                delta = rI + rJ - dist
                R = 1.0 / (1.0/rI + 1.0/rJ)
                a = math.sqrt(R * delta)
                area = math.pi * a * a

                idpairs.append([idxI, idxJ])
                unitvecs.append(unitvec)
                areas.append(area)
                deltas.append(delta)
                radpairs.append([rI, rJ])        

                nodesInElements.add(idxI) 
                nodesInElements.add(idxJ) 

    return idpairs,unitvecs,areas,deltas,radpairs,nodesInElements

def readFile(file):
    ids = []   #1-indexed
    rads = []
    centroids = []
    lineNum = 0

    xmin = None
    ymin = None
    zmin = None
    xmax = None
    ymax = None
    zmax = None

    with open(file,'r') as f:
        for line in f:
            lineNum += 1
#            print lineNum,line 
            if lineNum == 1:
                continue

            ls = line.strip().split(',')
            id = int(ls[0])
            vol = float(ls[1])
            xc = float(ls[8])
            yc = float(ls[9])
            zc = float(ls[10])

            xneg = float(ls[2])
            yneg = float(ls[3])
            zneg = float(ls[4])
            xpos = float(ls[5])
            ypos = float(ls[6])
            zpos = float(ls[7])

            if (xneg < xmin or xmin == None):
                xmin = xneg 
            if (yneg < ymin or ymin == None):
                ymin = yneg 
            if (zneg < zmin or zmin == None):
                zmin = zneg 
            if (xpos > xmax or xmax == None):
                xmax = xpos 
            if (ypos > ymax or ymax == None):
                ymax = ypos 
            if (zpos > zmax or zmax == None):
                zmax = zpos 
            
            rad = math.pow(vol*3.0/4.0, 1.0/3.0)
            
            ids.append(id)
            rads.append(rad)
            centroids.append([xc, yc, zc])
            limits = [ [xmin, xmax], [ymin, ymax], [zmin, zmax] ]
    return ids, rads, centroids, limits

def writePhysicalNodes(file, elemCount, verts, physical, nodesInElements):

    extraElem = 0

    for i in range(0, len(verts)):
        key = verts[i]

        #if node is not in an element, don't write it
        if (key not in nodesInElements):
            continue


        bc = physical[key]
        if (len(bc) == 0):
            continue
#        bcString = " ".join(str(x) for x in bc)
#        file.write('{} {} {} {} {}\n'.format(elemCount + extraElem, 15, len(bc), bcString, key) )   #"element" number, type 15=point, number of physical tags, physical ID, nodeID
        for item in bc:
            extraElem += 1
            file.write('{} {} {} {} {}\n'.format(elemCount + extraElem, 15, 1, item, key) )   #"element" number, type 15=point, physical ID, nodeID   ----- note limit of one physical tag per element!!!
    return 
 
def writeElementData(file, field, fieldname):
    #elemData as list, in order of ids:(id1,id2) vector
 
    file.write('$ElementData\n')
    file.write('{}\n'.format(1))    #one string tag
    file.write('"{}"\n'.format(fieldname))    #  the name of the view
    file.write('{}\n'.format(1))    #one real tag
    file.write('{}\n'.format(0.0))    #  the time value
    file.write('{}\n'.format(3))    #three integer tags
    file.write('{}\n'.format(0))    #  the timestep
    file.write('{}\n'.format(1))    #  number of components
    file.write('{}\n'.format(len(field)) )    #  count
    for i in range(0, len(field)):
        file.write('{} {}\n'.format(i+1, field[i]) )   #test data
    file.write('$EndElementData\n')
    return

def writePhysicalNames(file):
    file.write('$PhysicalNames\n')
    file.write('7\n')
    file.write('1 1 xNeg\n')
    file.write('2 2 xPos\n')
    file.write('3 3 yNeg\n')
    file.write('4 4 yPos\n')
    file.write('5 5 zNeg\n')
    file.write('6 6 zPos\n')
    file.write('7 7 disconnected\n')
    file.write('$EndPhysicalNames\n')
 
def writeNodeData(file, nodeIds, field, fieldname, nodesInElements):
    #nodeData as dict
#    print field   
#    print 'nodes: ',nodeIds    #contacts: [id1, id2]
    
    data0 = field[0]  #just to get data length
    if (isinstance(data0, list)):
        lenData = len(data0)
    else:
        lenData = 1

    file.write('$NodeData\n')
    file.write('{}\n'.format(1))    #one string tag
    file.write('"{}"\n'.format(fieldname))    #  the name of the view
    file.write('{}\n'.format(1))    #one real tag (time)
    file.write('{}\n'.format(0.0))    #  the time value
    file.write('{}\n'.format(3))    #three integer tags (timestep index, number of components, count
    file.write('{}\n'.format(0))    #  the timestep
    file.write('{}\n'.format(lenData))    #  number of components
    file.write('{}\n'.format(len(nodesInElements)))    #  count

    for i in range(0, len(nodeIds)):
        key = nodeIds[i]

        #if node not in any element, don't write
        if (key not in nodesInElements):
            continue

        val = field[i]
        if (lenData > 1):
            valstr = ' '.join([str(x) for x in val])
        else:
            valstr = str(val)
        file.write('{} {}\n'.format(key, valstr) )   #test data
    file.write('$EndNodeData\n')
    return
 
def writeMshFile(filename, verts, centroid, edges, physical, elemFields, elemData, nodeFields, nodeData, nodesInElements):

    physicalGroups = 1 

    print '# nodes = {}  ({} in elements)'.format(len(verts), len(nodesInElements))
    print '# elems = ',len(edges)

    if (filename == None):
        return
    with open(filename, 'w') as file:
        file.write('$MeshFormat\n')
        file.write('2.2 0 8\n')
        file.write('$EndMeshFormat\n')
        file.write('$Nodes\n')
        file.write('{}\n'.format(len(nodesInElements)))
        for i in range(0, len(verts)):
            vert = verts[i]
            if (vert not in nodesInElements):
                #if node is not in an element, don't write it.
                continue
            coord = centroid[i]
            file.write('{} {} {} {}\n'.format(vert, coord[0], coord[1], coord[2])) 
        file.write('$EndNodes\n')

        totalElements = len(edges)
        if (physicalGroups == 1):
#            totalElements += len(physical)
            for i in range(0, len(verts)):
                vert = verts[i]
                #only count if node is in an element
                if (vert not in nodesInElements):
                    continue
                totalElements += len(physical[vert])
            print '# physBC = ',totalElements-len(edges)
        file.write('$Elements\n')
        file.write('{}\n'.format(totalElements))
        for i in range(0, len(edges)):
            edge = edges[i]
            file.write('{} {} {} {} {}\n'.format(i+1, 1, 0, edge[0], edge[1])) #type: 1 = 2-node line, number of physical group tags = 0 
        if (physicalGroups == 1):
            writePhysicalNodes(file, len(edges), verts, physical, nodesInElements)
        file.write('$EndElements\n')

        writePhysicalNames(file)

        assert( len(elemFields) == len(elemData) )
        for i in range(0, len(elemFields)):
            writeElementData(file, elemData[i], elemFields[i]) 

        assert( len(nodeFields) == len(nodeData) )
        for i in range(0, len(nodeFields)):
            writeNodeData(file, verts, nodeData[i], nodeFields[i], nodesInElements) 
    return 

def applyDisconnected( disconnected, physical, code=7):
    newPhysical = dict(physical)
    
    for dxID in disconnected:
        if dxID in newPhysical:
            physIDs = newPhysical[dxID]
            physIDs.append(code)
            newPhysical[dxID] = physIDs
            print "physical node",dxID,'disconnected:',physIDs
        else:
            print " new node",dxID,'disconnected:',physIDs
            newPhysical[dxID] = [code]   
 
    return newPhysical

def main(files,tag,clustering=1):

    for file in files:
        base = os.path.basename(file)
        outFile = '.'.join(base.split('.')[:-1])+tag+'_contact.msh'
 
        print "*Reading ETHZ particle stats file ",file
        ids,rad,centroid,limits = readFile(file)
        physical = findPhysical(ids, rad, centroid, limits)

        print "*Finding contacts"
        idpair,unitvec,area,delta,radpair,nodesInElements = findContacts(ids, rad, centroid)

        if (clustering == 1):
            print "*Clustering network"
            nodes,groups,groupSize,nodeToIndex = calc_clusters(idpair)
            disconnected = identifyDisconnected(nodes, groups, groupSize, nodeToIndex)
            physical = applyDisconnected( disconnected, physical)

        print "*Writing msh file: ",outFile
        elemFields = ['delta','area']
        elemData = [area, delta]
        nodeFields = ['position','radius']
        nodeData = [centroid,rad]
        writeMshFile(outFile, ids, centroid, idpair, physical, elemFields, elemData, nodeFields, nodeData, nodesInElements)
    return

tag = ''
files = sys.argv[1:]
main(files,tag)
