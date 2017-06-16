#!/sw/bin/python

#python findContact.py post/dump104000.generation post/generation.104000.contact_dump 
#python findContact.py post/dump18000.generation post/generation.18000.contact_dump 
import sys
import math

#general:
#ITEM: ATOMS id type type x y z ix iy iz vx vy vz fx fy fz omegax omegay omegaz radius c_contactnum

#histogram:   --> cellDatatoPointData --> Histogram
#line thickness: ---> cellDatatoPointData ---> ExtractSurface --> Tube

def readGeneralFile(file,limits):
    radius = {}
    position = {}
    physical = {}

    xmin = limits[0][0]
    xmax = limits[0][1]
    ymin = limits[1][0]
    ymax = limits[1][1]
    zmin = limits[2][0]
    zmax = limits[2][1]

    #might consider a percentage-based BC

    atomFlag = -1
    with open(file,'r') as data:
        for line in data:
            ls = line.split()
            if (ls[0:2] == ['ITEM:','ATOMS']):
                atomFlag = 0
                continue
            if (atomFlag >= 0):
                atomFlag += 1
            if (atomFlag >= 1):
                ls = line.split()
                idx = int(ls[0])
                xc = float(ls[3])
                yc = float(ls[4])
                zc = float(ls[5])
                rad = float(ls[18])
                pos = [xc, yc, zc]

                radius[idx] = [rad]
                position[idx] = pos

                bc = []
                if xc - rad <= xmin: 
                    bc.append(1) 
                if xc + rad >= xmax:
                    bc.append(2) 
                if yc - rad <= ymin:
                    bc.append(3) 
                if yc + rad >= ymax:
                    bc.append(4) 
                if zc - rad <= zmin:
                    bc.append(5) 
                if zc + rad >= zmax: 
                    bc.append(6)
                physical[idx] = bc 

    return position,radius,physical

def readContactFile(file):
    pos = []
    limits = []
    ids = []
    forces = []
    deltas = []
    areas = []
    xmin = None
    xmax = None
    ymin = None
    ymax = None
    zmin = None
    zmax = None
    nAtom = None
    timestep = None
    boxBounds = -1
    atomFlag = -1
    tsFlag = -1
    with open(file,'r') as data:
        for line in data:
            if (boxBounds == -1):
                ls = line.split()
                if (ls[0:2] == ['ITEM:','TIMESTEP']):
                    tsFlag = 0
                    continue
                elif (ls[0:3] == ['ITEM:','NUMBER','OF']):
                    atomFlag = 0
                    continue
                elif (ls[0:3] == ['ITEM:','BOX','BOUNDS']):
                    boxBounds = 0
                    continue
                if (atomFlag == 0):
                    nAtom = int(line)
                    atomFlag = 1
                    continue
                if (tsFlag == 0):
                    timestep = int(line)
                    tsFlag = 1
                    continue
            elif (boxBounds == 0):
                ls = line.split()
                xmin = float(ls[0])
                xmax = float(ls[1])
                boxBounds += 1
                continue
            elif (boxBounds == 1):
                ls = line.split()
                ymin = float(ls[0])
                ymax = float(ls[1])
                boxBounds += 1
                continue
            elif (boxBounds == 2):
                ls = line.split()
                zmin = float(ls[0])
                zmax = float(ls[1])
                boxBounds += 1
                continue
            elif (boxBounds == 3):
                ls = line.split()
                if (ls[0:2] == ['ITEM:','ENTRIES']):
                    boxBounds += 1
                    continue
            elif (boxBounds == 4):
                ls = line.split()
                idx = int(ls[0])
                x1c = float(ls[1])
                y1c = float(ls[2])
                z1c = float(ls[3])
                x2c = float(ls[4])
                y2c = float(ls[5])
                z2c = float(ls[6])
                id1 = int(ls[7])
                id2 = int(ls[8])
                periodic = int(ls[9]) #1: periodic contact, 0: direct contact
                fx = float(ls[10])
                fy = float(ls[11])
                fz = float(ls[12])
                contactArea = float(ls[13])
                delta = float(ls[14]) #r1 + r2 - dist

                force = math.sqrt(fx*fx + fy*fy + fz*fz)

                #record quantities
                #if (1):
                if (periodic == 0):
                    ids.append([id1,id2])
                    pos.append([ [x1c, y1c, z1c], [x2c, y2c, z2c] ])
                    forces.append(force) #[fx, fy, fz])
                    deltas.append(delta)
                    areas.append(contactArea)
      
    print timestep
    print nAtom
    print xmin,xmax,ymin,ymax,zmin,zmax 
    limits = [ [xmin,xmax], [ymin,ymax], [zmin,zmax] ]
    return ids,pos,forces,areas,deltas,limits

def writePhysicalNodes(file, elemCount, verts, physical):

    keys = verts.keys()
    extraElem = 0

    for i in range(0, len(keys)):
        key = keys[i]
        bc = physical[key]
        for item in bc:
            extraElem += 1
            file.write('{} {} {} {} {}\n'.format(elemCount + extraElem, 15, 1, item, key) )   #"element" number, type 15=point, physical ID, nodeID
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
    file.write('6\n')
    file.write('1 1 xNeg\n')
    file.write('2 2 xPos\n')
    file.write('3 3 yNeg\n')
    file.write('4 4 yPos\n')
    file.write('5 5 zNeg\n')
    file.write('6 6 zPos\n')
    file.write('$EndPhysicalNames\n')
 
def writeNodeData(file, nodeIds, field, fieldname):
    #nodeData as dict
#    print field   
#    print 'nodes: ',nodeIds    #contacts: [id1, id2]
    
    data0 = field[nodeIds[0]]  #just to get data length

    file.write('$NodeData\n')
    file.write('{}\n'.format(1))    #one string tag
    file.write('"{}"\n'.format(fieldname))    #  the name of the view
    file.write('{}\n'.format(1))    #one real tag (time)
    file.write('{}\n'.format(0.0))    #  the time value
    file.write('{}\n'.format(3))    #three integer tags (timestep index, number of components, count
    file.write('{}\n'.format(0))    #  the timestep
    file.write('{}\n'.format(len(data0) ))    #  number of components
    file.write('{}\n'.format(len(nodeIds)))    #  count

    for i in range(0, len(nodeIds)):
        key = nodeIds[i]
        string = ' '.join(str(x) for x in field[key])
        file.write('{} {}\n'.format(key, string) )   #test data
    file.write('$EndNodeData\n')
    return
 
def writeMshFile(filename, verts, edges, physical, elemFields, elemData, nodeFields, nodeData):
#    print 'contacts:',len(edges)
#    print 'nodes:',len(verts.keys())
#    assert( len(ids) == len(verts.keys()) ) #should be the same thing   #definitely not

    physicalGroups = 1 

    keys = verts.keys() 
    print '# nodes = ',len(keys)
    print '# elems = ',len(edges)

    if (filename == None):
        return
    with open(filename, 'w') as file:
        file.write('$MeshFormat\n')
        file.write('2.2 0 8\n')
        file.write('$EndMeshFormat\n')
        file.write('$Nodes\n')
        file.write('{}\n'.format(len(keys)))
        for i in range(0, len(keys)):
            vert = keys[i]
            coord = verts[vert]
            file.write('{} {} {} {}\n'.format(vert, coord[0], coord[1], coord[2])) 
        file.write('$EndNodes\n')

        totalElements = len(edges)
        if (physicalGroups == 1):
            for i in range(0, len(keys)):
                vert = keys[i]
                totalElements += len(physical[vert])
            print '# physBC = ',totalElements-len(edges)
        file.write('$Elements\n')
        file.write('{}\n'.format(totalElements))
        for i in range(0, len(edges)):
            edge = edges[i]
            file.write('{} {} {} {} {}\n'.format(i+1, 1, 0, edge[0], edge[1])) #type: 1 = 2-node line, number of physical group tags = 0 
        if (physicalGroups == 1):
            writePhysicalNodes(file, len(edges), verts, physical)
        file.write('$EndElements\n')

        writePhysicalNames(file)

        assert( len(elemFields) == len(elemData) )
        for i in range(0, len(elemFields)):
            writeElementData(file, elemData[i], elemFields[i]) 

        assert( len(nodeFields) == len(nodeData) )
        for i in range(0, len(nodeFields)):
            writeNodeData(file, keys, nodeData[i], nodeFields[i]) 
        #do I need physical groups??? probably...
    return 

def getGeom(ids, pos):
    assert( len(ids) == len(pos) )
    verts = {}    #dict: id-coords
    edges = []    #id1, id2
    for i in range(0, len(ids)):
        id1 = ids[i][0] 
        id2 = ids[i][1]
        pos1 = pos[i][0] 
        pos2 = pos[i][1]

        #create verts and edges
        verts[id1] = pos1 
        verts[id2] = pos2

        if (id1 < id2):
            edges.append([id1, id2])
        else:
            edges.append([id2, id1])

    return verts,edges

def main(files):
    assert( len(files) % 2 == 0)

    for i in range(0, len(files)/2):
        generalFile = files[i]
        contactFile = files[i + len(files)/2]    
        outFile = '.'.join(contactFile.split('.')[:-1])+'_contact.msh'
 
        print "*Reading contact file ",contactFile
        ids,pos,force,area,delta,limits=readContactFile(contactFile)
        print "*Reading general file ",generalFile
        position,radius,physical=readGeneralFile(generalFile,limits)
        print "*Getting geometry"
        verts,edges = getGeom(ids,pos)
        print "*Writing msh file: ",outFile

        elemFields = ['force','delta','area']
        elemData = [force, area, delta]
        nodeFields = ['position','radius']
        nodeData = [position,radius]
        writeMshFile(outFile, verts, edges, physical, elemFields, elemData, nodeFields, nodeData)
    return

#general-output*.dump, contact_output*.dump (lengths must match)
files = sys.argv[1:]
main(files)

