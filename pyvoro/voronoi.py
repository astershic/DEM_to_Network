#!/sw/bin/python
import colorsys
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
import pyvoro
import scipy
import StringIO
import sys
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#from scipy.spatial import Voronoi

#https://github.com/JensWehner/KAI-Model/blob/master/lib/__geometry__.py

#http://math.lbl.gov/voro++/doc/voro++_overview.pdf
#http://math.lbl.gov/voro++/examples/radical/
#http://pydoc.net/Python/pyvoro/1.3.2/pyvoro/
#https://github.com/joe-jordan/pyvoro


class struct(object):
    pass


def calcPhysicalNode(verts, limits):

    physical = {}

    #is this node on the limits?
    xmin = limits[0][0]
    xmax = limits[0][1]
    ymin = limits[1][0]
    ymax = limits[1][1]
    zmin = limits[2][0]
    zmax = limits[2][1]

    for id in range(0, len(verts)):
        pt = verts[id]
        bc = []
        if np.allclose( pt[0], xmin):
            bc.append(1) 
        if np.allclose( pt[0], xmax):
            bc.append(2) 
        if np.allclose( pt[1], ymin):
            bc.append(3) 
        if np.allclose( pt[1], ymax):
            bc.append(4) 
        if np.allclose( pt[2], zmin):
            bc.append(5) 
        if np.allclose( pt[2], zmax):
            bc.append(6) 
        physical[id] = bc 

    return physical


def calcPhysicalEdge(edges, verts, limits):

    physical = {}

    #is this node on the limits?
    xmin = limits[0][0]
    xmax = limits[0][1]
    ymin = limits[1][0]
    ymax = limits[1][1]
    zmin = limits[2][0]
    zmax = limits[2][1]

    for id in range(0, len(edges)):
        edge = edges[id]
        points = [verts[edge[0]], verts[edge[1]]]
        bc = []
        if np.allclose( [pt[0] for pt in points], xmin):
            bc.append(1) 
        if np.allclose( [pt[0] for pt in points], xmax):
            bc.append(2) 
        if np.allclose( [pt[1] for pt in points], ymin):
            bc.append(3) 
        if np.allclose( [pt[1] for pt in points], ymax):
            bc.append(4) 
        if np.allclose( [pt[2] for pt in points], zmin):
            bc.append(5) 
        if np.allclose( [pt[2] for pt in points], zmax):
            bc.append(6) 
        physical[id] = bc 

    return physical
def readETHZ(filename, scale=1.0):
    limits = []
    points = []
    radii = []
    ids = []
    timestep = filename.split('/')[-1].split('.')[:-1]
    xMin = None
    xMax = None
    yMin = None
    yMax = None
    zMin = None
    zMax = None

    with open(filename, 'r') as file:
        lineNum = 0
        for line in file:
            lineNum += 1
            if (lineNum > 1):
                lineStrip = line.strip().split(',')
                id = int(lineStrip[0])
                vol = float(lineStrip[1])
                xmin = float(lineStrip[2]) * scale
                ymin = float(lineStrip[3]) * scale
                zmin = float(lineStrip[4]) * scale
                xmax = float(lineStrip[5]) * scale
                ymax = float(lineStrip[6]) * scale
                zmax = float(lineStrip[7]) * scale
                xc = float(lineStrip[8]) * scale
                yc = float(lineStrip[9]) * scale
                zc = float(lineStrip[10]) * scale

                rad = math.pow( (3.0 * vol) / (4.0 * math.pi), 1.0/3.0) * scale

                ids.append(id)
                points.append([xc,yc,zc])
                radii.append(rad)

                #limits
                if (lineNum == 2):
                    xMin = xmin
                    xMax = xmax
                    yMin = ymin
                    yMax = ymax
                    zMin = zmin
                    zMax = zmax
                if (xmin < xMin):
                    xMin = xmin
                if (xmax > xMax):
                    xMax = xmax
                if (ymin < yMin):
                    yMin = ymin
                if (ymax > yMax):
                    yMax = ymax
                if (zmin < zMin):
                    zMin = zmin
                if (zmax > zMax):
                    zMax = zmax
         
    limits = [[xMin, xMax], [yMin, yMax], [zMin, zMax]]

    return timestep, ids, points, radii, limits

def readLiggghts(filename):
    limits = []
    points = []
    radii = []
    ids = []
    timestep = None
    with open(filename, 'r') as file:
        inputSection = 0
        for line in file:
            if (inputSection == 3):
                lineStrip = line.strip().split()
                id = int(lineStrip[0])
                x = float(lineStrip[1])
                y = float(lineStrip[2])
                z = float(lineStrip[3])
                rad = float(lineStrip[4])
                ids.append(id)
                points.append([x,y,z])
                radii.append(rad)
                continue
            elif (inputSection == 2):
                if (line[0:11] == 'ITEM: ATOMS'):
                    inputSection = 3
                    continue
                if (len(limits) < 3):
                    lineStrip = line.strip().split()
                    lower = float(lineStrip[0])
                    upper = float(lineStrip[1])
                    limits.append([lower, upper])
                    continue
            elif (inputSection == 1):
                if (line[0:16] == 'ITEM: BOX BOUNDS'):
                    inputSection = 2
                    continue
                if (timestep == None):
                    timestep = int(line)
                    continue
            elif (inputSection == 0):
                if (line[0:14] == 'ITEM: TIMESTEP'):
                    inputSection = 1
                    continue
            else:
                assert(1 == 0)

    return timestep, ids, points, radii, limits

def writeVolumesFile(filename, ids, volumes):
    if (filename == None):
        return
    with open(filename, 'w') as file:
        file.write('#volumes: ID, volume\n') 
        for id,vol in zip(ids,volumes):
            file.write('{} {}\n'.format(id,vol)) 
    return

def writeFaceNormalsFile(filename, normals):
    if (filename == None):
        return
    with open(filename, 'w') as file:
        file.write('#face unit normals: nx, ny, nz\n') 
        for vec in normals:
            file.write('{} {} {}\n'.format(vec[0], vec[1], vec[2]))
    return


def writePhysicalNodes(file, elemCount, ids, physical):

    extraElem = 0

    for i in range(0, len(ids)):
        bc = physical[i]
        for item in bc:
            extraElem += 1
            file.write('{} {} {} {} {}\n'.format(elemCount + extraElem, 15, 1, item, i+1) )   #"element" number, type 15=point, physical ID, nodeID (1-indexing)
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

def writeMshFile(filename, verts, edges, limits):
    physicalGroups = 1

    physicalNode = calcPhysicalNode(verts, limits)
    physicalEdge = calcPhysicalEdge(edges, verts, limits)

    print '# nodes = ',len(verts)
    print '# elems = ',len(edges)

    if (filename == None):
        return
    with open(filename, 'w') as file:
        file.write('$MeshFormat\n') 
        file.write('2.2 0 8\n') 
        file.write('$EndMeshFormat\n') 
        file.write('$Nodes\n') 
        file.write('{}\n'.format(len(verts))) 
        for i in range(0, len(verts)):
            vert = verts[i]
            file.write('{} {} {} {}\n'.format(i+1, vert[0], vert[1], vert[2]))    #1-based indexing 
        file.write('$EndNodes\n') 
        file.write('$Elements\n')


        totalElements = len(edges)
        if (physicalGroups == 1):
            for i in range(0, len(verts)):
                totalElements += len(physicalNode[i])
            print '# physBC = ',totalElements-len(edges)
        
        file.write('{}\n'.format(totalElements))  #number of elements
        for i in range(0, len(edges)):
            edge = edges[i]
            if (physicalGroups == 1):
                bc = physicalEdge[i]
                if len(bc) == 0:
                    file.write('{} {} {} {} {}\n'.format(i+1, 1, 0, edge[0]+1, edge[1]+1)) #type: 1 = 2-node line, number of physical group tags = 0 
                else:
                    file.write('{} {} {} {} {} {}\n'.format(i+1, 1, len(bc), ' '.join([str(x) for x in bc]), edge[0]+1, edge[1]+1)) #type: 1 = 2-node line, number of physical group tags = 0 
            else:
                file.write('{} {} {} {} {}\n'.format(i+1, 1, 0, edge[0]+1, edge[1]+1)) #type: 1 = 2-node line, number of physical group tags = 0 
        if (physicalGroups == 1):
            writePhysicalNodes(file, len(edges), verts, physicalNode)
        file.write('$EndElements\n')

        writePhysicalNames(file) 
    return

def getNormal(faceVert, vertices):
    v0 = np.array( vertices[ faceVert[0] ] )
    v1 = np.array( vertices[ faceVert[1] ] )
    v2 = np.array( vertices[ faceVert[2] ] )
    normal = np.cross(v1-v0, v2-v0)
    normal = normal / np.linalg.norm(normal)
    return list(normal)

def findVert(vertex, uniqueVert):
    for i in range(0, len(uniqueVert)):
        if (np.allclose(list(vertex), uniqueVert[i])):
            return i
    assert(1 == 0)

def findVertSet(vertex, uniqueVert):
    for vert in uniqueVert: 
        if (np.allclose(list(vertex), vert)): 
            return vert
    return -1

def hashFn(vec):
    a = vec[0]
    b = vec[1]
    c = vec[2]
    dec = 4
    return hash(tuple([round(a,dec), round(b,dec), round(c,dec)])) 

def getVertices(cells):
    vertices = []
    vertSet = set()
    vertHashTable = dict()
    for i in range(0, len(cells)):
#        if (i % (len(cells)/10) == 0):
#            print '  {} of {}'.format(i, len(cells))
        cell = cells[i]
        verts = cell['vertices']
        cell['hashes'] = [] 
        for vert in verts:
            initLength = len(vertSet)

            hashVal = hashFn(vert)
            vertSet.add(hashVal)
 
            cell['hashes'].append(hashVal)
            #if new, add hash to dict
            if (len(vertSet) > initLength):
                assert( len(vertSet) == initLength + 1)
                vertices.append(tuple(vert))  #set version
                vertHashTable[hashVal] = initLength #0-based indexing, internally 
    return vertices,vertHashTable

def plotSphere(ax, col, centroid, rad):
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)*rad + centroid[0]
    y = np.sin(u)*np.sin(v)*rad + centroid[1]
    z = np.cos(v)*rad + centroid[2]
    ax.plot_wireframe(x, y, z, color=col)

def initPlot(N):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect("equal")
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    color = RGB_tuples

    return fig, ax, color

def atLimits(points,limits):
    #is this face on the limits?
    xmin = limits[0][0]
    xmax = limits[0][1]
    ymin = limits[1][0]
    ymax = limits[1][1]
    zmin = limits[2][0]
    zmax = limits[2][1]
    if np.allclose( [pt[0] for pt in points], xmin):
        return 1 
    if np.allclose( [pt[0] for pt in points], xmax):
        return 1 
    if np.allclose( [pt[1] for pt in points], ymin):
        return 1 
    if np.allclose( [pt[1] for pt in points], ymax):
        return 1 
    if np.allclose( [pt[2] for pt in points], zmin):
        return 1 
    if np.allclose( [pt[2] for pt in points], zmax):
        return 1 
    return 0

def plotCell(ax, col, cell, rad, limits, uniqueVert, vertHashTable):
    normals = []
    vertices = cell['vertices']
    verts = np.array(vertices)
    #plot vertices
    if (ax != None):
        ax.scatter(verts[:,0], verts[:,1], verts[:,2])
    #connect vertices with lines
    faces = cell['faces']

    edges = set()
 
#    for face in faces:
    for k in range(0, len(faces)):
        face = faces[k]
        faceVert = face['vertices']
        pverts = [zip(verts[faceVert,0], verts[faceVert,1], verts[faceVert,2])]

        normals.append(getNormal(faceVert, vertices))

        for i in range(0, len(faceVert)):
            j = (i+1) % len(faceVert) 
            m = faceVert[i]
            n = faceVert[j]
        #    if (ax != None):
#                ax.plot( verts[[m,n],0], verts[[m,n],1], verts[[m,n],2], color=col )

            #record edges
            if (1):   #all faces
#            if (not atLimits(pverts[0], limits)):    #EXCLUDE faces on limits
                everts = [zip(verts[[m,n],0], verts[[m,n],1], verts[[m,n],2])]
                if (1):  #all edges
#                if (not atLimits(everts[0], limits)):   #also exclude edges on limits
                    if (ax != None):
                        ax.plot( verts[[m,n],0], verts[[m,n],1], verts[[m,n],2], color=col )

                    #store edge (find global vertices)
                    mHash = cell['hashes'][m]
                    nHash = cell['hashes'][n]
                    mGlobal = vertHashTable[mHash] 
                    nGlobal = vertHashTable[nHash]
                    entry = [mGlobal, nGlobal]
                    entry.sort()
                    edges.add( tuple(entry) ) 
        #draw face
        #EXCLUDE on limits
        if (ax != None):
            if (not atLimits(pverts[0], limits)):
                ax.add_collection3d(Poly3DCollection(pverts, facecolors=col))

    #transcluent 3d object

    # draw sphere
    centroid = cell['original']
    if (ax != None):
        plotSphere(ax, col, centroid, rad)

    #output edges
    volume = cell['volume']
    return edges, volume, normals


def main(filename,plot='no'):
    varlist = struct()

    if (filename.split('.')[-1] == 'voro_input'):
        timestep,ids,points,radii,limits = readLiggghts(filename)
        periodic = [True, True, False]
    elif (filename.split('.')[-1] == 'txt'):
        timestep,ids,points,radii,limits = readETHZ(filename)
        periodic = [False,False,False]
    else:
        assert(1 == 0)

#    ids = ids[0:5500]
#    points = points[0:5500]
#    radii = radii[0:5500]

    #input
#    points = [[1.0, 2.0, 3.0], [4.0, 5.5, 6.0], [0, 5.0, 0] ] # point positions
#    limits = [[0.0, 10.0], [0.0, 10.0], [0.0, 10.0]] # limits
#    radii=[1.3,1.4, 2] # particle radii -- optional, and keyword-compatible arg.
    block_size = 4*max(radii) #2.0 # block size : maximum distance that two adjacent particles might be
    ids = [i for i in range(1,1+len(points))]
    base = '.'.join(os.path.basename(filename).split('.')[:-1])
    meshFilename = base+'_edges.msh'  #None 
    volFilename = base+'_vols.txt'    #None
    normsFilename = base+'_norms.txt' #None

    print "*Calculating voronoi tessellation for ",filename
    while True:
        try:
            varlist.cells = pyvoro.compute_voronoi( points, limits, block_size, radii=radii, periodic=periodic)
            break
        except:
            print "*Shrinking & Retrying"
            radii = list(np.array(radii) * 0.99)
 
    count = len(varlist.cells)

    #Count total number of vertices for check
    totalVert = 0
    for i in range(0,len(varlist.cells)):
        cell = varlist.cells[i]
        totalVert += len(cell['vertices'])

    print "*Post-processing voronoi tessellation for ",filename
    if (plot == 'yes'):
        fig,ax,color = initPlot( count )
    else:
        ax = None
        color = [None]*count

    #unique list of vertices
    print " *Getting unique vertices"
    uniqueVert,vertHashTable = getVertices(varlist.cells)

    assert(totalVert != len(uniqueVert) )

    #print varlist.cells
    print " *Getting edges"
    edges = set()
    normals = []
    volumes = []
    for i in range(0, count):
        if (i % (count/10) == 0):
            print ' c {} of {}'.format(i, count)
        cell = varlist.cells[i]
        rad = radii[i]
        col = color[i]
        cellEdge,cellVol,cellNorms = plotCell(ax,col,cell,rad,limits,uniqueVert,vertHashTable)
        edges = edges.union(cellEdge)
        normals.extend(cellNorms)
        volumes.append(cellVol)
#    print edges
    if (plot == 'yes'):
        plt.show()


    print "*Writing to files"
    #print global skeleton
    writeMshFile(meshFilename, uniqueVert, list(edges), limits)
 
    #write volumes to file
    writeVolumesFile(volFilename, ids, volumes)
  
    #write faceNormals to file
    writeFaceNormalsFile(normsFilename, normals)

#need to read in liggghts file, or be called by liggghts python
filename = sys.argv[1]
adjustList = sys.argv[2:]
#plot='yes'
plot='no'
main(filename,plot=plot)
