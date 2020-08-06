from ovito.data import *
import numpy as np
import sys

#---set molcule ID
ID=1

def vec_rot(rot,invec,outvec):
    outvec[0]=invec[0]*rot[0,0] +invec[1]*rot[1,0] +invec[2]*rot[2,0]
    outvec[1]=invec[0]*rot[0,1] +invec[1]*rot[1,1] +invec[2]*rot[2,1]
    outvec[2]=invec[0]*rot[0,2] +invec[1]*rot[1,2] +invec[2]*rot[2,2]

def modify(frame,data, odata):
    #---edit according to your environment
    path='C:/Users/murasima/Desktop/'
    filename='rotation.txt'
    path_to_file=path+filename
    rotdata=[]
    with open(path_to_file,'r') as f:
        il=0
        for line in f:
            line2=line.split()
            il += 1
            if il == 1:
                print(line)
            elif il == 2:
                print(line)
            elif il > 2:
                rotdata.append(line2)
    print ("=== %s ===" % frame)
    iframe=int(frame)
    rot=np.zeros((3,3))
    rot[0,0]=float(rotdata[iframe][1])
    rot[0,1]=float(rotdata[iframe][2])
    rot[0,2]=float(rotdata[iframe][3])
    rot[1,0]=float(rotdata[iframe][4])
    rot[1,1]=float(rotdata[iframe][5])
    rot[1,2]=float(rotdata[iframe][6])
    rot[2,0]=float(rotdata[iframe][7])
    rot[2,1]=float(rotdata[iframe][8])
    rot[2,2]=float(rotdata[iframe][9])

    if data.particles != None:
        cell=data.cell
        ocell=odata.cell_
        ipos=data.particles_['Position_']
        opos=odata.particles_['Position_']

        org=[cell[0][3],cell[1][3],cell[2][3]]
        e0=[cell[0][0],cell[1][0],cell[2][0]]
        e1=[cell[0][1],cell[1][1],cell[2][1]]
        e2=[cell[0][2],cell[1][2],cell[2][2]]
        er0=np.zeros(3)
        er1=np.zeros(3)
        er2=np.zeros(3)
        print(org)
        print(e0)
        print(e1)
        print(e2)

        # e0*rot
        vec_rot(rot,e0,er0)
        # e1*rot
        vec_rot(rot,e1,er1)
        # e2*rot
        vec_rot(rot,e2,er2)

        ocell[0][0]=er0[0]
        ocell[1][0]=er0[1]
        ocell[2][0]=er0[2]

        ocell[0][1]=er1[0]
        ocell[1][1]=er1[1]
        ocell[2][1]=er1[2]

        ocell[0][2]=er2[0]
        ocell[1][2]=er2[1]
        ocell[2][2]=er2[2]
        
        transparency=odata.particles_.create_property('Transparency')
        transparency[:]=0.9
        pid=odata.particles_['Particle Identifier']
        molid=odata.particles_['Molecule Identifier']
        M=max(molid[:])
        NM=molid.size

        flag=np.zeros(NM)
        
        for i in range(transparency.size):
            mid=molid[i]
            if mid == ID:
                transparency[i]=0.0

        r_tmp=np.zeros(3)
        r_rot=np.zeros(3)
        for i in range(NM):
            ii=molid[i]-1
            r_tmp[0]=ipos[i,0]
            r_tmp[1]=ipos[i,1]
            r_tmp[2]=ipos[i,2]
            vec_rot(rot,r_tmp,r_rot)
            opos[i,0]=r_rot[0]
            opos[i,1]=r_rot[1]
            opos[i,2]=r_rot[2]

        print("There are %i particles with the following properties:" % data.particles.count)
        for property_name in data.particles.keys():
            print("  '%s'" % property_name)
            
