step=100
frame0=str(step).zfill(3)
name='ring.N400'
print(frame0)
from ovito.io import import_file
filename=name+'.'+str(frame0)+".data"
pipeline=import_file(filename,atom_style="bond")

data=pipeline.compute()
mol=data.particles['Molecule Identifier']
import numpy as np
num_particles=mol.size
mol_max=np.max(mol)

rotdata=[]

with open('rotation.txt','r') as f:
    il=0
    for line in f:
        line2=line.split()
        il+=1
        if il >= 3:
            rotdata.append(line2)

def vec_rot(rot,invec,outvec):
    outvec[0]=invec[0]*rot[0,0] +invec[1]*rot[1,0] +invec[2]*rot[2,0]
    outvec[1]=invec[0]*rot[0,1] +invec[1]*rot[1,1] +invec[2]*rot[2,1]
    outvec[2]=invec[0]*rot[0,2] +invec[1]*rot[1,2] +invec[2]*rot[2,2]

def color_rgb(i,max):
    f=float(i)/float(max)
    r=0
    g=0
    b=0
    if f < 0.25:
        r=1.0
        g=4*f
        b=0.0
    elif f < 0.5:   
        r=1.0-4*(f-0.25)
        g=1.0
        b=0.0
    elif f < 0.75:
        r=0.0
        g=1.0
        b=4*(f-0.5)
    else:
        r=0.0
        g=1.0-4*(f-0.75)
        b=1.0
    return [r,g,b]

def modify(frame,data):
     color=data.particles_.create_property('Color')
     #trans=data.particles_.create_property('Transparency')
     for i in range(num_particles):
         mol_i=mol[i]
         #if mol_i == 1:
         #    trans[i]=0.0
         #else:
         #    trans[i]=0.9
         color[i]=color_rgb(mol_i,mol_max)
pipeline.modifiers.append(modify)


from ovito.modifiers import UnwrapTrajectoriesModifier
pipeline.modifiers.append(UnwrapTrajectoriesModifier())


def modify_pbc(frame,data):
    pid=data.particles['Particle Identifier']
    pos=data.particles_['Position_']
    opos=data.particles_['Position_']
    cell=data.cell
    com=np.zeros((num_particles,3))
    com2=np.zeros((num_particles,3))
    com_pbc=np.zeros((mol_max,3),dtype=int)
    com_pbc_s=np.zeros((mol_max,3))
    num_per_mol=np.zeros(mol_max,dtype=int)
    H=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            H[i,j]=cell[i,j]
    Hinv=np.linalg.inv(H)
    for i in range(num_particles):
        pid_i=pid[i]-1
        mol_i=mol[i]-1
        com[mol_i]+=pos[i]
        num_per_mol[mol_i]+=1
    for i in range(mol_max):
        com[i]/=num_per_mol[i]
    for i in range(mol_max):
        com_pbc[i]=np.ceil(np.dot(Hinv,com[i])-0.5)
    for i in range(num_particles):
        mol_i=mol[i]-1
        opos[i]=pos[i]-np.dot(H,com_pbc[mol_i])

pipeline.modifiers.append(modify_pbc)    


def modify_rotation(frame,data):
    cell=data.cell
    ocell=data.cell_
    ipos=data.particles_['Position_']
    opos=data.particles_['Position_']
    molid=data.particles['Molecule Identifier']
    NM=molid.size
    iframe=int(step)
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

    org=[cell[0][3],cell[1][3],cell[2][3]]
    e0=[cell[0][0],cell[1][0],cell[2][0]]
    e1=[cell[0][1],cell[1][1],cell[2][1]]
    e2=[cell[0][2],cell[1][2],cell[2][2]]
    er0=np.zeros(3)
    er1=np.zeros(3)
    er2=np.zeros(3)
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

    r_tmp=np.zeros(3)
    r_rot=np.zeros(3)
    for i in range(NM):
        r_tmp[0]=ipos[i,0]-org[0]
        r_tmp[1]=ipos[i,1]-org[1]
        r_tmp[2]=ipos[i,2]-org[2]
        vec_rot(rot,r_tmp,r_rot)
        opos[i,0]=r_rot[0]+org[0]
        opos[i,1]=r_rot[1]+org[1]
        opos[i,2]=r_rot[2]+org[2]
    
pipeline.modifiers.append(modify_rotation)

pipeline.add_to_scene()

from ovito.vis import ParticlesVis
particle_vis=pipeline.source.data.particles.vis
particle_vis.radius=0.4

from ovito.vis import BondsVis
bond_vis=pipeline.source.data.particles.bonds.vis
bond_vis.width=0.5
bond_vis.enabled=True

cell_vis=pipeline.source.data.cell.vis
cell_vis.rendering_color=(0.0,0.0,0.0)

from ovito.vis import CoordinateTripodOverlay
tripod = CoordinateTripodOverlay()
tripod.size = 0.07
tripod.offset_x = 0.02
tripod.offset_y = 0.02
lx=data.cell[0,0]
ly=data.cell[1,1]
lz=data.cell[2,2]

import math
from ovito.vis import Viewport
vp = Viewport(type=Viewport.Type.Perspective)
vp.camera_pos=(0,-420,20)
vp.camera_dir=(0,1,0)
vp.overlays.append(tripod)
from ovito.vis import TachyonRenderer
vp.render_image(size=(1600,1600),filename="figure."+str(frame0)+".png",background=(1.0,1.0,1.0),renderer=TachyonRenderer())

