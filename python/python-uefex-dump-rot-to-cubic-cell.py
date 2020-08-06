import numpy as np
import sys
import os
import math

#set density
density=0.85

def vec_rot(rot,invec,outvec):
    outvec[0]=invec[0]*rot[0,0] +invec[1]*rot[1,0] +invec[2]*rot[2,0]
    outvec[1]=invec[0]*rot[0,1] +invec[1]*rot[1,1] +invec[2]*rot[2,1]
    outvec[2]=invec[0]*rot[0,2] +invec[1]*rot[1,2] +invec[2]*rot[2,2]

    
file_dump=sys.argv[1]
filename_dump=str(file_dump)

file_rot=sys.argv[2]
filename_rot=str(file_rot)

rotdata=[]
iframe=0

w=open('rotation.'+filename_dump,'w')

with open(filename_rot,'r') as f:
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

with open(filename_dump,'r') as f:
    ick=0
    boxcount=0
    for line in f:
        line2=line.split()
        if ick == 0:
            if line2[1] == 'TIMESTEP':
                ick=1
                continue
            elif line2[1] == 'NUMBER':
                ick=2
                continue
            elif line2[1] == 'BOX':
                ick=3
                continue
            elif line2[1] == 'ATOMS':
                ick=4
                continue

        if ick == 1:
            print("Step="+line2[0])
            timestep=int(line2[0])
            ick=0

        if ick == 2:
            NM=int(line2[0])
            pos=np.zeros((NM,3))
            opos=np.zeros((NM,3))
            vel=np.zeros((NM,3))
            ovel=np.zeros((NM,3))
            force=np.zeros((NM,3))
            oforce=np.zeros((NM,3))
            pid=np.zeros(NM)
            mid=np.zeros(NM)
            tid=np.zeros(NM)
            
            # ghost particle
            gNM=int(0.5*NM)
            gopos=np.zeros((gNM,3))
            govel=np.zeros((gNM,3))
            goforce=np.zeros((gNM,3))
            gpid=np.zeros(gNM)
            gmid=np.zeros(gNM)
            gtid=np.zeros(gNM)
            gcount=0

            ii=0
            ick=0

        if ick == 3:
            if boxcount==0:
                xlo=float(line2[0])
                xhi=float(line2[1])
                xy=float(line2[2])
                boxcount += 1
            elif boxcount==1:
                ylo=float(line2[0])
                yhi=float(line2[1])
                xz=float(line2[2])
                boxcount += 1
            elif boxcount==2:
                zlo=float(line2[0])
                zhi=float(line2[1])
                yz=float(line2[2])
                boxcount = 0
                ick=0

        if ick == 4:
            pid[ii]=int(line2[0])
            mid[ii]=int(line2[1])
            tid[ii]=int(line2[2])
            pos[ii,0]=float(line2[3])
            pos[ii,1]=float(line2[4])
            pos[ii,2]=float(line2[5])
            vel[ii,0]=float(line2[6])
            vel[ii,1]=float(line2[7])
            vel[ii,2]=float(line2[8])
            force[ii,0]=float(line2[9])
            force[ii,1]=float(line2[10])
            force[ii,2]=float(line2[11])
            ii += 1
            if ii == NM:
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
                e0=[xhi-xlo-abs(xy)-abs(xz),0.0,0.0]
                e1=[xy,yhi-ylo-abs(yz),0.0]
                e2=[xz,yz,zhi-zlo]
                org=[0.0,0.0,0.0]
                er0=np.zeros(3)
                er1=np.zeros(3)
                er2=np.zeros(3)
                # e0*rot
                vec_rot(rot,e0,er0)
                # e1*rot
                vec_rot(rot,e1,er1)
                # e2*rot
                vec_rot(rot,e2,er2)
                V=NM/density
                L=V**(1./3.)
                flag=np.zeros(NM)
                r_tmp=np.zeros(3)
                r_rot=np.zeros(3)
                v_tmp=np.zeros(3)
                v_rot=np.zeros(3)
                f_tmp=np.zeros(3)
                f_rot=np.zeros(3)
                for i in range(NM):
                    j=mid[i]-1
                    for ilx in range(-1,2):
                        for ily in range(-1,2):
                            for ilz in range(-1,2):
                                r_tmp[0]=pos[i,0] + float(ilx)*e0[0] + float(ily)*e1[0] + float(ilz)*e2[0]
                                r_tmp[1]=pos[i,1] + float(ilx)*e0[1] + float(ily)*e1[1] + float(ilz)*e2[1]
                                r_tmp[2]=pos[i,2] + float(ilx)*e0[2] + float(ily)*e1[2] + float(ilz)*e2[2]
                                v_tmp[0]=vel[i,0]
                                v_tmp[1]=vel[i,1]
                                v_tmp[2]=vel[i,2]
                                f_tmp[0]=force[i,0]
                                f_tmp[1]=force[i,1]
                                f_tmp[2]=force[i,2]
                                vec_rot(rot,r_tmp,r_rot)
                                vec_rot(rot,v_tmp,v_rot)
                                vec_rot(rot,f_tmp,f_rot)

                                if abs(r_rot[0]) <= 0.5*L:
                                    if abs(r_rot[1]) <= 0.5*L:
                                        if abs(r_rot[2]) <= 0.5*L:
                                            if flag[i]<1:
                                                opos[i,0]=r_rot[0]
                                                opos[i,1]=r_rot[1]
                                                opos[i,2]=r_rot[2]
                                                ovel[i,0]=v_rot[0]
                                                ovel[i,1]=v_rot[1]
                                                ovel[i,2]=v_rot[2]
                                                oforce[i,0]=f_rot[0]
                                                oforce[i,1]=f_rot[1]
                                                oforce[i,2]=f_rot[2]
                                                flag[i]=1
                                            else:
                                                gpid[gcount]=pid[i]
                                                gmid[gcount]=mid[i]
                                                gtid[gcount]=tid[i]
                                                gopos[gcount,0]=r_rot[0]
                                                gopos[gcount,1]=r_rot[1]
                                                gopos[gcount,2]=r_rot[2]
                                                govel[gcount,0]=v_rot[0]
                                                govel[gcount,1]=v_rot[1]
                                                govel[gcount,2]=v_rot[2]
                                                goforce[gcount,0]=f_rot[0]
                                                goforce[gcount,1]=f_rot[1]
                                                goforce[gcount,2]=f_rot[2]
                                                gcount += 1
                                                
                ick=0
                iframe+=1
                w.write('ITEM: TIMESTEP\n')
                w.write(str(timestep)+"\n")
                w.write('ITEM: NUMBER OF ATOMS\n')
                iNM=np.sum(flag)
                NMnew=iNM+gcount
                print("Number of Particles="+str(int(NMnew)))
                w.write(str(NMnew)+"\n")
                w.write('ITEM: BOX BOUNDS xy xz yz pp pp pp\n')
                w.write(str(-0.5*L)+" "+str(0.5*L)+" "+str(0)+"\n")
                w.write(str(-0.5*L)+" "+str(0.5*L)+" "+str(0)+"\n")
                w.write(str(-0.5*L)+" "+str(0.5*L)+" "+str(0)+"\n")
                w.write('ITEM: ATOMS id mol type x y z vx vy vz fx fy fz\n')
                for i in range(NM):
                    if flag[i] > 0:
                        w.write(str(pid[i])+" "+str(mid[i])+" "+str(tid[i])+" "
                                +str(opos[i,0])+" "+str(opos[i,1])+" "+str(opos[i,2])+" "
                                +str(ovel[i,0])+" "+str(ovel[i,1])+" "+str(ovel[i,2])+" "
                                +str(oforce[i,0])+" "+str(oforce[i,1])+" "+str(oforce[i,2])+" "
                                +"\n")
                for i in range(gcount):
                    w.write(str(gpid[i])+" "+str(gmid[i])+" "+str(gtid[i])+" "
                            +str(gopos[i,0])+" "+str(gopos[i,1])+" "+str(gopos[i,2])+" "
                            +str(govel[i,0])+" "+str(govel[i,1])+" "+str(govel[i,2])+" "
                            +str(goforce[i,0])+" "+str(goforce[i,1])+" "+str(goforce[i,2])+" "
                            +"\n")

w.close()


