from ovito.data import *
import numpy as np
import sys
ID=1

def modify(frame,data, odata):
    if data.particles != None:
        molid=odata.particles_['Molecule Identifier']
        transparency=odata.particles_.create_property('Transparency')
        transparency[:]=0.9
        for i in range(transparency.size):
            mid=molid[i]
            if mid == ID:
                transparency[i]=0.0

        for property_name in data.particles.keys():
            print("  '%s'" % property_name)
            
