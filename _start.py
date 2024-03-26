from particle import *
from write import *
from tab_pot import *

e = np.linspace(6,10,7) # Ligand Receptor Interaction Energy used
f = np.linspace(.3,.8,7) # Fraction of the nanoparticle surface covered in ligands used 
sphereRad = 7  # Sphere Radius 
offset = 1  # Increase the distance between the nanoparticle and the membrane

core = np.asarray(sphere(sphereRad)) # Nx3 numpy array with nanoparticle coordinates 
memb = np.loadtxt('membrane.txt',skiprows=15,usecols=(0,1),max_rows=3) # Load membrane coordinates
Wc = 1.5 # Weeks-Chandler-Anderson factor
Epsilon = e[-1] # Ligand Receptor Interaction Energy 
binding_f = f[-1] # Fraction of the nanoparticle surface covered in ligands

# Position the nanoparticle in the center of the box, above the membrane
core[:,2] += sphereRad + memb[2,1] + offset
core[:,0] += (memb[0,1]-memb[0,0])/2
core[:,1] += (memb[1,1]-memb[1,0])/2

# Write the necessary LAMMPS files
writeParticle(core,[],binding_f,'Even',memb)
tabs(wc[REPLACE],Epsilon)
