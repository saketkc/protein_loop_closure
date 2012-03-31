import os.path as op
import numpy as np
from Bio.PDB import PDBIO
import matplotlib.pyplot as plt
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser()
filepath = '1G59.pdb'
structure = parser.get_structure(op.splitext(filepath)[0], filepath)
model = structure[0]
residues = model.get_residues()
residuecount=0
for residue in residues:
    if residue.id[0]!='W':
	residuecount+=1

print "Total Number of Residues: " + str(residuecount)
residue_for_perturbation = input("Enter Residue no. For Perturbation : ")

chain = model['A']
residue = chain[(' ', int(residue_for_perturbation), ' ')]

atom_N = residue['N']
atom_CA = residue['CA']
atom_C = residue['C']
#print chain.id
#atom_CA.coord[0] = atom_CA.coord[0] +10
#print residue.id[1], residue.resname
#print atom_N.name,atom_N.coord#, atom1.occupancy, atom1.bfactor
#print atom_CA.name,atom_CA.coord[0]#, atom2.occupancy, atom2.bfactor
#first_aa = chain_a[(' ', 1, ' ')]
x_C = atom_C.coord[0]
y_C = atom_C.coord[1]
z_C = atom_C.coord[2]

x_CA = atom_CA.coord[0]
y_CA = atom_CA.coord[1]
z_CA = atom_CA.coord[2]

x_N = atom_N.coord[0]
y_N = atom_N.coord[1]
z_N = atom_N.coord[2]

mean_C = [x_C,y_C,z_C]

covariance = [[1,0,0],[0,1,0],[0,0,1]]
C_x,C_y,C_z = np.random.multivariate_normal(mean_C,covariance,50000).T

for i in range(0,len(C_x)):
    atom_C.coord[0] = C_x[i]
    atom_C.coord[1] = C_y[i]
    atom_C.coord[2] = C_z[i]
    #w = PDBIO()
    #w.set_structure(structure)
    #w.save('1G59-r'+str(i)+'.pdb')
plt.plot(C_x,C_y,'x')

plt.axis('equal')
plt.show()
