from numpy import *
import os, sys

### In Old X section, bond distance in Bohr   1A = 1.88973 Bohr      ###
###                     angle in Radian 1 Degree = 0.0174533 Radian  ###
########################################################################

with open('%s/1.out'%sys.argv[1]) as f:
    lines = f.readlines()

fl = len(lines)

### keywords section ###
for i in range(fl):
    l = lines[i]
    if " #P " in l:
        if "------" not in lines[i+1]:
            v = (lines[i]+lines[i+1]).split()
        else:
            v = l.split()
        print(v)

### OPT step section for variables ############################################
step_num = []
eigen_l = []
eigen_vec = []
oldx_start = [] # For beginning of oldx newx
oldx_end = []   # For end of oldx newx
convg = []      # For max force/rms force/max displacement/rms displacement
                # with predicted energy change
scf = []
spin = []
dipole = []
for i in range(fl):
    if "Step number" in lines[i]:
        step_num.append(int(lines[i].split()[2]))
    elif "Eigenvalues" in lines[i]:
        eigen_l.append(i)
    elif "Eigenvectors required" in lines[i]:
        eigen_vec.append(i+1)
    elif "Variable" in lines[i]:
        oldx_start.append(i+2)
    elif "Item" in lines[i]:
        oldx_end.append(i-1)
        convg.append(i+1)
    elif "SCF Done:" in lines[i]:
        scf.append(i)
    elif "S**2 before annihilation" in lines[i]:
        spin.append(i)
    elif "Electric dipole moment (input orientation):" in lines[i]:
        dipole.append(i)
    else:
        continue
#print(oldx_start,len(oldx_start))
#print(oldx_end,len(oldx_end))
#print(len(step_num))
#print(len(eigen_l))  ### will be a lot
#print(len(eigen_vec))
#print(len(convg))
#print(len(scf))
#print(len(spin))

var_change = {}
for i in range(len(oldx_start)):
    if i <= len(step_num):
        snum = step_num[i]
    else:
        snum = step_num[-1]+1
    var_change[snum] = []
    for j in range(oldx_start[i],oldx_end[i]+1):
        v = lines[j].split()
        var_change[snum].append([float(v[1]),float(v[-2]),float(v[-1])])
print(var_change[step_num[-1]])        
