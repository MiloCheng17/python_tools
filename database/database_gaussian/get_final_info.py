from numpy import *
import os, sys

### In Old X section, bond distance in Bohr   1A = 1.88973 Bohr      ###
###                     angle in Radian 1 Degree = 0.0174533 Radian  ###
########################################################################

def process_one_outf(outf):
    with open('%s/1.out'%outf) as f:
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
            if 'Freq' in l or 'freq' in l and 'opt' not in l:
                freql = i
            #print(v)
    
    eigen_l = []
    eigen_vec = []
    oldx_start = [] # For beginning of oldx newx
    oldx_end = []   # For end of oldx newx
    convg = []      # For max force/rms force/max displacement/rms displacement
                    # with predicted energy change
    chag = []   # charge
    mult = []   # mulitplicity
    scfe = []   # SCF Done energy
    spin = []   # S**2 values 2
    dipl = []   # dipole moment
    frequency_l = []
    
    for i in range(freql,fl):
        ### Charge and Multiplicity as a single value ###
        if "Charge" in lines[i] and "Multiplicity" in lines[i]:
            v = lines[i].split()
            chag.append(int(v[2]))
            mult.append(int(v[-1]))
        ### Final optimized xyz coordinates in angstrom ###
        elif "Redundant internal coordinates found in file" in lines[i]:
            xyz_s = i+1
        elif "Recover connectivity data from disk" in lines[i]:
            xyz_e = i
        elif "SCF Done:" in lines[i]:
            v = lines[i].split()
            scfe.append(float(v[4]))
        elif "S**2 before annihilation" in lines[i]:
            v = lines[i].split()
            spin.append(float(v[3][:-1]))
            spin.append(float(v[5]))
        ### Then Atomic-Atomic Spin Densities ###
        ### Mulliken charges and spin densities with hydrogens summed into heavy atoms ###
        elif "Dipole moment (field-independent basis, Debye)" in lines[i]:
            dipl.append(float(lines[i+1].split()[-1]))
        elif "Frequencies --" in lines[i]:
            frequency_l.append(i)
        elif "Sum of electronic and zero-point Energies=" in lines[i]:
            zpe = [float(lines[i].split()[-1])]
        elif "Sum of electronic and thermal Energies=" in lines[i]:
            thermal_e = [float(lines[i].split()[-1])]
        elif "Sum of electronic and thermal Enthalpies=" in lines[i]:
            thermal_h = [float(lines[i].split()[-1])]
        elif "Sum of electronic and thermal Free Energies=" in lines[i]:
            thermal_g = [float(lines[i].split()[-1])]
        elif "Eigenvalues ---" in lines[i]:
            eigen_l.append(i)
        elif "Eigenvectors required" in lines[i]:
            eigen_vec.append(i+1)
        elif "Variable" in lines[i]:
            oldx_start.append(i+2)
        elif "Item" in lines[i]:
            oldx_end.append(i)
            convg.append(i+1)
        else:
            continue
    num_atoms = [xyz_e-1-xyz_s]
    #print(chag,mult,num_atoms,xyz_s,xyz_e,scfe,spin,dipl)
    #print(frequency_l[0],eigen_l[0],eigen_vec,oldx_start,oldx_end,convg)
    
    
    ### xyz coordinate in original coordinates from xyz_s, xyz_e ###
    ### has not used yet ###
    final_xyz = []
    for i in range(xyz_s, xyz_e):
        final_xyz.append(lines[i][:-1])
    
    ### Frequency information: Frequency, Red. massed, Frc consts, IR Inten for the first one ###
    freq_info  = []
    for i in range(4):
        v = lines[frequency_l[0]+i].split()
        freq_info.append(float(v[-3]))
    
    eigen_val = [float(lines[eigen_l[0]].split()[-5])]
    
    ### Old X New X in Bohr and Radians ###
    old_x = []
    new_x = []
    #print(oldx_start,oldx_end)
    for j in range(len(oldx_start)):
        for i in range(oldx_start[j], oldx_end[j]):
            v = lines[i].split()
            old_x.append(float(v[1]))
            new_x.append(float(v[-1]))
    
    convg_info = [] ### Maximun Force, RMS Force, Maximum Displacement, RMS Displacement ###
    for i in range(4):
        v = lines[convg[0]+i].split()
        convg_info.append(float(v[-3]))
        convg_info.append(v[-1])
    
    #print(chag,mult,num_atoms,scfe,spin,dipl,zpe,thermal_e,thermal_h,thermal_g)
    #print(final_xyz[0],final_xyz[-1],freq_info, eigen_val,old_x[0],old_x[-1],new_x[0],new_x[-1],convg_info)
    print_list = [outf]
    for print_item in (chag,mult,num_atoms,scfe,spin,dipl,zpe,thermal_e,thermal_h,thermal_g,freq_info,eigen_val,convg_info):
        for i in print_item:
            print_list.append(i)
    return print_list

list_file = genfromtxt('list_dir',dtype=str)
all_list = []
for filename in list_file:
    print_list = process_one_outf(filename)
    all_list.append(print_list)
f = open("all-final.csv",'w')
f.write('Path,Charge,Multiplicity,num_atoms,SCF energy,Spin1,Spin2,Dipole,ZPE,Thermal E,Thermal H,Thermal G,Frequency,RedMass,ForceConst,IR_Intern,Eigen_val,MaxForce,Met,RMSForce,Met,MaxDisp,Met,RMSDisp,Met\n')
for j in range(len(all_list)):
    for i in all_list[j]:
        f.write('%s,'%i)
    f.write('\n')
f.close()
