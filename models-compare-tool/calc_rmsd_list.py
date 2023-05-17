"""
This is a program written by Qianyi Cheng
at University of Memphis.
"""

import os, sys, re
import argparse
from numpy import *
from read_write_pdb import *
from scipy.linalg import norm
import pandas as pd


def get_res_atom(ref_pdb):
    pdb, res_info, tot_charge = read_pdb(ref_pdb.strip())
    res_atom = []
    xyz_atom = []
    for l in range(len(pdb)):
        atom = pdb[l]
        key = [atom[5].strip(),atom[6],atom[2].strip()]  ## key=[chainID,resID,atomNAME]
        res_atom.append(key)
        xyz = [atom[8],atom[9],atom[10]]
        xyz_atom.append(xyz)
    return res_atom, xyz_atom


def get_atom_idx(atom1,atom2): 
    idx1 = {}
    idx2 = {}
    for i in range(len(atom1)):
        key = (atom1[i][0],atom1[i][1])
        if key not in idx1.keys():
            idx1[key] = []
            idx2[key] = []
            if atom1[i] in atom2:
                idx1[key].append(i)
                idx2[key].append(atom2.index(atom1[i]))
        else:
            if atom1[i] in atom2:
                idx1[key].append(i)
                idx2[key].append(atom2.index(atom1[i]))

    return idx1,idx2


def calc_dist(a,b):  ## a and b a xyz coordinates array
    return sqrt( (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2 )

def calc_angle2(a,b,c):
    a = array(a)
    b = array(b)
    c = array(c)
    ba = a - b
    bc = c - b
    #cosine_angle = dot(ba,bc)/(norm(ba)*norm(bc))
    #angle = arccos(cosine_angle)
    cosine_angle = dot(ba,bc)
    sinine_angle = norm(cross(ba,bc))
    angle = arctan2(sinine_angle,cosine_angle)
    return angle*180/pi


def calc_dihedral2(a,b,c,d):
    a = array(a)
    b = array(b)
    c = array(c)
    d = array(d)
    ab = -1.0*(b - a)
    bc = c - b
    cd = d - c

    bc /= norm(bc)
    v = ab - dot(ab,bc)*bc
    w = cd - dot(cd,bc)*bc
    x = dot(v,w)
    y = dot(cross(bc,v),w)
    return degrees(arctan2(y,x))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compare optimized structures, save rmsd values to csv file')
    parser.add_argument('-pdbr', dest='pdbr', default=None, help='opt_pdbr')
    parser.add_argument('-pdbl', dest='pdbl', default=None, help='opt_pdbl')

    args = parser.parse_args()
    pdb2 = args.pdbr
    res_atom2, xyz_atom2 = get_res_atom(pdb2) 
    keys = []
    index_df = []
    for atom in read_pdb(pdb2)[0]:
        names = (atom[5].strip(),atom[6])
        if names not in keys:
            keys.append(names)
            index_df.append(atom[5].strip()+'-'+str(atom[6]))

    pdb_list = args.pdbl
    with open(pdb_list,'r') as f:
        lines = f.readlines()
    names = []
    for line in lines:
        names.append(line.split('/')[-1][:-1])

    num = len(lines)

    #rmsd_df = pd.DataFrame(0,columns=names,index=arange(len(keys)))#,index=index_df)
    rmsd_df = pd.DataFrame(0,columns=names,index=index_df)
    for i  in range(len(lines)):
        pdb1 = lines[i]
        res_atom1, xyz_atom1 = get_res_atom(pdb1) 

        idx1, idx2 = get_atom_idx(res_atom1,res_atom2)
#        print(idx1,idx2)
        rmsd = {}
        for key in sorted(idx1):
#            print(key)
            dist = []
            for el in range(len(idx2[key])):
#               print(res_atom1[idx1[key][el]],res_atom2[idx2[key][el]])
#                if res_atom1[idx1[key][el]][2][0] != 'H': ### Do not count H atoms
                dist.append(calc_dist(xyz_atom1[idx1[key][el]],xyz_atom2[idx2[key][el]])) 
            rmsd[key] = mean(array(dist))
        for key in keys:
            if key in rmsd.keys():
#                print(key, rmsd[key])
                rmsd_df.iloc[keys.index(key),i] = rmsd[key]
    rmsd_df.to_csv('rmsd_minus1.csv')


