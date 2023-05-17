"""
This is a program written by Qianyi Cheng
at University of Memphis.
"""

import os, sys, re
import argparse
import numpy as np
from read_write_pdb import *
from scipy.linalg import norm
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as mpl


def get_res_atom(ref_pdb):
    pdb, res_info, tot_charge = read_pdb(ref_pdb.strip())
    res_atom = []
    xyz_atom = []
    for l in range(len(pdb)):
        atom = pdb[l]
        key = [atom[5].strip(),str(atom[6]),atom[2].strip()]  ## key=[chainID,resID,atomNAME]
        res_atom.append(key)
        xyz = [atom[8],atom[9],atom[10]]
        xyz_atom.append(xyz)
    return res_atom, xyz_atom


def calc_dist(a,b):  ## a and b a xyz coordinates array
    return np.sqrt( (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2 )


def write_outf_1header(res_atom,dist_mat,out):
    ### write output file as csv file ###
    ### One header with A-1-CA joined by "-" ###
    header = []
    num_atom = len(res_atom)
    for i in range(num_atom):
        header.append('-'.join(res_atom[i]))
    outf = open(out,'w')
#    outf.write('name')
    for i in range(num_atom):
        outf.write(',%s'%header[i])
    outf.write('\n')
    for i in range(num_atom):
        outf.write('%s'%header[i])
        for j in range(num_atom):
            if j == i:
                outf.write(', 0.0')
            elif j < i:
                outf.write(', %f'%dist_mat[j,i])
            else:
                outf.write(', %f'%dist_mat[i,j])
        outf.write('\n')
    outf.close()                


def write_outf_2header(res_atom,dist_mat,out):
    ### write output file as csv file ###
    ### Two index column chain-resid, atomname ###
    header = []
    header2 = []
    num_atom = len(res_atom)
    for i in range(num_atom):
        header.append('-'.join(res_atom[i]))
        header2.append('-'.join(res_atom[i][:2]))
    outf = open(out,'w')
    outf.write('chain-resid,atom')
    for i in range(num_atom):
        outf.write(',%s'%header[i])
    outf.write('\n')
    for i in range(num_atom):
        outf.write('%s'%header2[i])
        outf.write(',%s'%res_atom[i][2])
        for j in range(num_atom):
            if j == i:
                outf.write(', 0.0')
            elif j < i:
                outf.write(', %f'%dist_mat[j,i])
            else:
                outf.write(', %f'%dist_mat[i,j])
        outf.write('\n')
    outf.close()                


def write_outf_3header(res_atom,dist_mat,out):
    ### write output file as csv file ###
    ### Three index columns chain, resid, atomname ###
    header = []
    num_atom = len(res_atom)
    for i in range(num_atom):
        header.append('-'.join(res_atom[i]))
    outf = open(out,'w')
    outf.write('chain,resid,atom')
    for i in range(num_atom):
        outf.write(',%s'%header[i])
    outf.write('\n')
    for i in range(num_atom):
        outf.write('%s'%res_atom[i][0])
        outf.write(',%s'%res_atom[i][1])
        outf.write(',%s'%res_atom[i][2])
        for j in range(num_atom):
            if j == i:
                outf.write(', 0.0')
            elif j < i:
                outf.write(', %f'%dist_mat[j,i])
            else:
                outf.write(', %f'%dist_mat[i,j])
        outf.write('\n')
    outf.close()                


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prepare template PDB files, write input files, save output PDB files in working directory')
    parser.add_argument('-pdb', dest='pdb', default=None, help='opt_pdb')

    args = parser.parse_args()
    pdb = args.pdb
    out = pdb.split('.')[0]+'-distmat.csv'
    

    res_atom, xyz_atom = get_res_atom(pdb) 
    num_atom = len(xyz_atom)
    dist_mat = np.zeros((num_atom, num_atom))
    for i in range(num_atom):
        for j in range(i+1,num_atom):
            if res_atom[i][:2] != res_atom[j][:2]:
                dist_mat[i,j] = calc_dist(xyz_atom[i],xyz_atom[j])
    
#    ### write output file as csv file ###
#    write_outf_1header(res_atom,dist_mat,out)
    write_outf_2header(res_atom,dist_mat,out)
