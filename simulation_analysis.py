import subprocess as sp
import os
import pandas as pd
import glob
from math import sqrt, pi
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from matplotlib.ticker import AutoMinorLocator
from itertools import cycle
from warnings import simplefilter
import psutil
matplotlib.use('Agg')
simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
plt.rcParams.update({'figure.max_open_warning': 0})

##############################################
# Highly Automatic Analysis of MD Simulation #
# for the hERG Channel                       #
##############################################
# Plot the following data:
# 1. Changes in area of filter (based on CA or O distances)
# 2. Changes in area of pore (based on CA or O distances)
# 3. Symmetry of Filter
# 4. Symmetry of Pore
# 5. Phi/Psi/Chi1/Chi2 angles of filter and pore residues
# 6. Numbers of water in filter and pore
#
# Requirement:
# VMD needs to be set so that it can be opened by invoking "vmd" in the terminal.
#
# Usage:
# 1. Place this script in directory containing simulation files.
# 2. Run this script (Python 3):
#   python3 simulation_analysis.py  -p [protein structure file (default=.psf in current dir)]
#                                   -d [simulation trajectory file (default=.dcd in current dir)]
#                                   -t [total simulation time of the whole trajectory (default = 1000)]
#                                   -e [analyze the simulation to this frame (default = -1 i.e. all)]
#                                   -s [split each plot into # of smaller plots covering different time periods (default = 1 i.e. do not split)]
#                                   --drug [set to the segname of the drug to analyze drug movement in the pore as opposed to ion movement in the filter]
#   Optional arguments:
#                                   -x [x label (default='Time (ns)')]
#                                   --runcommand [True/False, run VMD commands to generate input data, change if there is no need to calculate data again (default=True)]
#                                   --labelsize [label size (default=20)]

# Configurations
PDB_column_width_format = [(0, 4), (4, 11), (11, 16), (16, 20), (20, 22), (22, 26), (26, 38), (38, 46), (46, 54), (54, 60), (60, 66), (66, 90)]
parser = argparse.ArgumentParser(description='Construct intermolecular time series contact map')
parser.add_argument('-p', '--psf',
                    default=glob.glob('*.psf')[0],
                    dest='psf', action='store',
                    help='.psf file containing protein structural information')
parser.add_argument('-d', '--dcd',
                    default=glob.glob('*.dcd')[0],
                    dest='dcd', action='store',
                    help='.dcd file containing simulation trajectory (any trajectory format will also work)')
parser.add_argument('-t', '--time',
                    default=5000,
                    dest='time_total', action='store', type=float,
                    help='total simulation time of the full provided trajectory')
parser.add_argument('-e', '--end',
                    default=-1,
                    dest='end_frame', action='store', type=int,
                    help='analyze to which frame of the simulation')
parser.add_argument('--drug',
                    default=None,
                    dest='drug', action='store',
                    help='segname of the drug to analyze drug movement in the pore, if enabled (by inputting drug segname) will replace ion movement with drug movement')
parser.add_argument('-s', '--split',
                    default=1,
                    dest='split', action='store', type=int,
                    help='split each plot into # of smaller plots covering different time periods, useful for long simulations')
parser.add_argument('-x', '--xlabel',
                    default='Time (ns)',
                    dest='x_label', action='store',
                    help='x label')
parser.add_argument('--runcommand',
                    default=True,
                    dest='run_command', action='store',
                    help='run VMD commands to generate input data, only set to False if the script has already been ran at least once')
parser.add_argument('--labelsize',
                    default=20,
                    dest='size', action='store', type=float,
                    help='label font size (default = 20)')
arg = parser.parse_args()


def euclidean_distance(xyz_A, xyz_B):
    """Calculate Euclidean distance between point A and point B given x, y, z coordinates"""
    return sqrt(sum([(i - j)**2 for i, j in zip(xyz_A, xyz_B)]))


def closest(input_list, k):
    """Find closest number to k in list"""
    input_list = np.asarray(input_list)
    index = (np.abs(input_list - k)).argmin()
    return index, input_list[index]


def normalize_angles(angles, ax):
    """If there are too many angles values satisfying criteria below, add 360 to negative values"""
    if ((angles < -140).sum() + (angles > 140).sum()) > ((angles > -40).sum() & (angles < 40).sum()):
        ax.set(ylim=(0, 360))
        return np.where(angles < 0, angles + 360, angles)
    else:
        ax.set(ylim=(-180, 180))
        return angles


def split_data(data_set, current_part, num_parts=arg.split):
    """Split a large [data_set] into [num_parts], only consider data from the [current_part] (starting index = 0) for this run of the script"""
    if current_part >= num_parts:
        print('ERROR: When splitting data, [current_part] must be < [num_parts]. Note: for [current_part], 0 is the first part.')
        exit(1)
    return np.array_split(data_set, num_parts)[current_part]


def turn_y_axis_symmetric(ax):
    """Turn y axis symmetric by making the distance from 1 to max y limit same as the distance from 1 to min y limit"""
    max = ax.get_ylim()[1]
    min = ax.get_ylim()[0]
    if max - 1 > 1 - min:
        ax.set(ylim=(1 - (max - 1), max))
    else:
        ax.set(ylim=(min, 1 + (1 - min)))


# Define filter and pore labels
filter_labels = ['S624', 'V625', 'G626', 'F627', 'G628']
pore_labels = ['Y652', 'F656', 'S660']

# Conversion of 3-letter amino acid code to 1-letter code
aa_names = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
            'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
            'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
            'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
            'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

# Automatically determine input file name if given wildcard as input - will take first result that appears as input
if arg.psf.split('.')[0] == '*':
    arg.psf = glob.glob('*.' + arg.psf.split('.')[-1])[0]
if arg.dcd.split('.')[0] == '*':
    arg.dcd = glob.glob('*.' + arg.dcd.split('.')[-1])[0]

# Print input information
print('PSF          :', arg.psf)
print('DCD          :', arg.dcd)
print('Filter Res   :', filter_labels)
print('Pore Res     :', pore_labels)
name = arg.dcd.split('.')[0]

if arg.drug:
    print('Drug Seg     :', arg.drug, end='\n\n')

# Create folders to store data and output
os.makedirs('temp_pdb', exist_ok=True)
os.makedirs('saved_results', exist_ok=True)
os.makedirs('figures', exist_ok=True)

# Check required files
for script in ['calculate_dihedrals.tcl', 'dihedral_angles_atom_names.tcl', 'pore_water.tcl', 'prot_center.tcl', 'drug_movement.tcl']:
    if not os.path.exists(script):
        print('ERROR: Required script', script, 'not found in current directory.')
        exit(1)

# Load simulation trajectory and extract data
vmd_cmd_file = arg.dcd + '_vmd_cmd.tcl'

if arg.run_command == True:
    with open(vmd_cmd_file, 'w+') as f:
        # Load trajectory files
        f.write('mol new ' + arg.psf + ' type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n')
        f.write('mol addfile ' + arg.dcd + ' type ' + arg.dcd.split('.')[-1] + ' first 0 last ' + str(arg.end_frame) + ' step 1 filebonds 1 autobonds 1 waitfor all\n')

        # Center protein
        f.write('source prot_center.tcl\n')

        # Obtain S0, S4 and get water in pore and filter (this section needs to be first to analyze first frame)
        f.write('set PHE_O [atomselect top "segname PROA and sequence SVGF and resname PHE and name O" frame 0]\n')
        f.write('set PHE_O_z [$PHE_O get z]\n')
        f.write('set GLY_O [atomselect top "segname PROA and sequence SVGF and resname GLY and name O" frame 0]\n')
        f.write('set GLY_O_z [$GLY_O get z]\n')
        f.write('set S1 [expr {($PHE_O_z + $GLY_O_z) / 2}]\n')

        f.write('set SER_O [atomselect top "segname PROA and sequence SVGFG and resname SER and name O" frame 0]\n')
        f.write('set SER_O_z [$SER_O get z]\n')
        f.write('set SER_OG [atomselect top "segname PROA and sequence SVGFG and resname SER and name OG" frame 0]\n')
        f.write('set SER_OG_z [$SER_OG get z]\n')
        f.write('set S4 [expr {($SER_O_z + $SER_OG_z) / 2}]\n')
        #
        f.write('set LYS_CA [atomselect top "segname PROA and sequence QRL and resname LEU and name CA" frame 0]\n')
        f.write('set bottomPore [$LYS_CA get z]\n')
        f.write('set bottomPoreLocation [open "' + os.path.join('temp_pdb', name + '_bottom_pore.dat') + '" w]\n')
        f.write('puts $bottomPoreLocation "$bottomPore"\n')
        f.write('close $bottomPoreLocation\n')

        f.write('set outputname "' + os.path.join('temp_pdb', name+ '_water.dat') + '"\n')
        f.write('source pore_water.tcl\n')

        #  Obtain coordinates of atoms in filter and pore
        # f.write('set sel [atomselect top "(sequence SVGFG and name CA) or (sequence SVGFG and name O) or (sequence YASIFGNVS and resname TYR PHE SER and name CA)"]\n')
        f.write('set sel [atomselect top "(sequence SVGFG and name CA) or (sequence SVGFG and name O) or (((sequence YAS and resname TYR) or (sequence IFGNVS and resname PHE SER)) and name CA)"]\n')
        f.write('animate write pdb ' + os.path.join('temp_pdb', name + '_filter_and_pore_area.pdb') + ' skip 1 sel $sel\n')

        # Obtain S0-S4 K+ binding positions in filter
        f.write('set sel [atomselect top "segname PROA and sequence SVGFG and oxygen"]\n')
        f.write('animate write pdb ' + os.path.join('temp_pdb', name + '_S0-S4.pdb') + ' beg 0 end 0 skip 1 sel $sel\n')

        # Calculate phi/psi/chi1/chi2 angles
        f.write('source calculate_dihedrals.tcl\n')
        for segname in ['PROA', 'PROB', 'PROC', 'PROD']:
            f.write('print_sidechain_dihedrals ' + segname + ' "sequence SVGFG" top "' + os.path.join('temp_pdb', name + '_' + segname + '_filter_angles.dat') + '"\n')
            f.write('print_sidechain_dihedrals ' + segname + ' "sequence YASIFGNVS and resname TYR" top "' + os.path.join('temp_pdb', name + '_' + segname + '_poreY_angles.dat') + '"\n')
            f.write('print_sidechain_dihedrals ' + segname + ' "sequence IFGNVS and resname PHE SER" top "' + os.path.join('temp_pdb', name + '_' + segname + '_poreFS_angles.dat') + '"\n')

        # If applicable, get drug coordinates, else get K+ coordinates
        if arg.drug:
            f.write('set drug "' + arg.drug + '"\n')
            f.write('set outputname_max_z "' + os.path.join('temp_pdb', name + '_' + arg.drug + '_max_z.dat') + '"\n')
            f.write('set outputname_com_z "' + os.path.join('temp_pdb', name + '_' + arg.drug + '_com_z.dat') + '"\n')
            f.write('set outputname_min_z "' + os.path.join('temp_pdb', name + '_' + arg.drug + '_min_z.dat') + '"\n')
            f.write('source drug_movement.tcl\n')
        else:
            f.write('set sel [atomselect top "name POT"]\n')
            f.write('animate write pdb ' + os.path.join('temp_pdb', name + '_potassium_ions.pdb') + ' skip 1 sel $sel\n')

        f.write('exit')

    sp.call(['/bin/bash', '-i', '-c', 'vmd -dispdev text -e ' + vmd_cmd_file], stdin=sp.PIPE)

# Obtain number of atoms and atom list
num_atoms = 0
with open(os.path.join('temp_pdb', name + '_filter_and_pore_area.pdb')) as f:
    next(f)  # Skip header
    for line in f:
        if 'END' in line:
            break
        num_atoms += 1
atoms_list = pd.read_fwf(os.path.join('temp_pdb', name + '_filter_and_pore_area.pdb'), header=None, colspecs=PDB_column_width_format, skiprows=1, nrows=num_atoms)
num_atoms_in_a_segment = len(atoms_list[11]) // 4

# Set atom labels
atom_labels = []
opposing_atom_labels = []
for i in filter_labels:
    atom_labels.append(i + ' CA')
    atom_labels.append(i + ' O')
    opposing_atom_labels.append(i + ' CA{AC/BD}$')
    opposing_atom_labels.append(i + ' O{AC/BD}$')
for i in pore_labels:
    atom_labels.append(i + ' CA')
    opposing_atom_labels.append(i + ' O')

# for i in range(num_atoms_in_a_segment):
#     atom_labels.append(atoms_list[3][i] + (atoms_list[5][i] + 401).astype('str') + ' ' + atoms_list[2][i])
#     opposing_atom_labels.append(atoms_list[3][i] + (atoms_list[5][i] + 401).astype('str') + ' ' + atoms_list[2][i] + '(AC/BD)')

# Use data at initial frame to check and make sure obtained residues are correct and match with pre-defined labels before proceeding
for label, ref_residue in zip(atom_labels, atoms_list[3].values[:num_atoms_in_a_segment]):
    if label[0] != aa_names.get(ref_residue):
        print('ERROR: Data for filter or pore residues yield residue', aa_names.get(ref_residue), 'while', label[0], 'was expected in the provided labels. Check vmd_cmd_file section in this script.')
        exit(1)

# Obtain number of ions and initialize dataframe to store their coordinates
if not arg.drug:
    num_ions = 0
    with open(os.path.join('temp_pdb', name + '_potassium_ions.pdb')) as f:
        next(f)  # Skip header
        for line in f:
            if 'END' in line:
                break
            num_ions += 1
    ions_list = pd.read_fwf(os.path.join('temp_pdb', name + '_potassium_ions.pdb'), header=None, colspecs=PDB_column_width_format, skiprows=1, nrows=num_ions)
    ions_list = ions_list[ions_list[2] == 'POT']
    ions_list = (ions_list[2] + ions_list[1].astype('str')).values
    ions_x = pd.DataFrame(index=ions_list)
    ions_y = pd.DataFrame(index=ions_list)
    ions_z = pd.DataFrame(index=ions_list)

################################################################################################

# Obtain positions of K+ binding in the SF (S0, S1, S2, S3, S4) as well as location of the bottom of the pore (defined as "sequence QRL and resname LEU and name CA")
S_positions = pd.read_fwf(os.path.join('temp_pdb', name + '_S0-S4.pdb'), header=None, colspecs=PDB_column_width_format, skiprows=1)
S_positions.drop(S_positions.index[[-1]], inplace=True)
S0 = round((S_positions[8][5] + S_positions[8][4]) / 2, 3)
S1 = round((S_positions[8][4] + S_positions[8][3]) / 2, 3)
S2 = round((S_positions[8][3] + S_positions[8][2]) / 2, 3)
S3 = round((S_positions[8][2] + S_positions[8][1]) / 2, 3)
S4 = round((S_positions[8][1] + S_positions[8][0]) / 2, 3)
bottom_pore = float(open(os.path.join('temp_pdb', name + '_bottom_pore.dat'), 'r').readline())

################################################################################################

# Analyze distances between atoms in the filter and in the pore to calculate areas and symmetry ratios
frame_count = 0
area_list = []
symmetry_list = []
AC_list = []
BD_list = []
saved_results = [os.path.join('saved_results', name + '_filter_and_pore_area.csv'),
                 os.path.join('saved_results', name + '_symmetry.csv'),
                 os.path.join('saved_results', name + '_AC.csv'),
                 os.path.join('saved_results', name + '_BD.csv')]
print('>> Analyzing distances between atoms in the filter and in the pore...')
if not all(os.path.exists(result) for result in saved_results):
    for frame in pd.read_fwf(os.path.join('temp_pdb', name + '_filter_and_pore_area.pdb'), chunksize=num_atoms + 1, header=None, colspecs=PDB_column_width_format, skiprows=1):
        # Delete last line (containing just 'END'), reset index to 0, and name columns
        frame.drop(frame.index[[-1]], inplace=True)
        frame.reset_index(drop=True, inplace=True)
        frame.columns = ['Atom', 'Index', 'Type', 'Residue', 'Chain', 'ResNum', 'X', 'Y', 'Z', 'Occupancy', 'TempFactor', 'Segment']
        # Calculate distances and area formed by atoms in the filter with those in the opposing subunit
        filter_and_pore_area = []
        symmetry_ratio = []
        AC_values = []
        BD_values = []
        for i in range(num_atoms_in_a_segment):
            # print('AC', frame.loc[i, ['Type', 'Residue', 'Segment']].values, frame.loc[i + num_atoms_in_a_segment * 2, ['Type', 'Residue', 'Segment']].values)
            # print('BD', frame.loc[i + num_atoms_in_a_segment, ['Type', 'Residue', 'Segment']].values,  frame.loc[i + num_atoms_in_a_segment * 3, ['Type', 'Residue', 'Segment']].values)
            AC = euclidean_distance(frame.loc[i, ['X', 'Y', 'Z']].values, frame.loc[i + num_atoms_in_a_segment * 2, ['X', 'Y', 'Z']].values)
            BD = euclidean_distance(frame.loc[i + num_atoms_in_a_segment, ['X', 'Y', 'Z']].values, frame.loc[i + num_atoms_in_a_segment * 3, ['X', 'Y', 'Z']].values)
            filter_and_pore_area.append(pi * AC / 2 * BD / 2)
            symmetry_ratio.append(AC / BD)
            AC_values.append(AC)
            BD_values.append(BD)
        area_list.append(filter_and_pore_area)
        symmetry_list.append(symmetry_ratio)
        AC_list.append(AC_values)
        BD_list.append(BD_values)
        frame_count += 1
        if frame_count % 50 == 0:
            print('\rPROGRESS:     ', frame_count, end=' frames', flush=True)
    print('\rPROGRESS:     ', frame_count, end=' frames (DONE)\n', flush=True)
    # Concatenate to final results
    area = pd.DataFrame(np.transpose(area_list), index=atom_labels)
    symmetry = pd.DataFrame(np.transpose(symmetry_list), index=opposing_atom_labels)
    ACs = pd.DataFrame(np.transpose(AC_list), index=atom_labels)
    BDs = pd.DataFrame(np.transpose(BD_list), index=atom_labels)
    # Save results as csv
    area.to_csv(saved_results[0])
    symmetry.to_csv(saved_results[1])
    ACs.to_csv(saved_results[2])
    BDs.to_csv(saved_results[3])
    print('   Saved atom distances to', saved_results[0], 'at frame', frame_count)
    print('   Saved symmetry ratios to', saved_results[1])
    print('   Saved distances from subunits A-C to', saved_results[2])
    print('   Saved distances from subunits B-D to', saved_results[3])
else:
    print('   Previously saved results found. Skipping calculations!')
    area = pd.read_csv(saved_results[0], index_col=0)
    symmetry = pd.read_csv(saved_results[1], index_col=0)
    ACs = pd.read_csv(saved_results[2], index_col=0)
    BDs = pd.read_csv(saved_results[3], index_col=0)
print('   Memory usage: ', round(psutil.Process(os.getpid()).memory_percent(), 2), '%', sep='')
################################################################################################

# Analyze ion movement through the filter
frame_count = 0
if not arg.drug:
    saved_results = [os.path.join('saved_results', name + '_ions_x.csv'),
                     os.path.join('saved_results', name + '_ions_y.csv'),
                     os.path.join('saved_results', name + '_ions_z.csv')]
    print('>> Analyzing ion movement through the filter...')
    if not all(os.path.exists(result) for result in saved_results):
        for frame in pd.read_fwf(os.path.join('temp_pdb', name + '_potassium_ions.pdb'), chunksize=num_ions + 1, header=None, colspecs=PDB_column_width_format, skiprows=1):
            # Delete last line (containing just 'END'), reset index to 0, and name columns
            frame.drop(frame.index[[-1]], inplace=True)
            frame.reset_index(drop=True, inplace=True)
            frame.columns = ['Atom', 'Index', 'Type', 'Residue', 'Chain', 'ResNum', 'X', 'Y', 'Z', 'Occupancy', 'TempFactor', 'Segment']
            # Drop every row except those containing ion of interest
            frame = frame[frame['Residue'] == 'POT']
            # Record z values of each ion of interest to dataframe containing the results so far
            ions_x[frame_count] = frame['X'].values
            ions_y[frame_count] = frame['Y'].values
            ions_z[frame_count] = frame['Z'].values
            frame_count += 1
            if frame_count % 50 == 0:
                print('\rPROGRESS:     ', frame_count, end=' frames', flush=True)
        print('\rPROGRESS:     ', frame_count, end=' frames (DONE)\n', flush=True)
        # Save results as csv
        ions_x.to_csv(saved_results[0])
        ions_y.to_csv(saved_results[1])
        ions_z.to_csv(saved_results[2])
        print('   Saved ion x coordinates to', saved_results[0], 'at frame', frame_count)
        print('   Saved ion y coordinates to', saved_results[1])
        print('   Saved ion z coordinates to', saved_results[2])
    else:
        print('   Previously saved results found. Skipping calculations!')
        ions_x = pd.read_csv(saved_results[0], index_col=0)
        ions_y = pd.read_csv(saved_results[1], index_col=0)
        ions_z = pd.read_csv(saved_results[2], index_col=0)

################################################################################################

# Convert frame count to simulation time
print('>> Converting frame count to simulation time...')
datasets_list = [area, symmetry, ACs, BDs]
if not arg.drug:
    datasets_list += [ions_x, ions_y, ions_z]
for datasets in datasets_list:
    datasets.columns = [float(int(i) * arg.time_total / len(datasets.columns)) for i in datasets.columns]
# Set up time course
time = area.columns.values

################################################################################################

# Set up plots
sns.set_context('talk')
title_size = arg.size
label_size = arg.size
fig1, axes1 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filterCA_diameters', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig2, axes2 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filterO_diameters', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig3, axes3 = plt.subplots(4, 1, sharex='col', figsize=(19, 14), num='poreCA_diameters', gridspec_kw={'height_ratios': [1, 1, 1, 2]})
fig4, axes4 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filterCA_symmetry', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig5, axes5 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filterO_symmetry', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig6, axes6 = plt.subplots(4, 1, sharex='col', figsize=(19, 14), num='poreCA_symmetry', gridspec_kw={'height_ratios': [1, 1, 1, 2]})  # rows, columns
fig7, axes7 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filter_phi_AC', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig8, axes8 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filter_psi_AC', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig9, axes9 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filter_chi1_AC', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig10, axes10 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filter_chi2_AC', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig11, axes11 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filter_phi_BD', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig12, axes12 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filter_psi_BD', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig13, axes13 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filter_chi1_BD', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig14, axes14 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filter_chi2_BD', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig15, axes15 = plt.subplots(4, 1, sharex='col', figsize=(19, 14), num='pore_phi_AC', gridspec_kw={'height_ratios': [1, 1, 1, 2]})
fig16, axes16 = plt.subplots(4, 1, sharex='col', figsize=(19, 14), num='pore_psi_AC', gridspec_kw={'height_ratios': [1, 1, 1, 2]})
fig17, axes17 = plt.subplots(4, 1, sharex='col', figsize=(19, 14), num='pore_chi1_AC', gridspec_kw={'height_ratios': [1, 1, 1, 2]})
fig18, axes18 = plt.subplots(4, 1, sharex='col', figsize=(19, 14), num='pore_chi2_AC', gridspec_kw={'height_ratios': [1, 1, 1, 2]})
fig19, axes19 = plt.subplots(4, 1, sharex='col', figsize=(19, 14), num='pore_phi_BD', gridspec_kw={'height_ratios': [1, 1, 1, 2]})
fig20, axes20 = plt.subplots(4, 1, sharex='col', figsize=(19, 14), num='pore_psi_BD', gridspec_kw={'height_ratios': [1, 1, 1, 2]})
fig21, axes21 = plt.subplots(4, 1, sharex='col', figsize=(19, 14), num='pore_chi1_BD', gridspec_kw={'height_ratios': [1, 1, 1, 2]})
fig22, axes22 = plt.subplots(4, 1, sharex='col', figsize=(19, 14), num='pore_chi2_BD', gridspec_kw={'height_ratios': [1, 1, 1, 2]})
fig23, axes23 = plt.subplots(7, 1, sharex='col', figsize=(19, 23.3), num='water_filter_pore', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 1, 2]})
axes_list = [axes1, axes2, axes3, axes4, axes5, axes6, axes7, axes8, axes9, axes10, axes11, axes12, axes13, axes14, axes15, axes16, axes17, axes18, axes19, axes20, axes21, axes22, axes23]

print('   Memory usage: ', round(psutil.Process(os.getpid()).memory_percent(), 2), '%', sep='')

################################################################################################

# Change index order so filter residues would be on top while pore residues would be on the bottom
index_order = [8, 6, 4, 2, 0, 9, 7, 5, 3, 1, 10, 11, 12]  # first four are filter CA, next four are filter O, last three are pore CA
area = area.reindex(index=[area.index[i] for i in index_order])
symmetry = symmetry.reindex(index=[symmetry.index[i] for i in index_order])
ACs = ACs.reindex(index=[ACs.index[i] for i in index_order])
BDs = BDs.reindex(index=[BDs.index[i] for i in index_order])
# changes_in_area = area.sub(area.iloc[:, 0], axis=0)

# Plot areas and symmetry of filter and pore over time
print('>> Plotting areas and symmetry of filter and pore over time...')
for row, label in zip(range(5),  reversed(filter_labels)):
    # Areas
    # sns.lineplot(data=area.iloc[row], ax=axes1[row])
    # sns.lineplot(x=[time[0], time[-1]], y=[area.iloc[row][0]] * 2, color='navy', alpha=0.7, linestyle='--', ax=axes1[row])  # horizontal line to depict initial area
    axes1[row].stackplot(time, ACs.iloc[row].values, BDs.iloc[row].values, labels=['A-C', 'B-D'], colors=['#F3D3BD', '#5E5E5E'])
    axes1[0].legend(loc='lower right', ncol=2)
    axes1[row].set_ylabel(label + ' $C\\alpha$', fontsize=label_size)
    # Symmetry Ratios
    sns.lineplot(data=symmetry.iloc[row], ax=axes4[row], color='black')
    sns.lineplot(x=[time[0], time[-1]], y=[1] * 2, color='gray', alpha=0.7, linestyle='--', ax=axes4[row])  # horizontal line to depict perfect symmetry (1)
    axes4[row].set_ylabel(label + ' $C\\alpha_{AC/BD}$', fontsize=label_size)
    turn_y_axis_symmetric(axes4[row])
axes1[0].set_title('Combined Diameters of Filter Residues, $C\\alpha$ - ' + name, y=1.04, fontsize=title_size)
axes4[0].set_title('Ratios of Symmetry between Filter Residues, $C\\alpha$  - ' + name, y=1.04, fontsize=title_size)

for row, label in zip(range(5, 10),  reversed(filter_labels)):
    # Areas
    # sns.lineplot(data=area.iloc[row], ax=axes2[row - 5], color='red')
    # sns.lineplot(x=[time[0], time[-1]], y=[area.iloc[row][0]] * 2, color='firebrick', alpha=0.7, linestyle='--', ax=axes2[row - 5])  # create horizontal line to depict initial area
    axes2[row - 5].stackplot(time, ACs.iloc[row].values, BDs.iloc[row].values, labels=['A-C', 'B-D'], colors=['#A2A7A5', '#6D696A'])
    axes2[0].legend(loc='lower right', ncol=2)
    axes2[row - 5].set_ylabel(label + ' $O$', fontsize=label_size)
    # Symmetry Ratios
    sns.lineplot(data=symmetry.iloc[row], ax=axes5[row - 5], color='black')
    sns.lineplot(x=[time[0], time[-1]], y=[1] * 2, color='gray', alpha=0.7, linestyle='--', ax=axes5[row - 5])  # horizontal line to depict perfect symmetry (1)
    axes5[row - 5].set_ylabel(label + ' $O_{AC/BD}$', fontsize=label_size)
    turn_y_axis_symmetric(axes5[row - 5])
axes2[0].set_title('Combined Diameters of Filter Residues, $O$ - ' + name, y=1.04, fontsize=title_size)
axes5[0].set_title('Ratios of Symmetry between Filter Residues, $O$ - ' + name, y=1.04, fontsize=title_size)

for row, label in zip(range(10, 13), pore_labels):
    # Areas
    # sns.lineplot(data=area.iloc[row], ax=axes3[row - 10])
    # sns.lineplot(x=[time[0], time[-1]], y=[area.iloc[row][0]] * 2, color='navy', alpha=0.7, linestyle='--', ax=axes3[row - 10])  # create horizontal line to depict initial area
    axes3[row - 10].stackplot(time, ACs.iloc[row].values, BDs.iloc[row].values, labels=['A-C', 'B-D'], colors=['#B1E5F2', '#272635'])
    axes3[0].legend(loc='lower right', ncol=2)
    axes3[row - 10].set_ylabel(label + ' $C\\alpha$', fontsize=label_size)
    # Symmetry Ratios
    sns.lineplot(data=symmetry.iloc[row], ax=axes6[row - 10], color='black')
    sns.lineplot(x=[time[0], time[-1]], y=[1] * 2, color='seagreen', alpha=0.7, linestyle='--', ax=axes6[row - 10])  # horizontal line to depict perfect symmetry (1)
    axes6[row - 10].set_ylabel(label + ' $C\\alpha_{AC/BD}$', fontsize=label_size)
    turn_y_axis_symmetric(axes6[row - 10])
axes3[0].set_title('Combined Diameters of Pore Residues, $C\\alpha$  - ' + name, y=1.04, fontsize=title_size)
axes6[0].set_title('Ratios of Symmetry between Pore Residues, $C\\alpha$  - ' + name, y=1.04, fontsize=title_size)
print('   Memory usage: ', round(psutil.Process(os.getpid()).memory_percent(), 2), '%', sep='')

################################################################################################

# Plot phi/psi, chi1/chi2
print('>> Plotting phi, psi, chi1, chi2 angles...')
# Selectivity filter residues
angles_data_column_width_format = [(0, 5), (5, 12), (12, 17), (17, 23), (23, 29), (29, 35), (35, 40)]
for segname in ['PROA', 'PROB', 'PROC', 'PROD']:
    # Obtain angles data
    angles_data = pd.read_fwf(os.path.join('temp_pdb', name + '_' + segname + '_filter_angles.dat'), header=None, colspecs=angles_data_column_width_format)
    angles_data.columns = ['ResNum', 'Residue', 'Frame', 'Phi', 'Psi', 'Chi1', 'Chi2']
    # Use data at initial frame to check and make sure obtained residues are correct and match with pre-defined labels before proceeding
    filter_resnum = angles_data.loc[angles_data['Frame'] == 0, 'ResNum'].values
    for label, ref_residue in zip(filter_labels, angles_data.loc[angles_data['Frame'] == 0, 'Residue'].values):
        if label[0] != aa_names.get(ref_residue):
            print('ERROR: Data for filter angles yield residue', aa_names.get(ref_residue), 'while', label[0], 'was expected in the provided labels. Check vmd_cmd_file section in this script.')
            exit(1)
    # Plot angles, reversed to match visual orientation when viewed by Chimera or VMD
    for resnum, row, label in zip(reversed(filter_resnum), range(5), reversed(filter_labels)):
        if segname in ['PROA', 'PROC']:
            for angle_type, ax in zip(['Phi', 'Psi', 'Chi1', 'Chi2'], [axes7, axes8, axes9, axes10]):
                if segname == 'PROA':
                    axe = ax[row]
                    axe.set_ylabel(label + '$_A$', fontsize=label_size)
                    color = 'orangered'
                    axe.yaxis.label.set_color(color)
                    axe.tick_params(axis='y', colors=color)
                else:
                    axe = ax[row].twinx()
                    axe.set_ylabel(label + '$_C$', fontsize=label_size, rotation=270, labelpad=20)
                    color = 'steelblue'
                    axe.yaxis.label.set_color(color)
                    axe.tick_params(axis='y', colors=color)
                # If dataset is filled with NaN, skip
                if not np.isnan(angles_data.loc[angles_data['ResNum'] == resnum, angle_type].values).all():
                    sns.lineplot(x=time, y=normalize_angles(angles_data.loc[angles_data['ResNum'] == resnum, angle_type].values, axe), color=color, ax=axe)
        else:
            for angle_type, ax in zip(['Phi', 'Psi', 'Chi1', 'Chi2'], [axes11, axes12, axes13, axes14]):
                if segname == 'PROB':
                    axe = ax[row]
                    axe.set_ylabel(label + '$_B$', fontsize=label_size)
                    color = 'orangered'
                    axe.yaxis.label.set_color(color)
                    axe.tick_params(axis='y', colors=color)
                else:
                    axe = ax[row].twinx()
                    axe.set_ylabel(label + '$_D$', fontsize=label_size, rotation=270, labelpad=20)
                    color = 'steelblue'
                    axe.yaxis.label.set_color(color)
                    axe.tick_params(axis='y', colors=color)
                # If dataset is filled with NaN, skip
                if not np.isnan(angles_data.loc[angles_data['ResNum'] == resnum, angle_type].values).all():
                    sns.lineplot(x=time, y=normalize_angles(angles_data.loc[angles_data['ResNum'] == resnum, angle_type].values, axe), color=color, ax=axe)
# Set titles
axes7[0].set_title('Filter Phi Angles - ' + name, y=1.04, fontsize=title_size)
axes8[0].set_title('Filter Psi Angles - ' + name, y=1.04, fontsize=title_size)
axes9[0].set_title('Filter Chi1 Angles - ' + name, y=1.04, fontsize=title_size)
axes10[0].set_title('Filter Chi2 Angles - ' + name, y=1.04, fontsize=title_size)
axes11[0].set_title('Filter Phi Angles - ' + name, y=1.04, fontsize=title_size)
axes12[0].set_title('Filter Psi Angles - ' + name, y=1.04, fontsize=title_size)
axes13[0].set_title('Filter Chi1 Angles - ' + name, y=1.04, fontsize=title_size)
axes14[0].set_title('Filter Chi2 Angles - ' + name, y=1.04, fontsize=title_size)

# Pore residues
for segname in ['PROA', 'PROB', 'PROC', 'PROD']:
    # Obtain angles data and concatenate angles data
    Y_angles_data = pd.read_fwf(os.path.join('temp_pdb', name + '_' + segname + '_poreY_angles.dat'), header=None, colspecs=angles_data_column_width_format)
    FS_angles_data = pd.read_fwf(os.path.join('temp_pdb', name + '_' + segname + '_poreFS_angles.dat'), header=None, colspecs=angles_data_column_width_format)
    angles_data = pd.concat([Y_angles_data, FS_angles_data], ignore_index=True)
    angles_data.columns = ['ResNum', 'Residue', 'Frame', 'Phi', 'Psi', 'Chi1', 'Chi2']
    angles_data = angles_data.sort_values(['Frame', 'ResNum'], ascending=(True, True))
    angles_data.reset_index(drop=True, inplace=True)
    # Use data at initial frame to check and make sure obtained residues are correct and match with pre-defined labels before proceeding
    pore_resnum = angles_data.loc[angles_data['Frame'] == 0, 'ResNum'].values
    for label, ref_residue in zip(pore_labels, angles_data.loc[angles_data['Frame'] == 0, 'Residue'].values):
        if label[0] != aa_names.get(ref_residue):
            print('ERROR: Data for pore angles yield residue', aa_names.get(ref_residue), 'while', label[0], 'was expected in the provided labels. Check vmd_cmd_file section in this script.')
            exit(1)
    # Plot angles
    for resnum, row, label in zip(pore_resnum, range(3), pore_labels):
        if segname in ['PROA', 'PROC']:
            for angle_type, ax in zip(['Phi', 'Psi', 'Chi1', 'Chi2'], [axes15, axes16, axes17, axes18]):
                if segname == 'PROA':
                    axe = ax[row]
                    axe.set_ylabel(label + '$_A$', fontsize=label_size)
                    color = 'orangered'
                    axe.yaxis.label.set_color(color)
                    axe.tick_params(axis='y', colors=color)
                else:
                    axe = ax[row].twinx()
                    axe.set_ylabel(label + '$_C$', fontsize=label_size, rotation=270, labelpad=20)
                    color = 'steelblue'
                    axe.yaxis.label.set_color(color)
                    axe.tick_params(axis='y', colors=color)
                # If dataset is filled with NaN, skip
                if not np.isnan(angles_data.loc[angles_data['ResNum'] == resnum, angle_type].values).all():
                    sns.lineplot(x=time, y=normalize_angles(angles_data.loc[angles_data['ResNum'] == resnum, angle_type].values, axe), color=color, ax=axe)
        else:
            for angle_type, ax in zip(['Phi', 'Psi', 'Chi1', 'Chi2'], [axes19, axes20, axes21, axes22]):
                if segname == 'PROB':
                    axe = ax[row]
                    axe.set_ylabel(label + '$_B$', fontsize=label_size)
                    color = 'orangered'
                    axe.yaxis.label.set_color(color)
                    axe.tick_params(axis='y', colors=color)
                else:
                    axe = ax[row].twinx()
                    axe.set_ylabel(label + '$_D$', fontsize=label_size, rotation=270, labelpad=20)
                    color = 'steelblue'
                    axe.yaxis.label.set_color(color)
                    axe.tick_params(axis='y', colors=color)
                # If dataset is filled with NaN, skip
                if not np.isnan(angles_data.loc[angles_data['ResNum'] == resnum, angle_type].values).all():
                    sns.lineplot(x=time, y=normalize_angles(angles_data.loc[angles_data['ResNum'] == resnum, angle_type].values, axe), color=color, ax=axe)
# Set titles
axes15[0].set_title('Pore Phi Angles - ' + name, y=1.04, fontsize=title_size)
axes16[0].set_title('Pore Psi Angles - ' + name, y=1.04, fontsize=title_size)
axes17[0].set_title('Pore Chi1 Angles - ' + name, y=1.04, fontsize=title_size)
axes18[0].set_title('Pore Chi2 Angles - ' + name, y=1.04, fontsize=title_size)
axes19[0].set_title('Pore Phi Angles - ' + name, y=1.04, fontsize=title_size)
axes20[0].set_title('Pore Psi Angles - ' + name, y=1.04, fontsize=title_size)
axes21[0].set_title('Pore Chi1 Angles - ' + name, y=1.04, fontsize=title_size)
axes22[0].set_title('Pore Chi2 Angles - ' + name, y=1.04, fontsize=title_size)
print('   Memory usage: ', round(psutil.Process(os.getpid()).memory_percent(), 2), '%', sep='')

################################################################################################

# Plot water in pore and in filter
print('>> Plotting numbers of water molecules in pore and in filter...')
water_data = pd.read_csv(os.path.join('temp_pdb', name + '_water.dat'), header=None)
water_data.columns = ['Frame', 'Inside Filter', 'Behind Filter$_{A}$', 'Behind Filter$_{B}$', 'Behind Filter$_{C}$', 'Behind Filter$_{D}$', 'Inside Pore']

for row, label in zip(range(6), water_data.columns[1:]):
    sns.lineplot(x=time, y=water_data[label].values, color='darkviolet', ax=axes23[row])
    axes23[row].set_ylabel(label, fontsize=label_size)

axes23[0].set_title('Water Count - ' + name, y=1.04, fontsize=title_size)
print('   Memory usage: ', round(psutil.Process(os.getpid()).memory_percent(), 2), '%', sep='')

################################################################################################

# Plot drug movement in the pore if applicable, else plot ion movement
palette = cycle(sns.color_palette(n_colors=8))
if arg.drug:
    # Plot drug movement in the pore if specified in the input
    print('>> Plotting drug movement in the pore...')
    max_z_data = pd.read_csv(os.path.join('temp_pdb', name + '_' + arg.drug + '_max_z.dat'), header=None, delimiter=r'\s+', index_col=0)
    com_z_data = pd.read_csv(os.path.join('temp_pdb', name + '_' + arg.drug + '_com_z.dat'), header=None, delimiter=r'\s+', index_col=0)
    min_z_data = pd.read_csv(os.path.join('temp_pdb', name + '_' + arg.drug + '_min_z.dat'), header=None, delimiter=r'\s+', index_col=0)

    for (mol, max_z), (mol, com_z), (mol, min_z) in zip(max_z_data.iterrows(), com_z_data.iterrows(), min_z_data.iterrows()):
        # If dataset is filled with NaN, skip
        if not np.isnan(com_z.values).all():
            for ax in axes_list:
                ax[-1].plot(time, com_z.values)
                ax[-1].fill_between(time, min_z.values, max_z.values, alpha=0.2)
        print('\rPROGRESS:     ', mol, end=' molecules', flush=True)
    print('\rPROGRESS:     ', mol, end=' molecules (DONE)\n', flush=True)

    # Plot reference lines
    for ax in axes_list:
        sns.lineplot(x=[time[0], time[-1]], y=[S4] * 2, color='blue', alpha=0.7, linestyle='--', ax=ax[-1])
        sns.lineplot(x=time, y=bottom_pore, color='red', alpha=0.7, linestyle='--', ax=ax[-1])
        # Annotate reference positions, position of annotation is determined automatically
        ax[-1].annotate('S4', (time[len(time) // 20], S4 + 0.3), fontsize=label_size - 4, color='blue')
        ax[-1].annotate('S4', (time[int(len(time) // 20 * 18)], S4 + 0.3), fontsize=label_size - 4, color='blue')
        ax[-1].annotate('Pore Bottom', (time[len(time) // 20], bottom_pore + 0.3), fontsize=label_size - 4, color='red')
        ax[-1].annotate('Pore Bottom', (time[int(len(time) // 20 * 17)], bottom_pore + 0.3), fontsize=label_size - 4, color='red')
else:
    # Plot ion movement over time
    print('>> Plotting ion movement over time...')
    for (ion, x), (ion, y), (ion, z) in zip(ions_x.iterrows(), ions_y.iterrows(), ions_z.iterrows()):
        z_in_range = []
        time_in_range = []
        temp_z = []
        temp_t = []
        color = next(palette)

        for x, y, z, t in zip(x.values, y.values, z.values, time):
            if (sqrt(x**2 + y**2) <= 10 and S4 <= z <= S0) or (sqrt(x**2 + y**2) <= 15 and (S4 - 8 <= z < S4 or S0 < z <= S0 + 4)):
                if len(temp_z) != -1:
                    temp_z.append(z)
                    temp_t.append(t)
            else:
                if temp_z:
                    z_in_range.append(temp_z)
                    time_in_range.append(temp_t)
                    temp_z = []
                    temp_t = []
        if temp_z:
            z_in_range.append(temp_z)
            time_in_range.append(temp_t)
        if len(z_in_range) != 0:
            for t, z in zip(time_in_range, z_in_range):
                for ax in axes_list:
                    # ax[-1].axvspan(550, 720, color='yellow', alpha=0.2)
                    sns.lineplot(x=t, y=z, ax=ax[-1])

    # Plot reference lines
    for z, label in zip([S0, S1, S2, S3, S4], ['S0', 'S1', 'S2', 'S3', 'S4']):
        for ax in axes_list:
            sns.lineplot(x=[time[0], time[-1]], y=[z] * 2, color='blue', alpha=0.7, linestyle='--', ax=ax[-1])
            # Annotate S positions, position of annotation is determined automatically
            ax[-1].annotate(label, (time[len(time) // 20], z + 0.3), fontsize=label_size - 4, color='blue')
            ax[-1].annotate(label, (time[int(len(time) // 20 * 18)], z + 0.3), fontsize=label_size - 4, color='blue')

# Set labels
for ax in axes_list:
    if arg.drug:
        ax[-1].set_ylabel('Drug $z-pos.$ in Pore (Å)', fontsize=label_size)
    else:
        ax[-1].set_ylabel('$K^+ z-pos.$ (Å)', fontsize=label_size)
    ax[-1].set_xlabel(arg.x_label, fontsize=label_size)
    for i in ax:
        i.xaxis.set_minor_locator(AutoMinorLocator())
print('   Memory usage: ', round(psutil.Process(os.getpid()).memory_percent(), 2), '%', sep='')

################################################################################################

# Set limits
# Filter CA radii
# axes1[0].set(ylim=(20, 250))
# axes1[1].set(ylim=(20, 150))
# axes1[2].set(ylim=(20, 70))
# axes1[3].set(ylim=(40, 63))
# axes1[4].set(ylim=(50, 80))

# Filter O radii
# axes2[0].set(ylim=(30, 290))
# axes2[1].set(ylim=(0, 200))
# axes2[2].set(ylim=(10, 80))
# axes2[3].set(ylim=(10, 40))
# axes2[4].set(ylim=(10, 30))

# Filter CA symmetry
# axes4[0].set(ylim=(0.6, 2.5))
# axes4[1].set(ylim=(0.6, 1.5))
# axes4[2].set(ylim=(0.6, 1.5))
# axes4[3].set(ylim=(0.6, 1.5))
# axes4[4].set(ylim=(0.6, 1.5))

# Filter O symmetry
# axes5[0].set(ylim=(0.6, 2.5))
# axes5[1].set(ylim=(0.6, 1.8))
# axes5[2].set(ylim=(0.4, 2.5))
# axes5[3].set(ylim=(0.6, 2.5))
# axes5[4].set(ylim=(0.6, 1.5))

# # Pore radii
# axes3[0].set(ylim=(140, 290))
# axes3[1].set(ylim=(170, 270))
# axes3[2].set(ylim=(260, 410))
#

# Ions
if arg.drug:
    for ax in axes_list:
        ax[-1].set(ylim=(bottom_pore - 7, S4 + 3))
else:
    for ax in axes_list:
        ax[-1].set(ylim=(S4 - 7, S0 + 3))

# Water in filter and pore
axes23[0].set(ylim=(-0.1, 5.1))
axes23[1].set(ylim=(-0.1, 10.1))
axes23[2].set(ylim=(-0.1, 10.1))
axes23[3].set(ylim=(-0.1, 10.1))
axes23[4].set(ylim=(-0.1, 10.1))
axes23[5].set(ylim=(-0.1, 150.1))

################################################################################################

# Save figures
print('\nPlots are saved as:')
for figure in ['filterCA_diameters', 'filterO_diameters', 'poreCA_diameters',
               'filterCA_symmetry', 'filterO_symmetry', 'poreCA_symmetry',
               'filter_phi_AC', 'filter_psi_AC', 'filter_chi1_AC', 'filter_chi2_AC',
               'filter_phi_BD', 'filter_psi_BD', 'filter_chi1_BD', 'filter_chi2_BD',
               'pore_phi_AC', 'pore_psi_AC', 'pore_chi1_AC', 'pore_chi2_AC',
               'pore_phi_BD', 'pore_psi_BD', 'pore_chi1_BD', 'pore_chi2_BD',
               'water_filter_pore']:
    plt.figure(figure)
    plt.savefig(os.path.join('figures', name + '_' + figure + '.png'), dpi=300)
    print(os.path.join('figures', name + '_' + figure + '.png'))
    if arg.split > 1:
        # If specified in the input, split each plot into smaller plots covering different times
        for portion in range(arg.split):
            xmin = time[-1] // arg.split * portion
            if portion < arg.split - 1:
                xmax = time[-1] // arg.split * (portion + 1)
            else:
                xmax = time[-1]
            plt.xlim(xmin, xmax)
            plt.suptitle('Part ' + str(portion + 1) + '/' + str(arg.split) + ' (' + str(int(xmin)) + ' - ' + str(int(xmax)) + ' ns)', y=0.94)
            plt.title('Part ' + str(portion + 1) + '/' + str(arg.split))
            plt.savefig(os.path.join('figures', name + '_' + figure + '_part' + str(portion + 1) + '.png'), dpi=300)
            print(os.path.join('figures', name + '_' + figure + '_part' + str(portion + 1) + '.png'))

print("\n _._     _,-'\"\"`-._\n,-.`._,'(       |\\`-/|\n    `-.-' \\ )-`( , o o)\n          `-    \\`_`\"'-")
