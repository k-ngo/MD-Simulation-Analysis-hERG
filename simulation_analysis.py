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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patheffects as pe
from functools import reduce
matplotlib.use('Agg')
simplefilter(action='ignore', category=pd.errors.PerformanceWarning)  # ignore deprecated classes and functions warnings
simplefilter(action='ignore', category=FutureWarning)
simplefilter(action='ignore', category=UserWarning)
plt.rcParams.update({'figure.max_open_warning': 0})  # ignore warnings when creating multiple figures
sns.set_context('talk')


#######################################
# Automatic Analysis of MD Simulation #
# for the hERG Channel                #
#######################################
# Plot the following data:
# 1. Changes in filter diameters (based on CA or O distances)
# 2. Changes in pore diameters (based on CA or O distances)
# 3. Symmetry of filter
# 4. Symmetry of pore
# 5. Phi/Psi/Chi1/Chi2 angles of filter and pore residues
# 6. Numbers of water in filter and pore
# 7. Ion movement in the filter
# 8. Drug movement in the pore
#
# Requirements:
# Python3, VMD
#
# Usage:
# 1. Place this script in directory containing 1 protein structure file and 1 simulation trajectory file.
# 2. Run this script (Python 3):
#   python3 simulation_analysis.py  -p [protein structure file (default=.psf in current dir)]
#                                   -d [simulation trajectory file (default=.dcd in current dir)]
#                                   -t [total simulation time of the whole trajectory (default = 1000)]
#   Optional arguments:
#                                   --drug [set to the segname of the drug to analyze drug movement in the pore as opposed to ion movement in the filter]
#                                   -e [analyze the simulation to this frame (default = -1 i.e. all)]
#                                   -s [step to read trajectory file (default = 1 i.e. do not skip any frame, 2 to skip every 1 frame)]
#                                   --split [split each plot into # of smaller plots covering different time periods (default = 1 i.e. do not split)]
#                                   -x [x label (default='Time (ns)')]
#                                   --runcommand [True/False, run VMD commands to generate input data, change if there is no need to calculate data again (default=True)]
#                                   --labelsize [label size (default=20)]

# Configurations
PDB_column_width_format = [(0, 4), (4, 11), (11, 16), (16, 20), (20, 22), (22, 26), (26, 38), (38, 46), (46, 54), (54, 60), (60, 66), (66, 90)]
parser = argparse.ArgumentParser(description='hERG Simulation Analysis')
parser.add_argument('-p', '--psf',
                    default=glob.glob('*.psf')[0],
                    dest='psf', action='store',
                    help='.psf file containing protein structural information')
parser.add_argument('-d', '--dcd',
                    default=glob.glob('*.dcd')[0],
                    dest='dcd', action='store',
                    help='.dcd file containing simulation trajectory (any trajectory format will also work)')
parser.add_argument('-t', '--time',
                    default=777,
                    dest='time_total', action='store', type=float,
                    help='total simulation time of the full provided trajectory')
parser.add_argument('-e', '--end',
                    default=-1,
                    dest='end_frame', action='store', type=int,
                    help='analyze to which frame of the simulation')
parser.add_argument('-s', '--step',
                    default=1,
                    dest='step', action='store', type=int,
                    help='step used when loading trajectory (i.e., 1 to read every single frame, 2 to read every two frames...)')
parser.add_argument('--drug',
                    default=None,
                    dest='drug', action='store',
                    help='segname of the drug to analyze drug movement in the pore, if enabled (by inputting drug segname) will replace ion movement with drug movement')
parser.add_argument('--split',
                    default=1,
                    dest='split', action='store', type=int,
                    help='split each plot into # of smaller plots covering different time periods, useful for long simulations')
parser.add_argument('-x', '--xlabel',
                    default='Time (ns)',
                    dest='x_label', action='store',
                    help='label on x-axis')
parser.add_argument('--skipcommand',
                    dest='skip_command', action='store_true',
                    help='if toggled, skip running VMD commands to generate input data, only set if the script has already been ran at least once')
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
    return index


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


def create_labels(data):
    """Create list of labels, consisting of solely whole numbers, for plotting"""
    # Necessary because by default time labels on x-axis will be ugly (e.g., 137.13, 252.23, 356.45...),
    # this section is to make it to display only labels for every specified label_interval (e.g., 100, 200, 300... if label_interval = 100) through use of rounding
    if arg.time_total <= 100:
        label_interval = 5
    elif arg.time_total <= 300:
        label_interval = 15
    elif arg.time_total <= 500:
        label_interval = 25
    elif arg.time_total <= 1000:
        label_interval = 50
    elif arg.time_total <= 3000:
        label_interval = 150
    else:
        label_interval = 250
    label_list = [np.NaN] * len(data)
    label_locations = []
    time_labels = []
    for t in np.arange(0, int(arg.time_total), label_interval):
        index_of_closest_time_in_data = closest(time, t)
        label_list[index_of_closest_time_in_data] = round(t / label_interval) * label_interval
        time_labels.append(t)
        label_locations.append(index_of_closest_time_in_data)
    return label_list, label_locations, time_labels


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
name = '.'.join(arg.dcd.split('.')[:-1])

if arg.drug:
    print('Drug Seg     :', arg.drug, end='\n\n')

# Create folders to store data and output
os.makedirs('temp_pdb_' + name, exist_ok=True)
os.makedirs('saved_results_' + name, exist_ok=True)
os.makedirs('figures_' + name, exist_ok=True)

# Check required files
for script in ['calculate_dihedrals.tcl', 'dihedral_angles_atom_names.tcl', 'pore_water.tcl', 'prot_center.tcl', 'drug_movement.tcl', 'drug_binding_residues.tcl', 'hole2']:
    if not os.path.exists(os.path.join('aux_scripts', script)):
        print('ERROR: Required script', script, 'not found in current directory.')
        exit(1)

# Load simulation trajectory and extract data
vmd_cmd_file = name + '_vmd_cmd.tcl'

if not arg.skip_command:
    with open(vmd_cmd_file, 'w+') as f:
        # Load trajectory files
        f.write('mol new ' + arg.psf + ' type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n')
        f.write('mol addfile ' + arg.dcd + ' type ' + arg.dcd.split('.')[-1] + ' first 0 last ' + str(arg.end_frame) + ' step ' + str(arg.step) + ' filebonds 1 autobonds 1 waitfor all\n')

        # Center protein
        f.write('source aux_scripts/prot_center.tcl\n')

        # Obtain initial positions of K+ binding positions S0, S4 (this section needs to be first to analyze first frame)
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

        # Obtain S0-S4 K+ binding positions in filter
        f.write('set sel [atomselect top "segname PROA and sequence SVGFG and oxygen"]\n')
        f.write('animate write pdb ' + os.path.join('temp_pdb_' + name, 'S0-S4.pdb') + ' beg 0 end -1 skip 1 sel $sel\n')

        # Obtain positions of when the pore ends based on the last LYS residue
        f.write('set LYS_CA [atomselect top "segname PROA and sequence QRL and resname LEU and name CA"]\n')
        f.write('animate write pdb ' + os.path.join('temp_pdb_' + name, 'pore_end.pdb') + ' beg 0 end -1 skip 1 sel $LYS_CA\n')
        f.write('set LYS_CA [atomselect top "segname PROA and sequence QRL and resname LEU and name CA" frame 0]\n')
        f.write('set bottomPore [$LYS_CA get z]\n')

        # Obtain positions of drug binding pore residues (TYR, PHE, SER)
        f.write('set drugBindingPoreResLocation [open "' + os.path.join('temp_pdb_' + name, 'drugBindingPoreRes.dat') + '" w]\n')
        f.write('source aux_scripts/drug_binding_residues.tcl\n')

        if arg.drug:
            f.write('set drug "' + arg.drug + '"\n')

        # Obtain water in pore and filter
        f.write('set outputname "' + os.path.join('temp_pdb_' + name, 'water.dat') + '"\n')
        f.write('source aux_scripts/pore_water.tcl\n')

        # Obtain coordinates of atoms in filter and pore
        # f.write('set sel [atomselect top "(sequence SVGFG and name CA) or (sequence SVGFG and name O) or (sequence YASIFGNVS and resname TYR PHE SER and name CA)"]\n')
        f.write('set sel [atomselect top "(sequence SVGFG and name CA) or (sequence SVGFG and name O) or (((sequence YAS and resname TYR) or (sequence IFGNVS and resname PHE SER)) and name CA)"]\n')
        f.write('animate write pdb ' + os.path.join('temp_pdb_' + name, 'filter_and_pore_area.pdb') + ' skip 1 sel $sel\n')

        # Calculate phi/psi/chi1/chi2 angles
        f.write('source aux_scripts/calculate_dihedrals.tcl\n')
        f.write('source aux_scripts/dihedral_angles_atom_names.tcl\n')
        for segname in ['PROA', 'PROB', 'PROC', 'PROD']:
            f.write('print_sidechain_dihedrals ' + segname + ' "sequence SVGFG" top "' + os.path.join('temp_pdb_' + name, segname + '_filter_angles.dat') + '"\n')
            f.write('print_sidechain_dihedrals ' + segname + ' "sequence YASIFGNVS and resname TYR" top "' + os.path.join('temp_pdb_' + name, segname + '_poreY_angles.dat') + '"\n')
            f.write('print_sidechain_dihedrals ' + segname + ' "sequence IFGNVS and resname PHE SER" top "' + os.path.join('temp_pdb_' + name, segname + '_poreFS_angles.dat') + '"\n')

        # If applicable, get drug coordinates, else get K+ coordinates
        if arg.drug:
            f.write('set outputname_max_z "' + os.path.join('temp_pdb_' + name, arg.drug + '_max_z.dat') + '"\n')
            f.write('set outputname_com_z "' + os.path.join('temp_pdb_' + name, arg.drug + '_com_z.dat') + '"\n')
            f.write('set outputname_min_z "' + os.path.join('temp_pdb_' + name, arg.drug + '_min_z.dat') + '"\n')
            f.write('source aux_scripts/drug_movement.tcl\n')
        else:
            f.write('set sel [atomselect top {name "K\\+" or name "POT"}]\n')
            f.write('animate write pdb ' + os.path.join('temp_pdb_' + name, 'potassium_ions.pdb') + ' skip 1 sel $sel\n')

        f.write('exit')

    sp.call(['/bin/bash', '-i', '-c', 'vmd -dispdev text -e ' + vmd_cmd_file], stdin=sp.PIPE)

################################################################################################

# Obtain positions of K+ binding in the SF (S0, S1, S2, S3, S4) as well as location of the bottom of the pore (defined as "sequence QRL and resname LEU and name CA")
S0 = []
S1 = []
S2 = []
S3 = []
S4 = []
S5 = []
for frame in pd.read_fwf(os.path.join('temp_pdb_' + name, 'S0-S4.pdb'), chunksize=6 + 1, header=None, colspecs=PDB_column_width_format, skiprows=1):
    # Delete last line (containing just 'END'), reset index to 0, and name columns
    frame.drop(frame.index[[-1]], inplace=True)
    frame.reset_index(drop=True, inplace=True)
    # Calculate K+ binding positions and append to lists
    S0.append(round((frame[8][5] + frame[8][4]) / 2, 3))
    S1.append(round((frame[8][4] + frame[8][3]) / 2, 3))
    S2.append(round((frame[8][3] + frame[8][2]) / 2, 3))
    S3.append(round((frame[8][2] + frame[8][1]) / 2, 3))
    S4.append(round((frame[8][1] + frame[8][0]) / 2, 3))
S0_mean = np.mean(S0)
S1_mean = np.mean(S1)
S2_mean = np.mean(S2)
S3_mean = np.mean(S3)
S4_mean = np.mean(S4)

# Obtain pore end locations
pore_end = []
for frame in pd.read_fwf(os.path.join('temp_pdb_' + name, 'pore_end.pdb'), chunksize=1 + 1, header=None, colspecs=PDB_column_width_format, skiprows=1):
    # Delete last line (containing just 'END'), reset index to 0, and name columns
    frame.drop(frame.index[[-1]], inplace=True)
    frame.reset_index(drop=True, inplace=True)
    # Calculate pore end positions and append to lists
    pore_end.append(frame[8][0])
pore_end_mean = np.mean(pore_end)

# If applicable, obtain center of mass of drug binding residues in z as well as drugs
# Drug binding residues
drug_binding_res_data = pd.read_csv(os.path.join('temp_pdb_' + name, 'drugBindingPoreRes.dat'), header=None)
drug_binding_res_data.columns = ['Frame', 'Y652', 'F656', 'S660']
# Drugs
if arg.drug:
    max_z_data = pd.read_csv(os.path.join('temp_pdb_' + name, arg.drug + '_max_z.dat'), header=None, delimiter=r'\s+', index_col=0)
    com_z_data = pd.read_csv(os.path.join('temp_pdb_' + name, arg.drug + '_com_z.dat'), header=None, delimiter=r'\s+', index_col=0)
    min_z_data = pd.read_csv(os.path.join('temp_pdb_' + name, arg.drug + '_min_z.dat'), header=None, delimiter=r'\s+', index_col=0)

# Obtain number of atoms and atom list
num_atoms = 0
with open(os.path.join('temp_pdb_' + name, 'filter_and_pore_area.pdb')) as f:
    next(f)  # Skip header
    for line in f:
        if 'END' in line:
            break
        num_atoms += 1
atoms_list = pd.read_fwf(os.path.join('temp_pdb_' + name, 'filter_and_pore_area.pdb'), header=None, colspecs=PDB_column_width_format, skiprows=1, nrows=num_atoms)
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
    with open(os.path.join('temp_pdb_' + name, 'potassium_ions.pdb')) as f:
        next(f)  # Skip header
        for line in f:
            if 'END' in line:
                break
            num_ions += 1
    ions_list = pd.read_fwf(os.path.join('temp_pdb_' + name, 'potassium_ions.pdb'), header=None, colspecs=PDB_column_width_format, skiprows=1, nrows=num_ions)
    ions_list = ions_list[(ions_list[2] == 'K+') | (ions_list[2] == 'POT')]
    ions_list = (ions_list[2] + ions_list[1].astype('str')).values
    ions_x = pd.DataFrame(index=ions_list)
    ions_y = pd.DataFrame(index=ions_list)
    ions_z = pd.DataFrame(index=ions_list)

################################################################################################

# Analyze distances between atoms in the filter and in the pore to calculate areas and symmetry ratios
frame_count = 0
area_list = []
symmetry_list = []
AC_list = []
BD_list = []
saved_results = [os.path.join('saved_results_' + name, 'filter_and_pore_area.csv'),
                 os.path.join('saved_results_' + name, 'symmetry.csv'),
                 os.path.join('saved_results_' + name, 'AC.csv'),
                 os.path.join('saved_results_' + name, 'BD.csv')]
print('>> Analyzing distances between atoms in the filter and in the pore...')
if not all(os.path.exists(result) for result in saved_results):
    for frame in pd.read_fwf(os.path.join('temp_pdb_' + name, 'filter_and_pore_area.pdb'), chunksize=num_atoms + 1, header=None, colspecs=PDB_column_width_format, skiprows=1):
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
            print('\r   PROGRESS:     ', frame_count, end=' frames', flush=True)
    print('\r   PROGRESS:     ', frame_count, end=' frames (DONE)\n', flush=True)
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
################################################################################################

# Analyze ion movement through the filter
frame_count = 0
if not arg.drug:
    saved_results = [os.path.join('saved_results_' + name, 'ions_x.csv'),
                     os.path.join('saved_results_' + name, 'ions_y.csv'),
                     os.path.join('saved_results_' + name, 'ions_z.csv')]
    print('>> Analyzing ion movement through the filter...')
    if not all(os.path.exists(result) for result in saved_results):
        for frame in pd.read_fwf(os.path.join('temp_pdb_' + name, 'potassium_ions.pdb'), chunksize=num_ions + 1, header=None, colspecs=PDB_column_width_format, skiprows=1):
            # Delete last line (containing just 'END'), reset index to 0, and name columns
            frame.drop(frame.index[[-1]], inplace=True)
            frame.reset_index(drop=True, inplace=True)
            frame.columns = ['Atom', 'Index', 'Type', 'Residue', 'Chain', 'ResNum', 'X', 'Y', 'Z', 'Occupancy', 'TempFactor', 'Segment']
            # Drop every row except those containing ion of interest
            frame = frame[(frame['Residue'] == 'K+') | (frame['Residue'] == 'POT')]
            # Record z values of each ion of interest to dataframe containing the results so far
            ions_x[frame_count] = frame['X'].values
            ions_y[frame_count] = frame['Y'].values
            ions_z[frame_count] = frame['Z'].values
            frame_count += 1
            if frame_count % 50 == 0:
                print('\r   PROGRESS:     ', frame_count, end=' frames', flush=True)
        print('\r   PROGRESS:     ', frame_count, end=' frames (DONE)\n', flush=True)
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
title_size = arg.size
label_size = arg.size
fig1, axes1 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filterCA_diameters', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig2, axes2 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filterO_diameters', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig3, axes3 = plt.subplots(4, 1, sharex='col', figsize=(19, 14), num='poreCA_diameters', gridspec_kw={'height_ratios': [1, 1, 1, 2]})
fig4, axes4 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filterCA_symmetry', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig5, axes5 = plt.subplots(6, 1, sharex='col', figsize=(19, 20), num='filterO_symmetry', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 2]})  # rows, columns
fig6, axes6 = plt.subplots(4, 1, sharex='col', figsize=(19, 14), num='poreCA_symmetry', gridspec_kw={'height_ratios': [1, 1, 1, 2]})  # rows, columns
fig7, axes7 = plt.subplots(3, 1, figsize=(19, 17), num='phi', gridspec_kw={'height_ratios': [5, 3, 3]})  # rows, columns
fig8, axes8 = plt.subplots(3, 1, figsize=(19, 17), num='psi', gridspec_kw={'height_ratios': [5, 3, 3]})  # rows, columns
fig9, axes9 = plt.subplots(3, 1, figsize=(19, 13.9), num='chi1', gridspec_kw={'height_ratios': [3, 3, 3]})  # rows, columns
fig10, axes10 = plt.subplots(3, 1, figsize=(19, 10.8), num='chi2', gridspec_kw={'height_ratios': [1.3, 3, 3]})  # rows, columns
fig11, axes11 = plt.subplots(7, 1, sharex='col', figsize=(19, 23.3), num='water_filter_pore', gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 1, 2]})
axes_list = [axes1, axes2, axes3, axes4, axes5, axes6, axes7, axes8, axes9, axes10, axes11]
axes_list_heatmap = [axes7, axes8, axes9, axes10]
axes_list_noheatmap = [axes1, axes2, axes3, axes4, axes5, axes6, axes11]

################################################################################################

# Change index order so filter residues would be on top while pore residues would be on the bottom
index_order = [8, 6, 4, 2, 0, 9, 7, 5, 3, 1, 10, 11, 12]  # first four are filter CA, next four are filter O, last three are pore CA
area = area.reindex(index=[area.index[i] for i in index_order])
symmetry = symmetry.reindex(index=[symmetry.index[i] for i in index_order])
ACs = ACs.reindex(index=[ACs.index[i] for i in index_order])
BDs = BDs.reindex(index=[BDs.index[i] for i in index_order])
# changes_in_area = area.sub(area.iloc[:, 0], axis=0)

# Plot areas and symmetry of filter and pore over time
print('>> Plotting diameters and symmetry of filter and pore over time...')
for row, label in zip(range(5),  reversed(filter_labels)):
    # Areas
    # sns.lineplot(data=area.iloc[row], ax=axes1[row])
    # sns.lineplot(x=[time[0], time[-1]], y=[area.iloc[row][0]] * 2, color='navy', alpha=0.7, linestyle='--', ax=axes1[row])  # horizontal line to depict initial area
    axes1[row].stackplot(time, ACs.iloc[row].values, BDs.iloc[row].values, labels=['A to C', 'B to D'], colors=['#F3D3BD', '#5E5E5E'])
    # sns.lineplot(x=time, y=ACs.iloc[row].values, label='A to C', color='#F3D3BD', ax=axes1[row])
    # sns.lineplot(x=time, y=BDs.iloc[row].values, label='B to D', color='#5E5E5E', ax=axes1[row])
    axes1[0].legend(title='Subunit dist.', loc='center left', fancybox=True, bbox_to_anchor=(1, 0.5))
    axes1[row].set_ylabel(label, fontsize=label_size)
    # Symmetry Ratios
    sns.lineplot(data=symmetry.iloc[row], ax=axes4[row], color='black')
    sns.lineplot(x=[time[0], time[-1]], y=[1] * 2, color='gray', alpha=0.7, linestyle='--', ax=axes4[row])  # horizontal line to depict perfect symmetry (1)
    axes4[row].set_ylabel(label, fontsize=label_size)
    turn_y_axis_symmetric(axes4[row])
axes1[0].set_title('$C\\alpha$-$C\\alpha$ Distances Between Filter Residues - ' + name, y=1.04, fontsize=title_size)
axes4[0].set_title('$C\\alpha$-$C\\alpha$ Symmetry Between Filter Residues - ' + name, y=1.04, fontsize=title_size)

for row, label in zip(range(5, 10),  reversed(filter_labels)):
    # Areas
    # sns.lineplot(data=area.iloc[row], ax=axes2[row - 5], color='red')
    # sns.lineplot(x=[time[0], time[-1]], y=[area.iloc[row][0]] * 2, color='firebrick', alpha=0.7, linestyle='--', ax=axes2[row - 5])  # create horizontal line to depict initial area
    axes2[row - 5].stackplot(time, ACs.iloc[row].values, BDs.iloc[row].values, labels=['A to C', 'B to D'], colors=['#A2A7A5', '#6D696A'])
    # sns.lineplot(x=time, y=ACs.iloc[row].values, label='A to C', color='#A2A7A5', ax=axes2[row - 5])
    # sns.lineplot(x=time, y=BDs.iloc[row].values, label='B to D', color='#6D696A', ax=axes2[row - 5])
    axes2[0].legend(title='Subunit dist.', loc='center left', fancybox=True, bbox_to_anchor=(1, 0.5))
    axes2[row - 5].set_ylabel(label, fontsize=label_size)
    # Symmetry Ratios
    sns.lineplot(data=symmetry.iloc[row], ax=axes5[row - 5], color='black')
    sns.lineplot(x=[time[0], time[-1]], y=[1] * 2, color='gray', alpha=0.7, linestyle='--', ax=axes5[row - 5])  # horizontal line to depict perfect symmetry (1)
    axes5[row - 5].set_ylabel(label, fontsize=label_size)
    turn_y_axis_symmetric(axes5[row - 5])
axes2[0].set_title('Backbone Carbonyl Oxygen Distances Between Filter Residues - ' + name, y=1.04, fontsize=title_size)
axes5[0].set_title('Backbone Carbonyl Oxygen Symmetry Between Filter Residues - ' + name, y=1.04, fontsize=title_size)

for row, label in zip(range(10, 13), pore_labels):
    # Areas
    # sns.lineplot(data=area.iloc[row], ax=axes3[row - 10])
    # sns.lineplot(x=[time[0], time[-1]], y=[area.iloc[row][0]] * 2, color='navy', alpha=0.7, linestyle='--', ax=axes3[row - 10])  # create horizontal line to depict initial area
    axes3[row - 10].stackplot(time, ACs.iloc[row].values, BDs.iloc[row].values, labels=['A to C', 'B to D'], colors=['#B1E5F2', '#272635'])
    # sns.lineplot(x=time, y=ACs.iloc[row].values, label='A to C', color='#B1E5F2', ax=axes3[row - 10])
    # sns.lineplot(x=time, y=BDs.iloc[row].values, label='B to D', color='#272635', ax=axes3[row - 10])
    axes3[0].legend(title='Subunit dist.', loc='center left', fancybox=True, bbox_to_anchor=(1, 0.5))
    axes3[row - 10].set_ylabel(label, fontsize=label_size)
    # Symmetry Ratios
    sns.lineplot(data=symmetry.iloc[row], ax=axes6[row - 10], color='black')
    sns.lineplot(x=[time[0], time[-1]], y=[1] * 2, color='seagreen', alpha=0.7, linestyle='--', ax=axes6[row - 10])  # horizontal line to depict perfect symmetry (1)
    axes6[row - 10].set_ylabel(label, fontsize=label_size)
    turn_y_axis_symmetric(axes6[row - 10])
axes3[0].set_title('$C\\alpha$-$C\\alpha$ Distances Between Drug Binding Pore Residues - ' + name, y=1.04, fontsize=title_size)
axes6[0].set_title('$C\\alpha$-$C\\alpha$ Symmetry Between Drug Binding Pore Residues - ' + name, y=1.04, fontsize=title_size)

################################################################################################

# Plot phi/psi, chi1/chi2
print('>> Plotting phi, psi, chi1, chi2 angles...')
# Selectivity filter residues
angles_data_column_width_format = [(0, 5), (5, 12), (12, 17), (17, 23), (23, 29), (29, 35), (35, 40)]

phi = {}
psi = {}
chi1 = {}
chi2 = {}
for segname in ['PROA', 'PROB', 'PROC', 'PROD']:
    # Obtain angles data
    angles_data = pd.read_fwf(os.path.join('temp_pdb_' + name, segname + '_filter_angles.dat'), header=None, colspecs=angles_data_column_width_format)
    angles_data.columns = ['ResNum', 'Residue', 'Frame', 'Phi', 'Psi', 'Chi1', 'Chi2']
    # Use data at initial frame to check and make sure obtained residues are correct and match with pre-defined labels before proceeding
    filter_resnum = angles_data.loc[angles_data['Frame'] == 0, 'ResNum'].values
    for label, ref_residue in zip(filter_labels, angles_data.loc[angles_data['Frame'] == 0, 'Residue'].values):
        if label[0] != aa_names.get(ref_residue):
            print('ERROR: Data for filter angles yield residue', aa_names.get(ref_residue), 'while', label[0], 'was expected in the provided labels. Check vmd_cmd_file section in this script.')
            exit(1)
    # Acquire data and place them in a dictionary
    for angle_type, dataset in zip(['Phi', 'Psi', 'Chi1', 'Chi2'], [phi, psi, chi1, chi2]):
        for resnum, row, label in zip(filter_resnum, range(5), filter_labels):
            # If dataset is filled with NaN, skip plotting
            if not np.isnan(angles_data.loc[angles_data['ResNum'] == resnum, angle_type].values).all():
                current_residue = label + ':' + segname[-1]
                dataset[current_residue] = angles_data.loc[angles_data['ResNum'] == resnum, angle_type].values

# Plot angles
for angle_type, ax, fig, dataset in zip(['Phi', 'Psi', 'Chi1', 'Chi2'], [axes7, axes8, axes9, axes10], [fig7, fig8, fig9, fig10], [phi, psi, chi1, chi2]):
    dataset = pd.DataFrame(dataset, index=time)
    # Reindex so that the labels go 1.A, 1.B, 1.C, 2.A, 2.B... instead of 1.A, 2.A, 3.A, 1.B, 2.B,...
    dataset = dataset.reindex(sorted(dataset.columns, key=lambda x: x[1:-2], reverse=True), axis=1)  # reversed order to match visual orientation when viewed by Chimera or VMD
    my_map = sns.heatmap(dataset.T, vmin=-180, vmax=180,
                         xticklabels=False, yticklabels=1, ax=ax[0],
                         cmap=sns.color_palette('twilight_shifted', as_cmap=True),
                         cbar=False)
    # Add lines to separate 1.A, 1.B, 1.C into 1 group, 2.A, 2.B, 2.C into another group...
    ax[0].hlines(np.arange(4, len(dataset.columns.values), 4), *ax[0].get_xlim(), colors='white')
    # Add black frame around heatmap
    for _, spine in my_map.spines.items():
        spine.set_visible(True)
    print('\r   PROGRESS:     ', angle_type, end=' for SELECTIVITY FILTER', flush=True)
print('\r   PROGRESS:     ', angle_type, end=' for SELECTIVITY FILTER (DONE)\n', flush=True)

# Pore residues
phi = {}
psi = {}
chi1 = {}
chi2 = {}
for segname in ['PROA', 'PROB', 'PROC', 'PROD']:
    # Obtain angles data and concatenate angles data
    Y_angles_data = pd.read_fwf(os.path.join('temp_pdb_' + name, segname + '_poreY_angles.dat'), header=None, colspecs=angles_data_column_width_format)
    FS_angles_data = pd.read_fwf(os.path.join('temp_pdb_' + name, segname + '_poreFS_angles.dat'), header=None, colspecs=angles_data_column_width_format)
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
    # Acquire data and place them in a dictionary
    for angle_type, dataset in zip(['Phi', 'Psi', 'Chi1', 'Chi2'], [phi, psi, chi1, chi2]):
        for resnum, row, label in zip(pore_resnum, range(3), pore_labels):
            # If dataset is filled with NaN, skip plotting
            if not np.isnan(angles_data.loc[angles_data['ResNum'] == resnum, angle_type].values).all():
                current_residue = label + ':' + segname[-1]
                dataset[current_residue] = angles_data.loc[angles_data['ResNum'] == resnum, angle_type].values

# Plot angles
for angle_type, ax, dataset in zip(['Phi', 'Psi', 'Chi1', 'Chi2'], [axes7, axes8, axes9, axes10], [phi, psi, chi1, chi2]):
    dataset = pd.DataFrame(dataset, index=time)
    # Reindex so that the labels go 1.A, 1.B, 1.C, 2.A, 2.B... instead of 1.A, 2.A, 3.A, 1.B, 2.B,...
    dataset = dataset.reindex(sorted(dataset.columns, key=lambda x: x[1:-2]), axis=1)
    my_map = sns.heatmap(dataset.T, vmin=-180, vmax=180,
                         xticklabels=False, yticklabels=1, ax=ax[1],
                         cmap=sns.color_palette('twilight_shifted', as_cmap=True),
                         cbar=False)
    # Add lines to separate 1.A, 1.B, 1.C into 1 group, 2.A, 2.B, 2.C into another group...
    ax[1].hlines(np.arange(4, len(dataset.columns.values), 4), *ax[1].get_xlim(), colors='white')
    # Add black frame around heatmap
    for _, spine in my_map.spines.items():
        spine.set_visible(True)
    print('\r   PROGRESS:     ', angle_type, end=' for PORE', flush=True)
print('\r   PROGRESS:     ', angle_type, end=' for PORE (DONE)\n', flush=True)

# Set titles
axes7[0].set_title('Phi Angles of Filter & Drug Binding Pore Residues - ' + name, y=1.04, fontsize=title_size)
axes8[0].set_title('Psi Angles of Selectivity Filter & Drug Binding Pore Residues - ' + name, y=1.04, fontsize=title_size)
axes9[0].set_title('Chi1 Angles of Selectivity Filter & Drug Binding Pore Residues - ' + name, y=1.04, fontsize=title_size)
axes10[0].set_title('Chi2 Angles of Selectivity Filter & Drug Binding Pore Residues - ' + name, y=1.04, fontsize=title_size)

################################################################################################

# Plot water in pore and in filter
print('>> Plotting numbers of water molecules in pore and in filter...')
water_data = pd.read_csv(os.path.join('temp_pdb_' + name, 'water.dat'), header=None)
water_data.columns = ['Frame', 'Inside Filter', 'Behind Filter$_{A}$', 'Behind Filter$_{B}$', 'Behind Filter$_{C}$', 'Behind Filter$_{D}$', 'Inside Pore']

for row, label in zip(range(6), water_data.columns[1:]):
    sns.lineplot(x=time, y=water_data[label].values, color='darkviolet', ax=axes11[row])
    axes11[row].set_ylabel(label, fontsize=label_size)

axes11[0].set_title('Water Count - ' + name, y=1.04, fontsize=title_size)

################################################################################################

# Plot drug movement in the pore if applicable, else plot ion movement
palette = cycle(sns.color_palette(n_colors=8))
if arg.drug:
    # Plot drug movement in the pore if specified in the input
    print('>> Plotting drug movement in the pore...')

    for (mol, max_z), (mol, com_z), (mol, min_z) in zip(max_z_data.iterrows(), com_z_data.iterrows(), min_z_data.iterrows()):
        # If dataset is filled with NaN, skip
        if not np.isnan(com_z.values).all():
            for ax in axes_list:
                ax[-1].plot(time, com_z.values)
                ax[-1].fill_between(time, min_z.values, max_z.values, alpha=0.2)
        print('\r   PROGRESS:     ', mol, end=' molecules', flush=True)
    print('\r   PROGRESS:     ', mol, end=' molecules (DONE)\n', flush=True)

    # Plot reference lines
    for ax in axes_list:
        sns.lineplot(x=[time[0], time[-1]], y=[S4_mean] * 2, color='blue', alpha=0.7, linestyle='--', ax=ax[-1])
        sns.lineplot(x=time, y=pore_end_mean, color='red', alpha=0.7, linestyle='--', ax=ax[-1])
        # Annotate reference positions, position of annotation is determined automatically
        # ax[-1].text(arg.time_total // 20, S4_mean, r'${K^+} S4$', size=label_size - 4, color='blue', ha='center', va='center', path_effects=[pe.withStroke(linewidth=7, foreground='white')])
        ax[-1].text(arg.time_total // 20 * 19, S4_mean, r'${K^+} S4$', size=label_size - 4, color='blue', ha='center', va='center', path_effects=[pe.withStroke(linewidth=7, foreground='white')])
        # ax[-1].text(arg.time_total // 20, pore_end_mean, 'Pore End', size=label_size - 4, color='red', ha='center', va='center', path_effects=[pe.withStroke(linewidth=7, foreground='white')])
        ax[-1].text(arg.time_total // 20 * 19, pore_end_mean, 'Pore End', size=label_size - 4, color='red', ha='center', va='center', path_effects=[pe.withStroke(linewidth=7, foreground='white')])

        # Annotate Center of Mass in z of drug binding pore residues in the first frame
        ax[-1].text(arg.time_total // 20 * 19, drug_binding_res_data['Y652'][int(len(drug_binding_res_data['Y652']) // 20 * 19)], 'Y652', size=label_size - 4, color='darkslategray', ha='center', va='center', path_effects=[pe.withStroke(linewidth=7, foreground='white')])
        sns.lineplot(x=time, y=drug_binding_res_data['Y652'], color='darkslategray', alpha=0.7, linestyle='--', ax=ax[-1])

        ax[-1].text(arg.time_total // 20 * 19, drug_binding_res_data['F656'][int(len(drug_binding_res_data['F656']) // 20 * 19)], 'F656', size=label_size - 4, color='indianred', ha='center', va='center', path_effects=[pe.withStroke(linewidth=7, foreground='white')])
        sns.lineplot(x=time, y=drug_binding_res_data['F656'], color='indianred', alpha=0.7, linestyle='--', ax=ax[-1])

        ax[-1].text(arg.time_total // 20 * 19, drug_binding_res_data['S660'][int(len(drug_binding_res_data['S660']) // 20 * 19)], 'S660', size=label_size - 4, color='seagreen', ha='center', va='center', path_effects=[pe.withStroke(linewidth=7, foreground='white')])
        sns.lineplot(x=time, y=drug_binding_res_data['S660'], color='seagreen', alpha=0.7, linestyle='--', ax=ax[-1])

else:
    # Plot ion movement over time
    print('>> Plotting ion movement over time...')
    for (ion, x), (ion, y), (ion, z) in zip(ions_x.iterrows(), ions_y.iterrows(), ions_z.iterrows()):
        z_in_range = []
        time_in_range = []
        temp_z = []
        temp_t = []

        for x, y, z, t in zip(x.values, y.values, z.values, time):
            if (sqrt(x**2 + y**2) <= 10 and S4_mean <= z <= S0_mean) or (sqrt(x**2 + y**2) <= 15 and (S4_mean - 10 <= z < S4_mean or S0_mean < z <= S0_mean + 15)):
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
            color = next(palette)
            for t, z in zip(time_in_range, z_in_range):
                for ax in axes_list:
                    # ax[-1].axvspan(550, 720, color='yellow', alpha=0.2)
                    sns.lineplot(x=t, y=z, ax=ax[-1])

    # Plot reference lines
    for z, label in zip([S0_mean, S1_mean, S2_mean, S3_mean, S4_mean], ['S0', 'S1', 'S2', 'S3', 'S4']):
        for ax in axes_list:
            sns.lineplot(x=[time[0], time[-1]], y=[z] * 2, color='blue', alpha=0.7, linestyle='--', ax=ax[-1])
            # Annotate S positions, position of annotation is determined automatically
            ax[-1].text(arg.time_total // 20, z, label, size=label_size - 4, color='blue', ha='center', va='center', path_effects=[pe.withStroke(linewidth=7, foreground='white')])
            ax[-1].text(arg.time_total // 20 * 19, z, label, size=label_size - 4, color='blue', ha='center', va='center', path_effects=[pe.withStroke(linewidth=7, foreground='white')])

# Set labels
for ax in axes_list:
    if arg.drug:
        ax[-1].set_ylabel('Drug $z-pos.$ in Pore (Å)', fontsize=label_size)
    else:
        ax[-1].set_ylabel('$K^+ z-pos.$ (Å)', fontsize=label_size)
    ax[-1].set_xlabel(arg.x_label, fontsize=label_size)

for ax in axes_list_heatmap:
    ax[-1].xaxis.set_minor_locator(AutoMinorLocator())
    ax[-1].set(xlim=(0, time[-1]))
    # Color bar
    axins = inset_axes(ax[0],
                       width=0.3,
                       height=5,
                       loc='upper right',
                       bbox_to_anchor=(0.04, 0, 1, 1),
                       # (x0, y0, width, height) where (x0,y0) are the lower left corner coordinates of the bounding box
                       bbox_transform=ax[0].transAxes,
                       borderpad=0)
    cb1 = matplotlib.colorbar.ColorbarBase(axins, cmap=sns.color_palette('twilight_shifted', as_cmap=True),
                                           norm=matplotlib.colors.Normalize(vmin=-180, vmax=180),
                                           orientation='vertical')
    cb1.set_label('Angles (°)')
    # Generate X tick labels for heatmap (heatmap order is categorical, not numeric, so default xticks will be 103, 134, 162,... instead of 100, 200, 300,...)
    major_xticks = list(ax[-1].get_xticks(minor=False))[:-1]  # Skip last one because it is not plotted for some reason
    minor_xticks = list(ax[-1].get_xticks(minor=True))
    major_xtick_locations = []
    major_xtick_labels = []
    minor_xtick_locations = []
    # When I wrote this code, only God and I knew how it worked. Now, only God knows how it works.
    for t in major_xticks:
        correction_factor = 1 / len(time) * closest(time, t)
        time_index = round(len(time) * t / arg.time_total + correction_factor, 2)
        major_xtick_locations.append(time_index)
        major_xtick_labels.append(int(t))
    for t in minor_xticks:
        correction_factor = 1 / len(time) * closest(time, t)
        time_index = round(len(time) * t / arg.time_total + correction_factor, 2)
        minor_xtick_locations.append(time_index)
    # Set xticks for selectivity filter residues
    ax[0].set_xticks(major_xtick_locations)
    ax[0].set_xticklabels(major_xtick_labels)
    ax[0].set_xticks(minor_xtick_locations, minor=True)
    # Set xticks for pore residues
    ax[1].set_xticks(major_xtick_locations)
    ax[1].set_xticklabels(major_xtick_labels)
    ax[1].set_xticks(minor_xtick_locations, minor=True)

for ax in axes_list_noheatmap:
    for i in ax:
        i.xaxis.set_minor_locator(AutoMinorLocator())

################################################################################################

# Set limits
# Filter CA radii
# axes1[0].set(ylim=(0, 40))
# axes1[1].set(ylim=(0, 30))
# axes1[2].set(ylim=(0, 27))
# axes1[3].set(ylim=(0, 27))
# axes1[4].set(ylim=(0, 27))

# Filter O radii
# axes2[0].set(ylim=(0, 40))
# axes2[1].set(ylim=(0, 35))
# axes2[2].set(ylim=(0, 30))
# axes2[3].set(ylim=(0, 30))
# axes2[4].set(ylim=(0, 22))

# Filter CA symmetry
# axes4[0].set(ylim=(0.6, 3))
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
# axes3[0].set(ylim=(0, 50))
# axes3[1].set(ylim=(0, 50))
# axes3[2].set(ylim=(0, 50))

# Ions
if arg.drug:
    for ax in axes_list:
        ax[-1].set(ylim=(pore_end_mean - 7, S4_mean + 3))
else:
    for ax in axes_list:
        ax[-1].set(ylim=(S4_mean - 7, S0_mean + 3))

# Water in filter and pore
# axes11[0].set(ylim=(-0.4, 6.4))
# axes11[1].set(ylim=(-0.4, 10.4))
# axes11[2].set(ylim=(-0.4, 10.4))
# axes11[3].set(ylim=(-0.4, 10.4))
# axes11[4].set(ylim=(-0.4, 10.4))
# axes11[5].set(ylim=(-0.4, 150.4))

################################################################################################

# Save figures

print('\nPlots are saved as:')
for figure in ['filterCA_diameters', 'filterO_diameters', 'poreCA_diameters',
               'filterCA_symmetry', 'filterO_symmetry', 'poreCA_symmetry',
               'phi', 'psi', 'chi1', 'chi2',
               'water_filter_pore']:

    plt.figure(figure)

    plt.savefig(os.path.join('figures_' + name, figure + '.png'), dpi=300)
    print(os.path.join('figures_' + name, figure + '.png'))
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
            plt.savefig(os.path.join('figures_' + name, figure + '_part' + str(portion + 1) + '.png'), dpi=300)
            print(os.path.join('figures_' + name, figure + '_part' + str(portion + 1) + '.png'))

print("\n _._     _,-'\"\"`-._\n,-.`._,'(       |\\`-/|\n    `-.-' \\ )-`( , o o)\n          `-    \\`_`\"'-")
