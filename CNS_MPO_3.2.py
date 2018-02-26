#!/usr/bin/env python
# ----------------------------------------------------------------------#
# 			CNS MPO Score and Solubility forecaster index (3.2)			#
# --------------------------- Compatibility ----------------------------#
# 					Compatible with Python 2.7 & 3.5					#
# ----------------------Prerequiste installations-----------------------#
# 																		#
# 		- Chemaxon Marvin suite WITH license							#
# 				cxcalc module required to work 	(used for pKa pred)		#
# 					(make sure to add cxcalc to $PATH)					#
# by default : C:\Program Files (x86)\ChemAxon\MarvinBeans\bin			#
# How to use :															#
# 		- input the sdf file name when requested (.sdf included)		#
# 		- the script output an sdf with #name_out.sdf					#
# 			This sdf contains the fields with : 						#
# 						- CNS MPO Score									#
# 						- Solubility forecaster index (SFI)[optional]	#
# 						- bpKa,logD(7.4),logP,MW,HBD,#Ar,TPSA			#
# 						- All individual components of MPO Score		#
# ----------------------------------------------------------------------#
# (C) Dr. De Cesco Stephane - v3.2 - 22/03/2017							#
# ----------------------------------------------------------------------#
# -------------------------- MODULES IMPORT ----------------------------#
from __future__ import print_function

import argparse
import os
import os.path
import subprocess
import sys

# To remove danger of using input in python2 environments.
try:
    input = raw_input
except NameError:
    pass


# ----------------------------- Desc object -----------------------------#

def monotonic_score(value, lower, upper):
    # =======================================================================================
    #         upper
    # 1|---------
    #  |		  \
    #  |           \@value
    #  |			\
    # 0|_____________\______________________
    # 	             lower
    # Function to return a score between 0 and 1 depending of the value of the parameters.
    # 		 | 1 if value < upper
    # Score ={ f(upper and lower) if upper < value < lower
    # 		 | 0 if value > lower
    # =======================================================================================
    try:
        value = float(value)
    except ValueError as message:
        print(message)
        print(value)
    upper = float(upper)
    lower = float(lower)
    if value <= lower:
        score = 1
    elif value > upper:
        score = 0
    else:
        score = 1 - ((value - lower) * (1 / (upper - lower)))
    return score


def hump_score(value, low1, up1, up2, low2):
    # =======================================================================================
    #         	up1		 up2
    # 1|		  --------
    #  |		 / 		  \
    #  |        /      	   \
    #  |	   /			\
    # 0|______/______________\_______________
    # 	     low1			  low2
    # Function to return a score between 0 and 1 depending of the value of the parameters.
    # 		 | 0 if value < low1
    # 		 | f (low1 and up1) if low1 < value < up1
    # Score ={ 1 if up1 < value < up2
    # 		 | f (up2 and low2) if up2 < value < low2
    # 		 | 0 if value > lower
    # =======================================================================================

    value, low1, up1, up2, low2 = float(value), float(low1), float(up1), float(up2), float(low2)
    score = 0
    if value <= low1:
        score = 0
    elif up1 < value <= up2:
        score = 1
    elif value > low2:
        score = 0
    elif low1 < value <= up1:
        score = ((value - low1) * (1 / (up1 - low1)))
    elif up2 < value <= low2:
        score = 1 - ((value - up2) * (1 / (low2 - up2)))
    return score


class Desc:
    def __init__(self, mol):
        self.mol = mol  # stores the mol coordinate into .mol
        self.Prop = {}  # Create an empty dictionnary to store all the individual Properties
        self.get_param()
        if args.sfi:
            self.Prop['SFI'] = float(self.Prop['logD']) + float(self.Prop['ArRings'])  # Calc SFI
        self.calc_mpo_score()  # Calculate the MPO score
        # self.calc_mpo_area()

    def get_param(self):
        record_data = False
        coordinates = True
        new_mol = ''
        transcription = {'LOGD': 'logD', 'LOGP': 'logP', 'PSA': 'TPSA', 'MASS': 'MW', 'AROMATIC_RINGCOUNT': 'ArRings',
                         'pkacalculator': 'bpKa', 'DONOR_COUNT': 'HBD', 'CHARGE_DISTRIBUTION': 'charge_dist(7.4)'}

        for lines in self.mol.splitlines():
            if lines.startswith('>  <'):
                c_lines = lines.lstrip('>  <')
                chemprop = c_lines.rstrip('>')
                if chemprop in transcription.keys() or chemprop in transcription.values():
                    coordinates = False
            if coordinates:
                new_mol += lines + '\n'
            if record_data:
                value = str(lines)
                if chemprop == 'LOGD':
                    value = lines.split('\t')[1].strip()
                if chemprop == 'pkacalculator':
                    value = lines.split('\t')[0]
                    if value.strip() == "":
                        value = 0.0
                if chemprop == 'CHARGE_DISTRIBUTION':
                    value = lines.split('\t')[1].lstrip()
                if chemprop == 'Name':
                    chemprop = chemprop.lower()
                if chemprop in transcription.keys():
                    self.Prop[transcription[chemprop]] = value
                elif chemprop in transcription.values():
                    self.Prop[chemprop] = value
                record_data = False
            if lines.startswith('>  <'):
                lines = lines.lstrip('>  <')
                chemprop = lines.rstrip('>')
                record_data = True
        self.mol = new_mol

    def calc_mpo_score(self):  # call the monotonic or hump score function for each term with the boundaries and sum them at the
        #  end to populate the CNS MPO Score
        try:
            self.Prop['bpKaScore'] = float(monotonic_score(self.Prop['bpKa'], 8,
                                                           10))  # Todo : See with hump score for pKa (https://doi.org/10.1021/acs.jmedchem.6b01469)
            self.Prop['logPScore'] = float(monotonic_score(self.Prop['logP'], 3, 5))
            self.Prop['logDScore'] = float(monotonic_score(self.Prop['logD'], 2, 4))
            self.Prop['MWScore'] = float(monotonic_score(self.Prop['MW'], 360, 500))
            self.Prop['HBDScore'] = float(monotonic_score(self.Prop['HBD'], 0.5, 3.5))
            self.Prop['TPSAScore'] = float(hump_score(self.Prop['TPSA'], 20, 40, 90, 120))
            self.Prop['MPOScore'] = self.Prop['bpKaScore'] + self.Prop['logPScore'] + self.Prop['logDScore'] \
                                    + self.Prop['MWScore'] + self.Prop['HBDScore'] + self.Prop['TPSAScore']
            self.Prop['MPOScore_v2'] = self.Prop['bpKaScore'] + self.Prop['logPScore'] + self.Prop['MWScore'] + (
                2 * self.Prop['HBDScore']) + self.Prop['TPSAScore']
        except KeyError as missing:

            print(
                "The following parameter field [", missing,
                "] in the sdf file is missing or is not spelled as required "
                "(case sensitive):\n- \"bpKa\" (basic pKa)\n- \"logD\" ("
                "logD at pH=7.4)\n- \"logP\"\n- \"HBD\" (Hydrogen bond "
                "donor count)\n- \"TPSA\" (total polar surface area)\n- \"MW\" (molecular weight)\nMolecule skipped")
            self.Prop['bpKaScore'] = "Error"
            self.Prop['logPScore'] = "Error"
            self.Prop['logDScore'] = "Error"
            self.Prop['MWScore'] = "Error"
            self.Prop['HBDScore'] = "Error"
            self.Prop['TPSAScore'] = "Error"
            self.Prop['MPOScore'] = "Error"
            self.Prop['MPOScore_v2'] = "Error"


    def print_details(self, level):  # For visualisation in the terminal only
        if level == 0:
            print("======================")
            # noinspection PyBroadException
            try:
                print(self.Prop['name'])
            except:
                print("No name inputed")
            print('CNS MPO Score = ' + str(self.Prop['MPOScore']))
            # print('MPO percentage area = ', self.Prop['MPO_area'])
            print("======================")
        if level == 1:
            print("======================")
            # noinspection PyBroadException
            try:
                print(self.Prop['name'])
            except:
                print("No name inputed")
            for key, value in self.Prop.items():
                print(key + ": " + str(value))
            print("======================")

    def sdf_writer(self):
        mol = self.mol.replace("$$$$", '')
        mol = mol.rstrip('\n') + '\n'
        for key, value in self.Prop.items():
            mol = mol + "\n>  <" + key + ">\n"
            mol = mol + str(value) + "\n"
        mol += "\n$$$$\n"
        return mol

    def calc_mpo_area(self):
        values = [self.Prop['bpKaScore'], self.Prop['logPScore'], self.Prop['logDScore'], self.Prop['MWScore'],
                  self.Prop['HBDScore'], self.Prop['TPSAScore']]
        values.sort(reverse=True)
        areas = []
        for i in range(len(values)):
            if i == 5:
                area = ((values[i] * values[0]) / 2) * 0.866  # 0.866 = sin(60deg)
            else:
                area = ((values[i] * values[i + 1]) / 2) * 0.866  # 0.866 = sin(60deg)
            areas.append(area)
        total_area = (0.866 / 2) * 6
        fraction_area = sum(areas) / total_area
        self.Prop['MPO_area'] = round(fraction_area * 100, ndigits=2)


def calc_sdf_parser(sdf_param):
    molecule_list = []
    molecule = ''
    lines = sdf_param.splitlines(keepends=True)
    for mol_line in lines:
        molecule += mol_line
        if mol_line == "$$$$\n" or mol_line == '$$$$\r\n':
            molecule_list.append(molecule)
            molecule = ''
    return molecule_list


def sdf_parser(file):
    with open(file, 'r') as sdf_file:
        molecule_list = []
        molecule = ''
        lines = sdf_file.readlines()
        for mol_line in lines:
            molecule += mol_line
            if mol_line == "$$$$\n" or mol_line == '$$$$\r\n':
                molecule_list.append(molecule)
                molecule = ''
    return molecule_list


def calc_param(file):
    """
    :param file:
    :return returns a SDF file with all the calculated parameters:
    Uses cxcalc from chemaxon as chemical property calculator. Please refer to their website for installation and licensing
    """

    command_1 = 'cxcalc -i name -S pkacalculator -t basic -b 1 logD -H 7.4 logP donorcount polarsurfacearea mass ' \
                'aromaticringcount rotatablebondcount '
    if str(file).startswith('/'):
        command_1 += str(file)
    else:
        command_1 += os.getcwd() + '/' + str(file)

    sdf_file = subprocess.check_output([command_1], shell=True)

    sdf_file = sdf_file.decode()
    return sdf_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_file', help='Name of the input file [MENDATORY]', metavar='')
    parser.add_argument('-S', '--sfi',
                        help='Use this flag if you want to report SFI index (solubility Forecast index) the field '
                             '"ArRings" (count of '
                             'aromatics rings) becomes required in the sdf fields if the -C option is used (by '
                             'default option is off) ', action='store_true', default=False)
    parser.add_argument('-C', '--no_calc',
                        help='Use this flag to disable automatic calculation of chemical properties (need to be '
                             'provided in the sdf file then (beware the sdf fields need to have the same name as '
                             'following: bpKa (basic pKa), logD (logD at pH=7.4), logP, HBD (Hydrogen bond donor '
                             'count), TPSA (total polar surface area)'
                             ', MW (molecular weight)', action='store_true', default=False)
    parser.add_argument('-f', '--output_folder', type=str, help='create a folder to output the file', metavar='')
    parser.add_argument('-v', '--verbose', help='more details output in the terminal', action='store_true',
                        default=False)
    args = parser.parse_args()
    if not args.in_file:
        parser.print_help()
        sys.exit()
    while True:
        if args.in_file:
            if os.path.exists(args.in_file):
                break
            else:
                print('ERROR : file inputed as argument [-i] does not exist')
                sys.exit()
    if not args.no_calc:
        sdf_withparam = calc_param(args.in_file)
        mol_list = calc_sdf_parser(sdf_withparam)
    else:
        mol_list = sdf_parser(args.in_file)
    fpathname = args.in_file.split('/')  # parsing the file name to use it for the output
    fpath = '/'.join(fpathname[:-1]) + '/'
    if fpath == '/':
        fpath = ''
    fname = fpathname[-1].split('.')[0]
    mol_to_out = ''
    for m in mol_list:  # scanning through the molecules
        m = Desc(m)  # initializing the Desc object with all the parameters
        mol_to_out = mol_to_out + m.sdf_writer()
        if args.verbose:
            m.print_details(
                1)  # Print the details in the terminal. By default level = 0, enter Print_details(1) for more details
        else:
            m.print_details(0)
    if args.output_folder:
        if str(args.output_folder).startswith('/'):
            path_to_check = args.output_folder
        else:
            path_to_check = fpath + str(args.output_folder)
        if not os.path.exists(path_to_check):
            os.makedirs(path_to_check)
        with open(path_to_check + '/' + fname + "_out.sdf", 'w') as output:
            output.write(mol_to_out)
    else:
        with open(fpath + fname + "_out.sdf", 'w') as output:
            output.write(mol_to_out)
