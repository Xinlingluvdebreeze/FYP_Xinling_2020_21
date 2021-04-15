# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 12:57:17 2021

@author: Xinling
"""
import itertools as it
import os
import yaml      as ya
from   pathlib import Path
from itertools import permutations 

###############################################################################
#Non-Standard Imports
###############################################################################
try:
    import HiPAD.sort_align as san
except Exception as e:
    if Path.cwd() == Path(__file__).parent:
        import addpath
        import HiPAD.sort_align as san
    else:
        raise e
        
'''
scripts is a global variable where each key is protocol type and each 
value is the name of the file containing the appropriate opentrons script.
'''
###############################################################################
#Globals
###############################################################################
alphabets           = ['A','B','C','D','E','F','G','H']
position_in_96_well = ['A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12',
                       'B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12',
                       'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12',
                       'D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11','D12',
                       'E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12',
                       'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12',
                       'G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12',
                       'H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12',
                       ]

scripts             = {(2, 2): 'DNA_Assembly_2x2',
                       (2, 3): 'DNA_Assembly_2x3',
                       (1, 2, 2): 'DNA_Assembly_1x2x2',
                       (2, 2, 1): 'DNA_Assembly_2x2x1',
                       (1, 1, 3, 3): 'DNA_Assembly_1x1x3x3',
                       }

###############################################################################
#High-level Function
###############################################################################
'''
stage2ot is a high-level function which allows user to input stage and retrive 
script, source arrangement, output arrangement and volumes of reagents required

inputs2ot is a high-level function which sllows user to input all_inputs
(desired variant combinations) and retrieve the same info as stage2ot. 
Hence, yaml file is unnecessary
'''

def stage2ot(stage):
    
    all_inputs = get_inputs(stage)
    script, source_position, desired_variants_position, other_reagents_volume = inputs2ot(all_inputs)
    
    return script, source_position, desired_variants_position, other_reagents_volume

def inputs2ot(all_inputs):
    
    n_input_row     = check_n_positions(all_inputs)
    aligned_inputs  = align_inputs(all_inputs, n_input_row)
    source_position = map_to_source_well(aligned_inputs)
    protocol        = protocol_type(aligned_inputs)
    n_output        = total_output(protocol)
    script          = choose_script(protocol)
    reagents_fill_order       = reagents_fill(aligned_inputs, protocol, n_input_row)
    total_variants_list       = generate_variants_combinations(reagents_fill_order, aligned_inputs, protocol, n_input_row, n_output)
    output_position           = all_variants_position(n_output)
    desired_variants_position = desired_variant_output_well(all_inputs, total_variants_list, output_position)
    source_position, other_reagents_volume = optimal_reagents_volume(source_position, total_variants_list, n_input_row)
      
    return script, source_position, desired_variants_position, other_reagents_volume

###############################################################################
#Checks
###############################################################################
'''
check_stage function checks that recombination stage is not empty

check_n_positions function checks that the number of position is the same for each variant
'''
def check_stage(stage):
    if not stage:
        msg = f'Stage is empty and recombination step is not needed.'
        raise ValueError(msg)
        
def check_n_positions(all_inputs):

    n_positions = None
    for variant in all_inputs:
        if not n_positions:
            n_positions = len(variant)
        elif n_positions == len(variant):
            pass
        else:
            msg = f'The number of position is not the same for each variant.'
            raise NotImplementedError(msg)
    
    n_input_row = n_positions
    
    return n_input_row

###############################################################################
#Input Processing
############################################################################### 
'''
get_inputs function retrives all the inputs (desired variant combination)
from stage

align_inputs function removes repeating fragments in all_inputs and 
sorts them in ascending order, according to the length of fragment
'''
def get_inputs(stage): 
    all_inputs = [step['inputs'] for step in stage]        
    
    return all_inputs

def align_inputs(all_inputs, n_input_row):
    # Start off with first set of inputs in aligned_inputs
    step0_inputs   = all_inputs[0]
    aligned_inputs = [[step0_inputs[i]] for i in range(n_input_row)]
    indices        = list(range(len(aligned_inputs)))

    # Append fragments to each position in aligned_inputs if 
    # the fragment is not in it
    for position in all_inputs[1:]: 
        for i in indices:
            if position[i] not in aligned_inputs[i]:
                aligned_inputs[i].append(position[i])
            else: 
                pass
            
    # Finally, sort each lst in aligned_inputs so we no longer care about the 
    # initial order
    aligned_inputs = [sorted(lst) for lst in aligned_inputs]
    aligned_inputs = sorted(aligned_inputs, key=len)
    
    return aligned_inputs

###############################################################################
#Generate script to run in Opentron
###############################################################################
'''
protocol_type function generates the correct protocol

choose_script function tells user which script to run in Opentron
based on the protocol type
'''

def protocol_type(aligned_inputs):
    
    protocol_type = tuple([len(lst) for lst in aligned_inputs])
    
    return protocol_type

def choose_script(protocol_type):
    global scripts
    
    if protocol_type in scripts:
        script_info = scripts[protocol_type]
    else:
        raise NotImplementedError('No script can be mapped to yamldict')
    
    return script_info

###############################################################################
#Arrangement in source plate
###############################################################################
'''
2 96 well plates are used to contain sources
The first one is used to contain DNA fragments
The second plate is used to contain H2O, buffer and enzyme

map_to_source_well arranges the fragments in 96 well source plate

minimum_reagents_volume generates the minimum volume required for
each reagent used in DNA assembly
'''
def map_to_source_well(aligned_inputs):
    global alphabets
    
    n_fragments = 0
    for position in aligned_inputs:
        for fragment in position:
            n_fragments += 1
            
    source_position_list = []
    
    for i in range(n_fragments):
        source_position_list.append([])

    alphabet_count = 0
    number = 0
    for position in aligned_inputs:
        column_count = 1
        for fragment in position:
            source_position_list[number].append(fragment)
            source_position_list[number].append(alphabets[alphabet_count] + str(column_count))
            column_count += 1
            number += 1
        alphabet_count += 1
        
    return source_position_list

def optimal_reagents_volume(source_position_list, output_list, n_positions):
    
    min_DNA_volume = source_position_list
    
    # There needs to be min. 5uL of each reagent for aspiration
    total_buffer_volume = 5
    total_enzyme_volume = 5
    total_H2O_volume    = 5
    
    #Volume for each variant
    total_volume  = 20
    buffer_volume = 2
    enzyme_volume = 2
    reagent_no    = 0
    for source_info in source_position_list:
        DNA_volume = 5                 # minimum volume to aspirate is 2uL
        for variant in output_list:
            for fragment in variant:
                if source_info[0] == fragment:
                    if n_positions<=5:  # if no. of DNA fragments per well is 5 or less
                        DNA = 2         # 2uL of each DNA fragment is dispensed
                    else:               # for 6 or more fragments
                        DNA = 1         # 1uL of each DNA fragment is dispensed
                    DNA_volume += DNA
        min_DNA_volume[reagent_no].append(str(DNA_volume)+'uL')
        reagent_no += 1
        
    H2O_volume = total_volume - DNA*n_positions - buffer_volume - enzyme_volume

    #Calculate total volume of H2O, Buffer and Enzymes required
    for i in range(len(output_list)):
        total_buffer_volume += buffer_volume
        total_enzyme_volume += enzyme_volume
        total_H2O_volume    += H2O_volume
        
    # Arrangement of other reagents in 24 source well
    other_reagents_volume = [['H20', 'H1',str(total_H2O_volume) + 'uL'], ['Buffer', 'H2',str(total_buffer_volume) + 'uL'],['Enzyme', 'H3',str(total_enzyme_volume) + 'uL']]

    return min_DNA_volume, other_reagents_volume

###############################################################################
#Arrangement in output plate
###############################################################################
'''
total_output generate the total number of variant combinations

reagents_fill determines how a reagent from a specific row is filled
E.g 4*3*2
Each A reagent will repeat 3*2 time before the next A reagent is filled
Each B reagent will repeat 2 times before the next B reagent is filled,
repeat the addition of B reagents till n_output is maximum
For the last row, addition of reagent at least once before the next reagent
is filled, till n_output is maximum
Use reagents_fill_order to append fragment in output_list

generate_variants_combination generates all variant combinations

all_variants_position arranges all variants in 96 well output

desired_variant_output_well tells user the location of the desired variants
e.g plate, well, variant combination
'''
def total_output(protocol):
    n_output = 1
    
    #Determine number of output
    for fragment_variance in protocol:
        n_output = n_output*fragment_variance
        
    return n_output

def reagents_fill(aligned_inputs, protocol, n_input_row):
    reagents_fill_order = []
    start = 0
    
    for i in range(n_input_row-1):
        fill = 1
        for i in range(start, n_input_row-1):
            fill = fill*protocol[i+1]
        start += 1
        reagents_fill_order.append(fill)

        
    reagents_fill_order.append(1) 
    
    return reagents_fill_order

def generate_variants_combinations(reagents_fill_order, aligned_inputs, protocol, n_input_row, n_output): 
    output_list = []
    
    for i in range(n_output):
        output_list.append([])
    
    for i in range(n_input_row):
        count = 0
        max_fill = 0
        change_fragment = 0
        while count<n_output:
            if max_fill<reagents_fill_order[i]:
                output_list[count].append(aligned_inputs[i][change_fragment])
                count += 1
                max_fill += 1
            elif change_fragment<protocol[i]-1: # change_fragment should not exceed fragments in respective row
                change_fragment += 1
                max_fill = 0
            else:
                change_fragment = 0
                max_fill = 0
    
    return output_list

def all_variants_position(n_output):
    global position_in_96_well

    plate_num = 1
    position_count = 0
    position_list = []
    for i in range(n_output):
        position_list.append([])
    
    for i in range(n_output):
        if position_count>95:
            plate_num = round(i +1/96 +0.5)
            position_count = i-96
        position_list[i].append('Plate ' + str(plate_num))
        position_list[i].append(position_in_96_well[i])
        position_count += 1
    
    return position_list
      
def desired_variant_output_well(all_inputs, output_list, position_list):
    
    desired_position_list = dict()

    # a copy of all_inputs to check match with desired variant
    # and remove it from all_inputs_check to prevent repetition
    all_inputs_check = all_inputs
    
    #Generate plasmid list to match plasmid number to output variant
    plasmids_list = []

    for i in range(len(all_inputs)):
        plasmids_list.append('Plasmid ' + str(i+1))
    
    # Match desired plasmid to position
    count = 0
    plasmid_num = 0
    plasmid_count = 0
    for variant in output_list:
        perm = permutations(variant)
        for i in list(perm):
            if list(i) in all_inputs:
                for inputs in all_inputs:
                    if inputs == list(i):
                        plasmid_num = plasmid_count
                plasmid_count += 1
                desired_position_list[str(plasmids_list[plasmid_num])]=str(position_list[count]),str(list(i))
                all_inputs_check.remove(list(i))
        count += 1
        
    return desired_position_list

if __name__ == '__main__':
    yamldict = san.read_yaml('steps_3.yaml')
    stages = san.sort_align(yamldict)
    stage = stages['Recombine'][0]
    
    check_stage(stage)
    
    assert inputs2ot([['2','1'],['3','1'],['2','4'],['5','4']])
    assert inputs2ot([['1','4','1','5'],['2','4','2','5'],['3','4','3','5']])
    assert inputs2ot([['GFP','KM','15A'],['GFP','AMP','15A'],['RFP','KM','15A'],['RFP','AMP','15A']])
    
    #script, source_position, desired_variants_position, other_reagents_volume = inputs2ot([['GFP','KM','15A'],['GFP','AMP','15A'],['RFP','KM','15A'],['RFP','AMP','15A']])
    script, source_position, desired_variants_position, other_reagents_volume = stage2ot(stage)
    
    print(f'Script:\n{script}\n')
    print(f'Source position:\n{source_position}\n')
    print(f'Desired variants position:\n{desired_variants_position}\n')
    print(f'Volume of other reagents:\n{other_reagents_volume}\n')
    
