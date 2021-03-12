# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 16:09:12 2021

@author: Xinling
"""

import itertools as it
import os
import yaml      as ya
from   pathlib import Path
import math

#You need this for accessing HiPAD from this folder
import addpath

#I have written a function for sorting the steps in this module
import HiPAD.sort_align as san

###############################################################################
#Globals
###############################################################################
alphabets = ['A','B','C','D','E','F','G','H']
position_in_96_well = ['A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12',
                       'B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12',
                       'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12',
                       'D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11','D12',
                       'E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12',
                       'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12',
                       'G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12',
                       'H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12']
    
def checks(stage):
    #First, make sure the stage is not empty
    if not stage:
        msg = f'Stage is empty and recombination step is not needed.'
        raise ValueError(msg)

    #We want to identify the type of protocol e.g. (1, 2, 2) or (2, 2, 2) etc.
    #This assumes that the number of positions is the SAME for each variant
    #Check this now
    n_positions = None
    for step in stage:
        if not n_positions:
            n_positions = len(step['inputs'])
        elif n_positions == len(step['inputs']):
            pass
        else:
            msg = f'The number of position is not the same for each variant.'
            raise NotImplementedError(msg)
    
    return n_positions
    
    #Extracts inputs        
def get_inputs(stage):
    all_inputs = [step['inputs'] for step in stage]        
    
    return all_inputs
    
def align_inputs(all_inputs):
    n_positions = 0
    for position in all_inputs:
        max_positions = len(position)
        if max_positions>n_positions:
            n_positions = max_positions  #caution, raise error or n_positions differ?
 
    #Start off with first set of inputs in aligned_inputs
    step0_inputs   = all_inputs[0]
    aligned_inputs = [[step0_inputs[i]] for i in range(n_positions)]
    indices = list(range(len(aligned_inputs)))

    #Append fragments to each position in aligned_inputs if 
    # the fragment is not in it
    for position in all_inputs[1:]:
        
        for i in indices:
            if position[i] not in aligned_inputs[i]:
                aligned_inputs[i].append(position[i])
            else: 
                pass
            
    print(f'aligned_inputs before: {aligned_inputs}\n')
    
    #Finally, sort each lst in aligned_inputs so we no longer care about the 
    #initial order
    aligned_inputs = [sorted(lst) for lst in aligned_inputs]
    aligned_inputs = sorted(aligned_inputs, key=len)
    

    return aligned_inputs

def protocol_type(aligned_inputs):
    
    protocol_type = tuple([len(lst) for lst in aligned_inputs])
    
    return protocol_type
 
def choose_script(protocol_type):
    global scripts
    script_info = scripts[protocol_type]     #add scripts later
    
    return script_info
    
    #edit the code in this function to return instead of printing
    #ensure column_count doesn't exceed 12
def map_to_source_well(aligned_inputs):
    global alphabets
    alphabet_count = 0
    for position in aligned_inputs:
        column_count = 1
        for fragment in position:
            print(f'Place {fragment} in source well {alphabets[alphabet_count]}{column_count}')
            column_count+=1
        alphabet_count+=1
        
def map_to_output_well(aligned_inputs,protocol):
    global alphabets
    global position_in_96_well
    
    #No. of rows of input in source plate
    n_input_row = len(protocol)
    n_output = 1
    
    #Determine number of output
    for fragment_variance in protocol:
        n_output = n_output*fragment_variance
    
    output_list = []
    for i in range(n_output):
        output_list.append([])
    
    #Determine the no. of times a reagent from a specific row is filled
    #E.g 4*3*2
    #Each A reagent will repeat 3*2 time before the next A reagent is filled
    #Each B reagent will repeat 2 times before the next B reagent is filled,
    #repeat the addition of B reagents till n_output is maximum
    #For the last row, addition of reagent at least once before the next reagent
    #is filled, til n_output is maximum
    reagents_fill_order = []
    start = 0
    for i in range(n_input_row-1):
        fill = 1
        for i in range(start,n_input_row-1):
            fill = fill*protocol[i+1]
        start+=1
        reagents_fill_order.append(fill)

        
    reagents_fill_order.append(1) 
    print('Reagent fill order:')
    print(reagents_fill_order)
    print('')
    
    #Use reagents_fill_order to append fragment in output_list
    

    for i in range(n_input_row):
        count = 0
        max_fill = 0
        change_fragment=0
        while count<n_output:
            if max_fill<reagents_fill_order[i]:
                output_list[count].append(aligned_inputs[i][change_fragment])
                count+=1
                max_fill+=1
            elif change_fragment<protocol[i]-1: #change_fragment should not excess fragments in respective row
                change_fragment+=1
                max_fill = 0
            else:
                change_fragment=0
                max_fill = 0
   

    print('output list')
    print(output_list)
    print('')
    
    
    #Determine plate number and arrange output_list in rows (check if this is correct later)
    plate_num = 1
    position_count = 0
    position_list = []
    for i in range(n_output):
        position_list.append([])
    
    for i in range(n_output):
        if position_count>95:
            plate_num = math.ceil(i+1/96)
            position_count=i-96
        position_list[i].append('Plate ' + str(plate_num))
        position_list[i].append(position_in_96_well[i])
        position_count+=1
     
        
    print('position list')
    print(position_list)
    print('')
    
    #Match output list to plasmid number, need to edit this later
    track = 0
    for variant in output_list:
        if variant in all_inputs:
            position_list[track].append(output_list[track])
        track+=1
        
    print('position list after')
    print(position_list)
    print('')
    
if __name__ == '__main__':
    #Use the function I have written to sort the steps
    yamldict = san.read_yaml('steps_1.yaml')
    stages   = san.sort_align(yamldict)#, _store=lambda temp, step: temp.append(step['step_num']))

    #For now let's start small
    #Imagine you are just working on a single stage
    #DEFINITION: A stage is a group of steps that can be done in parallel
    #And if they can be done in parallel, you can use opentrons.
    stage = stages['Recombine'][0]
    #print('stage')
    #print(stage)
    #print('')
    
    checks(stage)

    all_inputs = get_inputs(stage)
    
    #Test cases
    #all_inputs = [['1','4','1','5'],['2','4','2','5'],['3','4','3','5']]
    #all_inputs = [['2','1'],['3','1'],['2','4'],['5','4']]
    print(f'all_inputs: {all_inputs}\n')
    
    aligned_inputs = align_inputs(all_inputs)
    print(f'aligned_inputs after: {aligned_inputs}\n')
    
    protocol = protocol_type(aligned_inputs)
    print(f'protocol: {protocol}\n')
    
    print('Source plate')
    map_to_source_well(aligned_inputs)
    print('')
    
    print('Output plate')
    map_to_output_well(aligned_inputs,protocol)
    
