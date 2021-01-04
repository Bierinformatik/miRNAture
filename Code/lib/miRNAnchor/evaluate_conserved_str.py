#!/usr/bin/env python3

import sys
import re
import os
import RNA

from Bio import AlignIO

def evaluate_if_no_structure(str):
    str_len = len(str)
    count_gaps = 0
    for i in str:
        match_gap = re.match(r'\.', i)
        if match_gap:
            count_gaps = count_gaps + 1
    if str_len == count_gaps:
        return 1
    return 0


def get_final_block(query_string, start, stop_symbol='Def'):
    block = query_string[start:]
    end = start
    for i in range(0, len(block) -1):
    #for character in list(block):
        if stop_symbol != 'Def':  # Here, to clean the first part
            if block[i] == stop_symbol:  # corresponds to the stop symbol
                pre_block = query_string[start:end]  # ((..((.(((........
                # Here clean the last points that correspond to the loop
                cleaned_block = clean_block_tail(pre_block, '(')
                return cleaned_block
            # Here include ( and .
            end = end + 1
        else:  # Here is the last block only clean # )))....)......
            cleaned_block = clean_block_tail(block, ')')
            return cleaned_block


def clean_block_tail(strg_to_clean, symbol):
    count = 0
    to_delete = "."
    for char in reversed(list(strg_to_clean)):
        if char == to_delete and char != symbol:  # Identify last one symbol
            count = count - 1
        else:
            if count == 0: #There is not elements to remove
                cleaned_block = strg_to_clean[:]
                return cleaned_block
            cleaned_block = strg_to_clean[:count]
            break
    return cleaned_block


def get_block_coord(complete_str, pattern):
    start = complete_str.find(pattern)
    length = len(pattern)
    end = (start + length) - 1  # Lookin for array length, this is 0-based
    return(start, end)


def generate_tuples_positions_matching(block5, block3):
    block5_list = list(block5)
    block3_list = list(block3)  # Reversed just to take advantage of index
    block3_list.reverse()
    temporal_positions = []
    j = 0
    for i in range(len(block5_list)):
        if j > (len(block3) - 1):
            break
        if block5_list[i] == '(' and block3_list[j] == ')':  # ( )
            k = len(block3_list) - j  # Reverse coordinates len(block) - pos
            positions = (i, k)
            temporal_positions.append(positions)
            j = j + 1
        elif block5_list[i] == '(' and block3_list[j] != ')':  # ( .
            for l in range(j,(len(block3) - 1),1):
                if block5_list[i] == '(' and block3_list[l] == ')':
                    k = len(block3_list) - l  # Reverse coord len(block)- pos
                    positions = (i, k)
                    temporal_positions.append(positions)
                    j = j + 1
                    break
                j = j + 1
                continue
        elif block5_list[i] != '(' and block3_list[j] != ')':  # . .
            j = j + 1
        elif block5_list[i] != '(' and block3_list[j] == ')':  # . )
            continue
        else:
            print("Position ", i, " and ", j, " does not fit")
    return temporal_positions



def count_matching_sites(block, pattern):
    block_as_list = list(block)
    count = 0
    for i in block_as_list:
        if i == pattern:
            count = count + 1
    return count


def analize_consensus(str_string, family):
    eval_empty = evaluate_if_no_structure(str_string)
    length_str = len(str_string)
    if eval_empty == 1:
        print("No_valid_NoStr")
        #print(family + " does not have a valid Secondary structure:\n"
        #      + str_string + "\n")
        #outNoValidFile = open("no_valid_str.txt", "a+")
        #outNoValidFile.write(family + "\n")
        sys.exit()
    # We defined 5 blocks on a conserved str from miRNAs:
    # 5' complement, miR, loop, miR* and complement 3'
    mature5pStart = str_string.find('(')  # (((....((....))..)))......
    mature5pEndblock = get_final_block(str_string, mature5pStart, ')')
    (start_mature5p, end_mature5p) = get_block_coord(str_string,
                                                     mature5pEndblock)
    # ))..)))........
    mature3pStart = str_string.find(')')
    mature3pEnd = get_final_block(str_string, mature3pStart)
    (start_mature3p, end_mature3p) = get_block_coord(str_string, mature3pEnd)
    # Complements
    # From the defined block, the first '.'
    if mature5pStart > 0 and start_mature5p > 0:
        complement5pStart = str_string[:mature5pStart - 1].find('.')
        complement5pEnd = start_mature5p - 1
    else:
        complement5pStart = 0
        complement5pEnd = 0
    if end_mature3p < length_str:
        complement3pStart = end_mature3p + 1  # From def. block,the first '.'
        complement3pEnd = (len(str_string) - 1)  # Array length 0-based
    else:
        complement3pStart = 0
        complement3pEnd = 0
    # Loop
    startLoop = end_mature5p + 1
    endLoop = start_mature3p - 1
    block5pComplement = (complement5pStart, complement5pEnd)
    mature5block = (start_mature5p, end_mature5p)
    loop = (startLoop, endLoop)
    mature3block = (start_mature3p, end_mature3p)
    block3pComplement = (complement3pStart, complement3pEnd)
    return (block5pComplement, mature5block, loop, mature3block,
            block3pComplement)


def convert_matching_positions(matching_positions, factor5, factor3):
    new_positions = []
    for i in matching_positions:
        new5 = i[0] + factor5
        new3 = i[1] + factor3
        pair = (new5, new3)
        new_positions.append(pair)
    return new_positions


def validate_folding_miRNA(block5, block3):
    number_match_5p = count_matching_sites(block5, '(')
    number_match_3p = count_matching_sites(block3, ')')
    if number_match_5p > 20 and number_match_3p > 20:
        if number_match_5p == number_match_3p:  # Test if steem loop present
            #temporal = generate_tuples_positions_matching(block5, block3)
            temporal = "Valid_All"
        else:
            temporal = "No_valid_No_match"
    else:
        temporal = "No_valid_Short"
    return temporal


if __name__ == '__main__':
    seq = sys.argv[1]
    if os.path.exists(seq):
        name = os.path.basename(seq)
        family = name.split(".")[0]
    else:
        #print("The required file: " + seq + " did not exists!")
        sys.exit()
    align = AlignIO.read(seq, "stockholm")
    # Obtain information from secondary structure
    structure = align.column_annotations['secondary_structure']  # Obtain SStr
    structure = structure[9:-9]
    # Obtain blocks from consensus miRNA
    (block5pComplement, mature5block, loop, mature3block,
     block3pComplement) = analize_consensus(structure, family)
    block5mature = structure[mature5block[0]:mature5block[1] + 1]
    block3mature = structure[mature3block[0]:mature3block[1] + 1]
    # Look mature matching and retrieve positions but not updated
    # respect to miRNA context
    matching_positions = validate_folding_miRNA(block5mature, block3mature)
    # Only looking matching for the mature blocks
    # Convert coordinates to real Secondary structure coordinates
    #start5 = mature5block[0]
    #start3 = mature3block[0]
    #matching_positions_context = convert_matching_positions(matching_positions,
    #                                                        start5, start3)
    #print(family, structure, matching_positions + "\n" + block5mature + "\n" + block3mature)
    #print(family + "\t" + structure + "\t" + matching_positions)
    print(matching_positions)
