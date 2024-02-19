#!/usr/bin/env python3
import sys, os
import argparse
import random

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Given an input file, extract a user specified number of random entries from this file and write to an output file")
    parser.add_argument('input', help = 'The input file', type = str)
    parser.add_argument('entries', help = 'The number of entries you want from the input file', type = int)
    parser.add_argument('output', help = 'The output file', type = str)

    args = parser.parse_args()
    input = args.input
    entries = args.entries
    output = args.output

    with open(input, 'r') as FI:
        input_list = [line.rstrip('\n') for line in FI.readlines()]

    rand_entries_list = list()
    i = 1
    while i <= entries:
        rand_entry = random.choice(input_list)
        if rand_entry not in rand_entries_list:
            rand_entries_list.append(rand_entry)
            i += 1

    with open(output, 'w') as FO:
        for entry in rand_entries_list:
            FO.write(f'{entry}\n')
