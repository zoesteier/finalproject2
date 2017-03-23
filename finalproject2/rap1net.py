# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 11:05:47 2017

@author: Zoë
"""
#%% This neural network will learn to identify the RAP1 transcription factor binding site from dna sequences.


import neuralnet
import os
import numpy as np

# Read in positive sequences and test sequences from the data folder.
# Generate a list of sequences as strings. Each sequence is 17 bp.
def readPosTest(name):
    with open(os.path.join("..","data",name)) as f:
        #os.path.join("..","data","name") # .. goes up one folder
        inseqs = f.readlines()       
    seqlist = [] # a list to contain the sequences
    
    #row = row number, line is a string
    for row, line in enumerate(inseqs):
        tfbs = line.strip()
        seqlist.append(tfbs)    
    return seqlist

    
# Read in negative sequences from the data folder.
# Generate a list of sequences as strings. Each sequence is 1000 bp.
def readNeg(name):
    with open(os.path.join("..","data",name)) as f:
        #os.path.join("..","data","name") # .. goes up one folder
        inseqs = f.readlines()
        
    seqlist = [] # a list to contain the sequences
    
    seq = '' # initialize sequence string
    for row, text in enumerate(inseqs):
        if text[0] == '>': # start new sequence
            if seq != '': # if not the first sequence in the file
                seqlist.append(seq)
            seq = '' # re-initialize sequence string
        else: # not a new sequence
            line = text.strip()
            seq += line        
    return seqlist

# Create the lists of sequences
Posseqs = readPosTest('rap1-lieb-positives.txt')
Testseqs = readPosTest('rap1-lieb-test.txt')
Negseqs = readNeg('yeast-upstream-1k-negative.fa')


# Convert DNA sequences to a format that can be used in the neural network. Express each base as follows:
#    A = 0001
#    C = 0010
#    G = 0100
#    T = 1000
# With this method, every four digits represents one base pair and the hamming distance between each base is equal.
# A 17 bp sequence will be converted to 4*17=68 digits.

def convertDNAinput(sequence):
    """ Input a DNA sequence as a string. Output a string with the sequence mapped to binary."""
    sequence = sequence.replace('A', '0001')
    sequence = sequence.replace('C', '0010')
    sequence = sequence.replace('G', '0100')
    sequence = sequence.replace('T', '1000')
    return sequence
    

# Neural network for RAP1 will have 68 inputs (one for each digit in the input) and one output (0 is not a binding site, 1 is a binding site) with some number of hidden nodes.

