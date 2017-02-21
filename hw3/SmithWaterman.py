# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 11:24:32 2017

@author: Zoë
"""
import os
import glob
import numpy as np


#%%
# This function takes a score matrix in the data folder as an input and returns a numpy array containing the score of each amino acid pair

def readScoreMatrix(matrix): # matrix is name of matrix as a string
    with open(os.path.join("..","data",matrix)) as f:
        #os.path.join("..","data","BLOSUM50") # .. goes up one folder
        #print(f.readlines())
        input = f.readlines()
        #print(input[6:])
        
    matorder = input[6].strip().split() # list of amino acids, same for all matrices
    #print(matorder)
    #print(len(matorder)) # there are 20 amino acids and 4 additional characters
    
    scorematrix = np.zeros((len(matorder),len(matorder)))
    
    #row = 0 # initialize row number
    for row, matline in enumerate(input[7:]):
        #print(matline)
        line = matline.strip().split()
        #print(line)
        for column in range(len(line)):
            val = int(line[column])
            scorematrix[row, column] = val
    
    return scorematrix, matorder
    
BLOSUM50, matorder = readScoreMatrix('BLOSUM50')
BLOSUM62, matorder = readScoreMatrix('BLOSUM62')
PAM100, matorder = readScoreMatrix('PAM100')
PAM250, matorder = readScoreMatrix('PAM250')
#%%
# look up scores using a dictionary
# key is amino acid (string), value is position in the matrix (int)
scoredict = {}
for aa in range(len(matorder)):
    scoredict[matorder[aa]] = aa

# Look up the score for two amino acids
# Read in a sequence

#%%
# Read in sequences
def read_sequence(filepath):
    """
    Read in a single sequence given a .fa file

    """
    basename = os.path.basename(filepath)
    name = os.path.splitext(basename)[0]
    
    with open(filepath, "r") as f:
        input = f.readlines()
    seq = ''
    for row, text in enumerate(input[1:]):
        line = text.strip()
        seq += line
    #print(seq)
    #print(type(seq))
    return seq, name
    

dir = os.path.join("C:\\Users\\Zoë\\Documents\\GitHub\\hw3\\data\\sequences")
#files = glob.glob(dir + '/*.fa')
sequences = [] # a list of all sequences as strings
sequencenames = [] # a list of the names of each sequence (in order)
# iterate over each .fa file in the sequences directory
for filepath in glob.iglob(os.path.join(dir, "*.fa")):
    sequence, sequencename = read_sequence(filepath)
    sequences.append(sequence)
    #sequences.append(read_sequence(filepath))
    sequencenames.append(sequencename)
    # names are a list of strings, e.g. ['prot-0004','prot-0008',...]

print("Read in %d sequences"%len(sequences))  
#print(sequencenames)  
#print(sequences)
#print(sequences[0])

#%%
# Extract score from two characters
# To extract a score from the matrix, use score = BLOSUM50[scoredict[seq1[i]]][scoredict[seq2[j]]]

seq1 = sequences[0]
seq2 = sequences[1]

seq1char = scoredict[seq1[1]]
#print(seq1char)
seq2char = scoredict[seq2[1]]
#print(seq2char)
score = BLOSUM50[seq1char][seq2char]
#print(score)
# test: 0, 0: S, A score = 1
seq1char = scoredict[seq1[0]]
#print(seq1char)
seq2char = scoredict[seq2[0]]
#print(seq2char)
score = BLOSUM50[seq1char][seq2char]
#score1 = BLOSUM50[scoredict[seq1[0]]][scoredict[seq2[0]]]
#print(score)
#print(score1)
#print(score)
# test: 1, 1: L, N score = -4

#%%
# Smith Waterman algorithm
## Smith-Waterman algorithm

#!/Users/Rad/anaconda/bin/python
# (c) 2013 Ryan Boehning

#%% Implement Smith-Waterman algorithm
'''A Python implementation of the Smith-Waterman algorithm for local alignment
of nucleotide sequences.
'''

import argparse
import os
import re
import sys
import unittest


# Use scores from the scoring matrix.
SCORES = BLOSUM50
#gapopening = 10
#gapextension = 1
seq1     = None
seq2     = None
lastmove = None

def SmithWaterman(seq1, seq2):
#def main():
#    try:
#        seq1, seq2 = parse_cmd_line()
##        print(seq1)
##        print(seq2)
#    except ValueError as err:
#        print('error:', err)
#        return

    # Create a matrix to store scores (first sequence vertical, second horizontal, with an extra row and column for the gap character)
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    # Initialize the scoring matrix.
    score_matrix, start_pos, pointer_matrix = create_score_matrix(rows, cols, seq1, seq2)
    #print(score_matrix[start_pos[0]][start_pos[1]]) #confirm that find max score
    finalscore = score_matrix[start_pos[0]][start_pos[1]]

    # Traceback. Find the optimal path through the scoring matrix. This path
    # corresponds to the optimal local sequence alignment.
    seq1_aligned, seq2_aligned = traceback(score_matrix, start_pos, seq1, seq2, pointer_matrix)
    assert len(seq1_aligned) == len(seq2_aligned), 'aligned strings are not the same size'

#    # Print the results. The printing follows the format of BLAST results
#    # as closely as possible.
#    alignment_str, idents, gaps, mismatches = alignment_string(seq1_aligned, seq2_aligned)
#    alength = len(seq1_aligned)
#    print()
#    print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents,
#          alength, idents / alength, gaps, alength, gaps / alength))
#    print()
#    for i in range(0, alength, 60):
#        seq1_slice = seq1_aligned[i:i+60]
#        print('Query  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
#        print('             {0}'.format(alignment_str[i:i+60]))
#        seq2_slice = seq2_aligned[i:i+60]
#        print('Sbjct  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
#        print()
    
    return(finalscore)

def parse_cmd_line():
    '''Parse the command line arguments.
    Create a help menu, take input from the command line, and validate the
    input by ensuring it does not contain invalid characters (i.e. characters
    that aren't the bases A, C, G, or T).
    '''

    seq1 = sequences[0]
    seq2 = sequences[1]
    
# Example sequences
#    seq1 = "ATAGACGACATACAGACAGCATACAGACAGCATACAGA"
#    seq2 = "TTTAGCATGCGCATATCAGCAATACAGACAGATACG"
    return seq1, seq2


def create_score_matrix(rows, cols, seq1, seq2):
    '''Create a matrix of scores representing trial alignments of the two sequences.
    Sequence alignment can be treated as a graph search problem. This function
    creates a graph (2D matrix) of scores, which are based on trial alignments
    of different amino acids. The path with the highest cummulative score is the
    best alignment.
    '''
    #score_matrix = [[0 for col in range(cols)] for row in range(rows)]
    score_matrix = np.zeros((rows,cols)) # use numpy array for faster runtime
    pointer_matrix = np.zeros((rows,cols)) # store where each score came from in order to trace back the path later

    # Fill the scoring matrix.
    max_score = 0
    max_pos   = None    # The row and column of the highest score in matrix.
    
    # Calculate the score at each potential alignment
#    print(rows)
#    print(cols)
    for i in range(1, rows-1):
        for j in range(1, cols-1):
            #print([i,j])
            score = calc_score(score_matrix, i, j, seq1, seq2)
            if score > max_score:
                max_score = score
                max_pos   = (i, j)

            score_matrix[i][j] = score
            pointer_matrix[i][j] = lastmove

    assert max_pos is not None, 'the x, y position with the highest score was not found'
    #print(pointer_matrix)
    
    return score_matrix, max_pos, pointer_matrix


def calc_score(matrix, x, y, seq1, seq2):
    '''Calculate score for a given x, y position in the scoring matrix.
    The score is based on the up, left, and upper-left neighbors.
    '''
    
    #similarity = match if seq1[x - 1] == seq2[y - 1] else mismatch
    
    # The score at the new x,y position based on the SCORES matrix (e.g. BLOSUM50)
#    print(seq1)
#    print(SCORES[scoredict[seq1[x]]][scoredict[seq2[y]]])

    
    newposscore = SCORES[scoredict[seq1[x].upper()]][scoredict[seq2[y].upper()]]
    # use string.upper() to convert all AAs to upper case so their keys are recognized in the dictionary

    diag_score = matrix[x - 1][y - 1] + newposscore
    
    global lastmove # keep track of last move between rounds so we know whether the gap has already been opened
    if lastmove == 2: 
        up_score   = matrix[x - 1][y] - gapextension # extend gap vertically
        left_score = matrix[x][y - 1] - gapopening # open new gap horizontally
    elif lastmove == 3:
        up_score   = matrix[x - 1][y] - gapopening # open new gap vertically
        left_score = matrix[x][y - 1] - gapextension # extend gap horizontally
    else: # last move was diagonal or 0, so need to open gap
        up_score   = matrix[x - 1][y] - gapopening
        left_score = matrix[x][y - 1] - gapopening
    
    scorechoices = [0, diag_score, up_score, left_score]
    lastmove = scorechoices.index(max(scorechoices)) # update last move
    
    return max(0, diag_score, up_score, left_score)


def traceback(score_matrix, start_pos, seq1, seq2, pointer_matrix):
    '''Find the optimal path through the matrix.
    This function traces a path from the bottom-right to the top-left corner of
    the scoring matrix. Each move corresponds to a match, mismatch, or gap in one
    or both of the sequences being aligned. Moves are determined by the score of
    three adjacent squares: the upper square, the left square, and the diagonal
    upper-left square.
    WHAT EACH MOVE REPRESENTS
        diagonal: match/mismatch
        up:       gap in sequence 2
        left:     gap in sequence 1
    '''

    END, DIAG, UP, LEFT = range(4)
    aligned_seq1 = []
    aligned_seq2 = []
    x, y         = start_pos # position in the score matrix with the highest score
    # is the starting point for the alignment
    #move         = next_move(score_matrix, x, y) # This function traces the highest score through the matrix. Instead use the pointer matrix from where moves actually came from
    
    move = pointer_matrix[x][y]
    while move != END:
        if move == DIAG:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append(seq2[y - 1])
            x -= 1
            y -= 1
        elif move == UP:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append('-')
            x -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[y - 1])
            y -= 1

        move = next_move(score_matrix, x, y)

    aligned_seq1.append(seq1[x - 1])
    aligned_seq2.append(seq2[y - 1])

    # aligned from end to start, now reverse order of sequence to be from start to end
    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


def next_move(score_matrix, x, y):
    # I did not use this implementation (taken from another source), and instead created a pointer matrix keeping track of where each score came from to use for the traceback
    diag = score_matrix[x - 1][y - 1]
    up   = score_matrix[x - 1][y]
    left = score_matrix[x][y - 1]
    if diag >= up and diag >= left:     # Tie goes to the DIAG move.
        return 1 if diag != 0 else 0    # 1 signals a DIAG move. 0 signals the end.
    elif up > diag and up >= left:      # Tie goes to UP move.
        return 2 if up != 0 else 0      # UP move or end.
    elif left > diag and left > up:
        return 3 if left != 0 else 0    # LEFT move or end.
    else:
        # Execution should not reach here.
        raise ValueError('invalid move during traceback')


def alignment_string(aligned_seq1, aligned_seq2):
    '''Construct a special string showing identities, gaps, and mismatches.
    This string is printed between the two aligned sequences and shows the
    identities (|), gaps (-), and mismatches (:). As the string is constructed,
    it also counts number of identities, gaps, and mismatches and returns the
    counts along with the alignment string.
    AAGGATGCCTCAAATCGATCT-TTTTCTTGG-
    ::||::::::||:|::::::: |:  :||:|   <-- alignment string
    CTGGTACTTGCAGAGAAGGGGGTA--ATTTGG
    '''
    # Build the string as a list of characters to avoid costly string
    # concatenation.
    idents, gaps, mismatches = 0, 0, 0
    alignment_string = []
    for base1, base2 in zip(aligned_seq1, aligned_seq2):
        # zip function iterats over the i'th character in seq1 and seq2
        if base1 == base2:
            alignment_string.append('|')
            idents += 1
        elif '-' in (base1, base2):
            alignment_string.append(' ')
            gaps += 1
        else:
            alignment_string.append(':')
            mismatches += 1

    return ''.join(alignment_string), idents, gaps, mismatches


def print_matrix(matrix):
    '''Print the scoring matrix.
    ex:
    0   0   0   0   0   0
    0   2   1   2   1   2
    0   1   1   1   1   1
    0   0   3   2   3   2
    0   2   2   5   4   5
    0   1   4   4   7   6
    '''
    for row in matrix:
        for col in row:
            print('{0:>4}'.format(col))
        print()


class ScoreMatrixTest(unittest.TestCase):
    '''Compare the matrix produced by create_score_matrix() with a known matrix.'''
    def test_matrix(self):
        # From Wikipedia (en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
        #                -   A   C   A   C   A   C   T   A
        known_matrix = [[0,  0,  0,  0,  0,  0,  0,  0,  0],  # -
                        [0,  2,  1,  2,  1,  2,  1,  0,  2],  # A
                        [0,  1,  1,  1,  1,  1,  1,  0,  1],  # G
                        [0,  0,  3,  2,  3,  2,  3,  2,  1],  # C
                        [0,  2,  2,  5,  4,  5,  4,  3,  4],  # A
                        [0,  1,  4,  4,  7,  6,  7,  6,  5],  # C
                        [0,  2,  3,  6,  6,  9,  8,  7,  8],  # A
                        [0,  1,  4,  5,  8,  8, 11, 10,  9],  # C
                        [0,  2,  3,  6,  7, 10, 10, 10, 12]]  # A

        global seq1, seq2
        seq1 = 'AGCACACA'
        seq2 = 'ACACACTA'
        rows = len(seq1) + 1
        cols = len(seq2) + 1

        matrix_to_test, max_pos = create_score_matrix(rows, cols)
        self.assertEqual(known_matrix, matrix_to_test)


#if __name__ == '__main__':
#    sys.exit(main())

## Try an alignment
#seq1 = sequences[0]
#seq2 = sequences[1]
#final = SmithWaterman(seq1, seq2)
#print(final)

#%% 
# Question 1:
# #''' 
#    1. align all sequences and collect final score for each alignment >
#    2. separate scores for positive pairs and negative pairs >
#    3. find the threshold score at which 70% of true positives are above that score >
#    3. calculate false positive rate and true positive rate >
#    4. Vary gap openening (1 to 20) and gap extension( 1 to 5) with BLOSUM50 to find best gap penalty combination. >'''
 
def readSeqPairs(filepath):
# First need to generate a list of positive and negative sequences to align
    
    with open(filepath, "r") as f:
        input = f.readlines()
        # print(len(input)) # there are 50 pospairs and 50 negpairs
    seqpairs = [[]]*(len(input)) # a list of sequence pairs, each pair is a list of two names
    # characters 10-18 are seq1name, 33-41 are seq2name
    for row, text in enumerate(input[:]):
        seq1name = text[10:19]
        seq2name = text[33:42]
        seqpairs[row] = [seq1name, seq2name]

    return seqpairs

# Make lists of positive and negative sequence pairs
#Pospairs = readSeqPairs('C:\\Users\\Zoë\\Documents\\GitHub\\hw3\\data\\Pospairs.txt')
#Negpairs = readSeqPairs('C:\\Users\\Zoë\\Documents\\GitHub\\hw3\\data\\Negpairs.txt')
#print(Negpairs)

def alignSeqPairs(pairlist):
    '''Align sequence pairs contained in a list. Output a list of scores from those sequences.'''
    scorelist = np.zeros(len(pairlist)) # store scores in array, initialize to 0
    scorelistnorm = np.zeros(len(pairlist)) # store scores in array, initialize to 0
    
    # Get the sequence from the list of names
    for i in range(len(pairlist)):
        #print(pairlist[i])
        seq1name = pairlist[i][0]
        seq2name = pairlist[i][1]
        seq1 = sequences[sequencenames.index(seq1name)]
        seq2 = sequences[sequencenames.index(seq2name)]
        
        # Run the Smith-Waterman algorithm for sequence pair i
        score = SmithWaterman(seq1, seq2)
        scorenormalized = score/min(len(seq1), len(seq2)) # normalize score by dividing by the length of the shorter sequence
        
        scorelist[i] = score
        scorelistnorm[i] = scorenormalized
        
#    return scorelist
#    print(scorelist)
#    print(scorelistnorm)
    return scorelistnorm
    
## Make sorted lists of scores for positive and negative pairs
#Posscores = np.sort(alignSeqPairs(Pospairs))
#Negscores = np.sort(alignSeqPairs(Negpairs))

## Find the threshold value that 70% of positive scores fall above
#threshold = Posscores[int(np.ceil(0.3*len(Posscores)))]
## Find index of the first neg pair score above the threshold
#indabovethresh = np.where(Negscores > threshold)[0][0]
## Calculate false positive rate
#fprate = (len(Negscores) - indabovethresh)/len(Negscores)



# What is the best (lowest) false positive rate you can achieve by varying gap opening/extension penalties with BLOSUM50? Vary gap openening (1 to 20) and gap extension( 1 to 5) with BLOSUM50 to find best gap penalty combination.
def findBestPenalties():
    '''To find the best gap opening and extension penalites, vary the penalties and track lowest false positive rate.
    '''
    lowestfprate = 1 # initialize false positive rate (highest: 100% fp)
    
    # Make lists of positive and negative sequence pairs
    Pospairs = readSeqPairs('C:\\Users\\Zoë\\Documents\\GitHub\\hw3\\data\\Pospairs.txt')
    Negpairs = readSeqPairs('C:\\Users\\Zoë\\Documents\\GitHub\\hw3\\data\\Negpairs.txt')
    
#    global gapopening
#    global gapextension
    
    for gapopening in range(1, 21):
        for gapextension in range(1,6):
            print('Testing open = ', gapopening, 'and extension = ', gapextension)
            # Make sorted lists of scores for positive and negative pairs
            Posscores = np.sort(alignSeqPairs(Pospairs))
            Negscores = np.sort(alignSeqPairs(Negpairs))
            
            # Find the threshold value that 70% of positive scores fall above
            threshold = Posscores[int(np.ceil(0.3*len(Posscores)))]
            # Find index of the first neg pair score above the threshold
            indabovethresh = np.where(Negscores > threshold)[0][0]
            # Calculate false positive rate
            fprate = (len(Negscores) - indabovethresh)/len(Negscores)
            print('False positive rate = ', fprate)
            
            if fprate < lowestfprate:
                lowestfprate = fprate
                bestopen = gapopening
                bestextension = gapextension
            
    return bestopen, bestextension, lowestfprate
        
#gapopen, gapextension, fprate = findBestPenalties()
#print(gapopen), best = 9
#print(gapextension), best =  4
#print(fprate), best =  0.16

#%%
# Question 2:
    '''
    Plot ROC curve using gap penalties determined in question 1. Which scoring matrix performs best? 
    1. Find false positive rate at true positive rate = 0.7 for each scoring matrix >
    2. Which matrix performs best? >
    3. Create ROC curve (TP vs FP) showing plot from each matrix >
    '''

from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

# example for sklearn roc_curve
#actual = [1,1,1,0,0,0]
#predictions = [0.9,0.9,0.9,0.1,0.1,0.1]

gapopening = 9
gapextension = 4

def findBestMatrix():
    # Find the score matrix with the lowest false positive rate
    
    # Use gap scores determined above
    global SCORES
    SCORE_matrices = [BLOSUM50, BLOSUM62, PAM100, PAM250]
    
    # Find the matrix with lowest FP rate at TP rate = 0.7
    lowestfprate = 1 # initialize false positive rate (highest: 100% fp)
    
    # Make lists of positive and negative sequence pairs
    Pospairs = readSeqPairs('C:\\Users\\Zoë\\Documents\\GitHub\\hw3\\data\\Pospairs.txt')
    Negpairs = readSeqPairs('C:\\Users\\Zoë\\Documents\\GitHub\\hw3\\data\\Negpairs.txt')
    
    # store scores from all score matrices in an array
    allposscores = np.zeros((len(SCORE_matrices), 1))
    allnegscores = np.zeros((len(SCORE_matrices), 1))
    
    for i in range(len(SCORE_matrices)):
        SCORES = SCORE_matrices[i]
        print('Testing matrix: ', i)
        # Make sorted lists of scores for positive and negative pairs
        Posscores = np.sort(alignSeqPairs(Pospairs))
        Negscores = np.sort(alignSeqPairs(Negpairs))
        
        # add scores to allscores 
        allposscores[i] = Posscores
        allnegscores[i] = Negscores
        
        # Find the threshold value that 70% of positive scores fall above
        threshold = Posscores[int(np.ceil(0.3*len(Posscores)))]
        # Find index of the first neg pair score above the threshold
        indabovethresh = np.where(Negscores > threshold)[0][0]
        # Calculate false positive rate
        fprate = (len(Negscores) - indabovethresh)/len(Negscores)
        print('False positive rate = ', fprate)
        
        if fprate < lowestfprate:
            lowestfprate = fprate
            bestmatrixindex = i
            
    return lowestfprate, bestmatrixindex

#fp, bestmat = findBestMatrix(), fp = 0.12 for bestmat = 2 (PAM100)

def plotROC():
    # Plot ROC curve with values from each score matrix
    
    global SCORES
    SCORE_matrices = [BLOSUM50, BLOSUM62, PAM100, PAM250]
    
    # Make lists of positive and negative sequence pairs
    Pospairs = readSeqPairs('C:\\Users\\Zoë\\Documents\\GitHub\\hw3\\data\\Pospairs.txt')
    Negpairs = readSeqPairs('C:\\Users\\Zoë\\Documents\\GitHub\\hw3\\data\\Negpairs.txt')
    
    # store scores from all score matrices in an array
    allposscores = np.zeros((len(SCORE_matrices), len(Pospairs)))
    allnegscores = np.zeros((len(SCORE_matrices), len(Pospairs)))
    
    for i in range(len(SCORE_matrices)):
        SCORES = SCORE_matrices[i]
        print('Testing matrix: ', i)
        
        # Make sorted lists of scores for positive and negative pairs
        Posscores = np.sort(alignSeqPairs(Pospairs))
        Negscores = np.sort(alignSeqPairs(Negpairs))
        
        # add scores to allscores 
        allposscores[i] = Posscores
        allnegscores[i] = Negscores
        
    # concatenate all positive and negative scores to use as sklearn ROC input
    actualpos = np.ones(len(Posscores)) # true positive scores are "1'
    actualneg = np.zeros(len(Negscores)) # true negative scores are '0'
    actual = np.hstack((actualpos,actualneg)) # array of [1,1,1,1,...0,0,0]
    
    predictions = np.hstack((allposscores, allnegscores)) # predicted values are the scores output by Smith-Waterman, concatenate to [pos scores, ...neg scores] where each row is the result from one score matrix
    
#    FPrates = np.zeros((len(SCORE_matrices), len(actual))) # store all false positives
#    TPrates = np.zeros((len(SCORE_matrices), len(actual))) # store all true positives
    FPrates = [[]]*len(SCORE_matrices)
    TPrates = [[]]*len(SCORE_matrices)
    AUCs = np.zeros((len(SCORE_matrices), 1)) # store all AUC values (area under the ROC curve)
    
    for i in range(len(SCORE_matrices)):
        false_positive_rate, true_positive_rate, thresholds = roc_curve(actual, predictions[i])
        roc_auc = auc(false_positive_rate, true_positive_rate)
        
        # collect values
        FPrates[i] = false_positive_rate
        TPrates[i] = true_positive_rate
        AUCs[i] = roc_auc
    
    plt.title('Receiver Operating Characteristic (ROC Curve)')
    plt.plot(FPrates[0], TPrates[0], label = ('BLOSUM50, AUC = %0.2f'% AUCs[0]))
    plt.plot(FPrates[1], TPrates[1], label = ('BLOSUM62, AUC = %0.2f'% AUCs[1]))
    plt.plot(FPrates[2], TPrates[2], label = ('PAM100, AUC = %0.2f'% AUCs[2]))
    plt.plot(FPrates[3], TPrates[3], label = ('PAM250, AUC = %0.2f'% AUCs[3]))
    #label='AUC = %0.2f'% roc_auc)
    plt.legend(loc='lower right')
    plt.plot([0,1],[0,1],'r--')
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    return
    
plotROC()

#%%
# Question 3: Normalize scores

    

        