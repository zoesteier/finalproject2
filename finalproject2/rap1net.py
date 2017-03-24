# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 11:05:47 2017

@author: Zoë
"""
#%% This neural network will learn to identify the RAP1 transcription factor binding site from dna sequences.


import neuralnet as nn
import os
import numpy as np
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

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

def convertDNA(sequence):
    """ Input a DNA sequence as a string. Output a string with the sequence mapped to binary."""
    sequence = sequence.replace('A', '0001')
    sequence = sequence.replace('C', '0010')
    sequence = sequence.replace('G', '0100')
    sequence = sequence.replace('T', '1000')
    return sequence
  
# Express converted DNA strings as column vectors so they can be used as inputs to the neural network.
def DNAtoInput(convertedStr):
    """ Input converted DNA sequence as a string and output a column vector."""
    assert len(convertedStr) == 68 # need 68 nodes to input to network
    inarray = np.zeros([68,1])
    for i, dig in enumerate(convertedStr): # i is row number, dig is digit
        inarray[i, 0] = dig    
    return inarray
    
#print(Posseqs[4])
#tryin = DNAtoInput(convertDNA(Posseqs[4]))
#print(np.shape(tryin))

def seqListtoMatInput(inseqlist):
    """Convert a list of sequences into a matrix input where each column is one sequence encoded as a 68x1 vector."""
    inmatrix = np.zeros((68,len(inseqlist)))
    for i, inseq in enumerate(inseqlist): # i = example index, inseq is string
        colinput = DNAtoInput(convertDNA(inseq)) # one input as a column vector
        inmatrix[:,i] = colinput[:,0]
    return inmatrix

# Divide positive and negative sequences for cross-validation
def divideData(pos, neg, k):
    """Input lists of all positive and negative sequences and number of sections for k-fold cross-validation. Output a list of k lists containing sets of data."""
    
    # Use all of the positive data and a subset of the negative data to result in a 1:1 ratio of positive:negative in each training and testing set.
    
    avg = len(pos) / float(k) # average number of sequences per set
    posdivided = []
    negdivided = []
    last = 0.0

    while last < len(pos):
        
        # make positive sequence list
        subset = pos[int(last):int(last + avg)] 
        # append list of sequences of the average length
        posdivided.append(subset)
        
        # make negative sequence lists the same length as positives for a 1:1 ratio
        # randomly choose negative sequences of length of subset
        negstart = np.random.randint(0, len(neg) - len(subset)) # index of first sequence
        negsubset = [] # store the subset of negative sequences
        # choose a random starting location to extract a 17bp neg sequence
        for i in range(len(subset)):
            sequence = negstart + i # index of sequence in neg list
            seqstart = np.random.randint(0, len(neg[sequence]) - 17)
            subseq = neg[sequence][seqstart:(seqstart+17)]
            negsubset.append(subseq) # append subseq to this set
        negdivided.append(negsubset) # append this set to list      
  
        last += avg    
    
    return posdivided, negdivided

#DivP, DivN = divideData(Posseqs, Negseqs, 1)
    
#%%
## Neural network for RAP1 will have 68 inputs (one for each digit in the input) and one output (0 is not a binding site, 1 is a binding site) with some number of hidden nodes.

# Make a neural network for RAP1 learning based on the autoencoder in neuralnet.py
def rap1net(shape, inseqlist, outbindvec):
    """Learns a neural network from inputs as sequences as a list of strings and outputs as a list of numbers where 0 means the sequence is not a RAP1 binding site and 1 means it is a binding site. Also output the network following learning.""" 
    # Make three layers for the network defined in shape.
    # shape = np.array([68,7,1]) # number of nodes in each layer of the network
    # Define the layers: Layer(nodes, inputLayer, nextLayerNodes = None)
    inputs = nn.Layer(shape[0], None, shape[1])
    hidden = nn.Layer(shape[1], inputs, shape[2])
    outputs = nn.Layer(shape[2], hidden)    
    net = np.array([inputs, hidden, outputs]) # this is the network
 
    # Convert a list of sequences into a matrix input where each column is one sequence encoded as a 68x1 vector.
    inmatrix = seqListtoMatInput(inseqlist)
    
    # repeat steps of gradient descent
    for i in range(10000):
        if i%10 == 0:
            print(i)
        cost, finalactivation = nn.gradientdescent(net, inmatrix, outbindvec, alpha = 0.5, weightdecay = 0.1)
  
    return finalactivation, net

# Try on one set
#DivP, DivN = divideData(Posseqs, Negseqs, 5)
#testin = np.hstack((np.array(DivP[0]), np.array(DivN[0])))
#testoutP = np.ones((1, len(DivP[0])))
#testoutN = np.zeros((1,len(DivN[0])))
#testout = np.hstack((testoutP, testoutN))
#testrap = rap1net(np.array([68,7,1]), testin, testout)

def forwardPass(network, x):
    """Given a column vector input (x), calculate output from a neural network."""
    
    activations = [None]*(len(network)) # store activations, each item in list is a column vector of activations for that layer
    for layer in range(len(network)):
        if layer == 0: # for the input layer
            activations[layer] = x
        else: # for all other layers after input
            activations[layer] = network[layer].activation(network[layer-1],activations[layer-1])
    a = activations[-1] # final activations are the output of the network
    return a

def crossValidate(pos, neg, k):
    """ Perform k-fold cross validation."""
    DivP, DivN = divideData(pos, neg, k) # divide pos and neg data into k sets
    predictoutput = [] # store all predicted values generated from the test data
    actualoutput = [] # store all true outputs in a list [1,1,1,...,0,0,0]
    
    for i in range(k): # iterate through each set as the test set
        print('Testing k = ')
        print(i)
        print('')
        
        # create input sequence list by concatenating pos and neg sequences
        # test set input array
        testin = np.hstack((np.array(DivP[i]), np.array(DivN[i])))
        # test set output array
        testoutP = np.ones((1, len(DivP[i])))
        testoutN = np.zeros((1,len(DivN[i])))
        testout = np.hstack((testoutP, testoutN)) # this is an nparray
        
        # training set input
        traininP = []
        traininN = []
        for j in range(k): # add each set except the test set to the training set
            if j !=i:
                traininP.extend(DivP[j])
                traininN.extend(DivN[j])

        trainin = np.hstack((np.array(traininP), np.array(traininN)))
        # training set output
        trainoutP = np.ones((1, len(traininP)))
        trainoutN = np.zeros((1,len(traininN)))
        trainout = np.hstack((trainoutP, trainoutN))
        
        # train the network on the training set, rap1net(shape, inseqlist, outbindvec)
        # this initializes a new network to train from scratch for each training set
        raptrainactivation, learnednet = rap1net(np.array([68,7,1]), trainin, trainout)
        
        
        # test on the test set        
        testinmatrix = seqListtoMatInput(testin) # convert seq list input to matrix
        # get the output from the network for each input in the matrix
        m = np.shape(testinmatrix)[1] # number of columns in the input matrix.
        predictout = [] # store output values in list

        for example in range(m): # loop through each input (one column of testinmatrix)
            x = np.zeros([68, 1])
            x[:,0] = testinmatrix[:,example]    
            raptestactivation = forwardPass(learnednet, x)
            val = raptestactivation[0,0]
            #print(np.shape(raptestactivation))
            predictout.append(val) # store the output from this example
        
        # collect results from this test set (for this value of k)
        predictoutput.extend(predictout) # append predicted outputs for the test set
        actualoutput.extend(testout.tolist()[0]) # append true values for the test set
        
    return predictoutput, actualoutput

#%%
#predict, actual = crossValidate(Posseqs, Negseqs,10)

def plotROC(actual, predictions):
    """Use sklearn's ROC function to plot ROC curves."""

#    FPrates = [[]]*len(SCORE_matrices)
#    TPrates = [[]]*len(SCORE_matrices)
#    AUCs = np.zeros((len(SCORE_matrices), 1)) # store all AUC values (area under the ROC curve)
    
    # Use sklearn's ROC curve function
#    for i in range(len(SCORE_matrices)):
    FPR, TPR, thresholds = roc_curve(actual, predictions)
    roc_auc = auc(FPR, TPR)
    print(roc_auc)
        
#        # collect values
#        FPrates[i] = false_positive_rate
#        TPrates[i] = true_positive_rate
#        AUCs[i] = roc_auc

    # plot ROC curve
    plt.title('ROC curve')
    plt.plot(FPR, TPR, label = ('k = 10, AUC = %0.5f'% roc_auc))
    plt.legend(loc='lower right')
    plt.plot([0,1],[0,1],'r--')
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')

    return FPR, TPR
    
#Predict, Actual = crossValidate(Posseqs, Negseqs, 10)
#%%
#plotROC(Actual, Predict)

# Evaluate test data set
def evaluateTestData(testseqlist, poslist, neglist):
    """ Evaluate the likelihood of a binding site in test data sequences."""
        
    # train the network on all available known data with best parameters  
    DivP1, DivN1 = divideData(Posseqs, Negseqs, 1) # use all training data in 1 set 
    # training set input
    trainin = np.hstack((np.array(DivP1[0]), np.array(DivN1[0])))
    # training set output
    trainoutP = np.ones((1, len(DivP1[0])))
    trainoutN = np.zeros((1,len(DivN1[0])))
    trainout = np.hstack((trainoutP, trainoutN))
    # train the network on the training set, rap1net(shape, inseqlist, outbindvec)
    # this initializes a new network to train from scratch for each training set
    raptrainactivation, bestnet = rap1net(np.array([68,7,1]), trainin, trainout)
    
    # test on the test set  
    testinmatrix = seqListtoMatInput(testseqlist) # convert test sequences to matrix input      
    # get the output from the network for each input in the matrix
    m = np.shape(testinmatrix)[1] # number of columns in the input matrix.
    predictout = [] # store output values in list

    for example in range(m): # loop through each input (one column of testinmatrix)
        x = np.zeros([68, 1])
        x[:,0] = testinmatrix[:,example]    
        raptestactivation = forwardPass(bestnet, x)
        val = raptestactivation[0,0]
        #print(np.shape(raptestactivation))
        predictout.append(val) # store the output from this example in a list
        
    # save the outputs in the correct format in a text file
    file = open('TestOutput.txt','w') 
    for seq in range(len(testseqlist)):
        file.write('\t'.join((str(testseqlist[seq]),str(predictout[seq]))))
        file.write('\n')     
    file.close()   
        
    return predictout

TestPrediction = evaluateTestData(Testseqs, Posseqs, Negseqs)
