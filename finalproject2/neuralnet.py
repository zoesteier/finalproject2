# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 15:24:36 2017

@author: Zoë
"""

import numpy as np
#%%
# Implement an 8x3x8 autoencoder. This neural network should take a matrix input and return the same matrix as an output.

# First, represent the neural network as a list of layers, where each layer in the network is represented as a class with a weight matrix, bias vector, activation function, function's derivative

class Layer:
    """
    A class to represent a layer of nodes in a neural network.
    """
    
    def __init__(self, nodes, inputLayer):
        self.nodes = nodes # number of nodes in this layer
        self.inputLayer = inputLayer # the previous layer of the network
        self.weights = self.computeWeights() # weight matrix
        self.bias = self.computeBias() # bias vector
        self.activation = self.activationfunction(self.inputLayer) # output vector of activations, one for each node
    
    def computeWeights(self):
        """ Initialize weight matrix with small random values generated from a normal distribution with mean = 0 and sd = 0.01**2 in order to break symmetry."""
        if type(self.inputLayer) == np.ndarray: # input layer is array of inputs
            innodes = len(self.inputLayer)
        else: # input layer is class Layer (a hidden layer)
            innodes = len(self.inputLayer.activation) # number of nodes in input
        weights = np.zeros((innodes,self.nodes)) # weight matrix
        for i in range(innodes): # rows are input nodes
            for j in range(self.nodes): # columns are current layer nodes
                
                # for testing, stop use random seed
                np.random.seed(3)
                weights[i][j] = np.random.normal(0,0.01**2) 
        #self.weights = weights
        return weights
    
    def computeBias(self):
        """ Initialize bias vector with small random values generated from a normal distribution with mean = 0 and sd = 0.01**2 in order to break symmetry."""     
        bias = np.zeros(self.nodes) # vector length is number of nodes in current layer
        for i in range(self.nodes): 
            # for testing, stop use random seed
            np.random.seed(3)
            bias[i] = np.random.normal(0,0.01**2) 
        #self.bias = bias
        return bias
        
    def activationfunction(self, inputLayer):
        """ Use the sigmoid function as the activation function for forward propagation. f(x) = 1/(1 + e^-x). Derivative of the sigmoid function is f'(x) = f(x)*(1 - f(x))."""
        if type(inputLayer) == np.ndarray: # input layer is array of inputs
            x = inputLayer
        else: # input layer is class Layer (a hidden layer)
            x = inputLayer.activation # array of activation functions coming from input
        z = np.dot(x,self.weights) + self.bias # z = x*W + b, total weighted sum of inputs
        activation = 1/(1 + np.exp(-z)) # the sigmoid function
        return activation
        
    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return "{0} {1}".format(self.nodes, self.activation)

#%%
# Backpropagation algorithm to update the weight matrix.

def costfunction(outputLayer):
    """ Use a linear regression cost function. J(W,b,x,y) = 0.5*(h(x) - y(x))**2 where x,y are true input and answer from training set and h is the predicted answer."""
    totalcost = None
    return totalcost
               
        
def autoencoder():
    net = 9 #for testing
    
    # Make three layers for the network
    inputs = np.array([1,0,0]) # input layer with 8 nodes
    hidden = Layer(2, inputs) # hidden layer with 3 nodes, takes in inputs layer
    outputs = Layer(1, hidden) # output layer with 8 nodes, takes in hidden layer
    
    # Make the network as an array of layers
    network = np.array([inputs, hidden, outputs])
    print(network[0])
    print(network[1].weights)
    print(network[1])
    print(network[2].weights)
    print(network[2])
    
    # Use the idenity matrix as an input to test the autoencoder.
    identityinput3 = np.identity(3)
    identityinput8 = np.identity(8)
    return net

autoencoder()