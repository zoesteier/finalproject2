# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 15:24:36 2017

@author: Zoë
"""

import numpy as np
#%%
# Implement an 8x3x8 autoencoder. This neural network should take a matrix input and return the same matrix as an output.

# First, represent the neural network as a list of layers, where each layer in the network is represented as a class with a weight matrix, bias vector, activation function, function's derivative.

class Layer:
    """
    A class to represent a layer of nodes in a neural network.
    """
    
    def __init__(self, nodes, inputLayer, nextLayerNodes = None):
        self.nodes = nodes # number of nodes in this layer
        self.inputLayer = inputLayer # the previous layer of the network
        self.nextLayerNodes = nextLayerNodes # number of nodes in next layer
        self.weights = self.computeWeights() # weight matrix
        self.bias = self.computeBias() # bias vector
        #self.activation = self.activationfunction(self.inputLayer) # output vector of activations, one for each node
    
    def computeWeights(self):
        """ Initialize weight matrix with small random values generated from a normal distribution with mean = 0 and sd = 0.01**2 in order to break symmetry."""
#        if type(self.inputLayer) == np.ndarray: # input layer is array of inputs
#            innodes = len(self.inputLayer)
#        else: # input layer is class Layer (a hidden layer)
#            innodes = len(self.inputLayer.activation) # number of nodes in input
#        weights = np.zeros((innodes,self.nodes)) # weight matrix
#        for i in range(innodes): # rows are input nodes
#            for j in range(self.nodes): # columns are current layer nodes
#                
#                # for testing, stop use random seed
#                np.random.seed(3)
#                weights[i][j] = np.random.normal(0,0.01**2) 
        # For weights associated with the previous layer
        if self.nextLayerNodes == None:
            weights = None # output layer has no weights
        else: # current layer is not the output (i.e it is the input or hidden layer)
            nextnodes = self.nextLayerNodes
            weights = np.zeros((self.nodes,nextnodes)) # weight matrix
            for i in range(self.nodes): # rows are input nodes
                for j in range(nextnodes): # columns are next layer nodes
                    # for testing, stop use random seed
                    #np.random.seed(3)
                    weights[i][j] = np.random.normal(0,0.01**2) 
            self.weights = weights
        return weights
    
    def computeBias(self):
        """ Initialize bias vector with small random values generated from a normal distribution with mean = 0 and sd = 0.01**2 in order to break symmetry."""     
#        bias = np.zeros(self.nodes) # vector length is number of nodes in current layer
#        for i in range(self.nodes): 
#            # for testing, stop use random seed
#            np.random.seed(3)
#            bias[i] = np.random.normal(0,0.01**2) 

        # For bias associated with next layer.
        if self.nextLayerNodes == None:
            bias = None # output layer has no weights
        else: # current layer is not the output (i.e it is the input or hidden layer)
            bias = np.zeros([self.nextLayerNodes,1]) # vector length is number of nodes in next layer
            for i in range(self.nextLayerNodes): 
                # for testing, stop use random seed
                #np.random.seed(3)
                bias[i] = np.random.normal(0,0.01**2)

        return bias
        
    def activation(self, inputLayer, x):
        """ Compute the activation function for a given array of inputs. Use the sigmoid function as the activation function for forward propagation. f(x) = 1/(1 + e^-x). Derivative of the sigmoid function is f'(x) = f(x)*(1 - f(x))."""
        # If input is a Layer, not an array:
#        if type(inputLayer) == np.ndarray: # input layer is array of inputs
#            activation = np.transpose(np.array([inputLayer]))
#        else: # input layer is class Layer (a hidden layer)
#            #x = inputLayer.activation # array of activation functions coming from input
            
        # Rework activation function so that it evaluates the activation function for a given input layer. The current input (x) is the activation of the previous layer.
        # x is an array containing either the inputs to the network or the activation array of the previous layer.
    
    # weights and bias are from a layer are the output weights and bias of that layer. Calculate the input to the next layer using the weights and bias of the previous layer.

        z = np.dot(np.transpose(inputLayer.weights),x) + inputLayer.bias # z = Wtranspose*x + b, total weighted sum of inputs

        activation = 1/(1 + np.exp(-z)) # the sigmoid function
        return activation
        
        
    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return "{0} {1}".format(self.nodes, self.activation)

#%%
# Backpropagation algorithm to calculate delta values.
               
def backpropagation(network, x, y):
    """ Calculate delta values of a network for array inputs x and outputs y.
    Return gradw and gradb of the cost function as lists of matrices for each layer.
    """
    # Initialize lists to store gradwJ(W,b,x,y) = d(l+1)*a(l)t and gradbJ(W,b,x,y) = d(l+1) with each item in list corresponding to one layer from layer 0 to second to last (nl-1)
    gradw = [None]*(len(network)-1)
    gradb = [None]*(len(network)-1)
    
    # for output layer, d = -(y-a)*f'(z), f'(z) = a*(1-a)
    #a = network[-1].activation # get the activation of the output layer
    
    # Rework activation function. Feedforward pass to calculate activation.
    # Forward pass
    activations = [None]*(len(network)) # store activations, each item in list is a column vector of activations for that layer
    for layer in range(len(network)):
        if layer == 0: # for the input layer
            activations[layer] = x
        else: # for all other layers after input
            activations[layer] = network[layer].activation(network[layer-1],activations[layer-1])
    
    # Backward pass
    #a = network[-1].activationfunction(network[-1].inputLayer)
    # Rework for activation update
    a = activations[-1]
    endact = activations[-1] # store the activations of last layer
    fprime = a*(1-a) # f'(z) = a*(1-a)
    d = -(y-a)*fprime
    
    # for every other layer from the second to last to first, d = W*d(l+1)*f'(z)
    for layer in range(len(network)-2,-1,-1): # count backwards from second to last layer to the second layer (i.e. for 3 layers, this is just the hidden layer)
        #a = network[layer].activation # get the activation of the current layer
        # Rework a
        a = activations[layer]
        
        # calculate grad of current layer once we've calculated delta of next layer
        gradw[layer] = np.dot(a,np.transpose(d))
        gradb[layer] = d
        if layer == 0: # don't calculate delta for the first layer
            break
        
        # update fprime and delta for the current layer
        fprime = a*(1-a) # f'(z) = a*(1-a)
        d = np.dot(network[layer].weights,d)*fprime
        
        # collect the activations of the last layer for each input by returning a
        
    return gradw, gradb, endact
    

def gradientdescent(network, xmatrix, ymatrix, alpha = 0.5, weightdecay = 0.9):
    """ Use a linear regression cost function. J(W,b,x,y) = 0.5*(h(x) - y(x))**2 where x,y are true input and answer from training set and h is the predicted answer."""
    
    # set parameters
#    alpha = 0.5 # learning rate (start with 0.05, could try 0.5 or 0.1)
#    weightdecay = .9 # lambda (weight decay parameter, start with 0.9)
    
    # Initialize delW and delB with zeros of same dimensions as weight/bias at each layer.
    delW = [None]*(len(network)-1) # list to store weight matrix for each layer
    delB = [None]*(len(network)-1) # list to store bias vector for each layer
    for layer in range(len(network)-1): # for each layer except the last
        currentnodes = network[layer].nodes
        nextnodes = network[layer].nextLayerNodes
        delW[layer] = np.zeros([currentnodes, nextnodes])
        delB[layer] = np.zeros([nextnodes, 1])
        
    # Calculate delW and delB by adding the grad values for each input
    nodesperinput = np.shape(xmatrix)[0] # number of rows in input matrix
    nodesperoutput = np.shape(ymatrix)[0] # number of rows in output matrix
    assert nodesperinput == network[0].nodes # make sure input is correct size
    assert nodesperoutput == network[-1].nodes # make sure output is correct size
    m = np.shape(xmatrix)[1] # number of columns in the input matrix.
    assert m == np.shape(ymatrix)[1] # make sure same number of inputs and outputs
    
    # Store output values to compare to correct y values. Create a matrix of zeros with same shape as ymatrix.
    finalactivation = np.zeros_like(ymatrix)

    # Each column is one set of inputs.
    for i in range(m): # loop through each input
        x = np.zeros([nodesperinput, 1])
        y = np.zeros([nodesperoutput, 1])
        x[:,0] = xmatrix[:,i]
        y[:,0] = ymatrix[:,i]
        
        # calculate the deltas using backpropagation
        # actxy is the evaluated output of a given x,y input
        gradW, gradB, actxy = backpropagation(network, x, y)
        #print(x)
        #print(actxy)
        # update output values with new actxy values
#        
#        print(np.shape(actxy))
#        print(np.shape(finalactivation))
#        print(np.shape(ymatrix))
        #
        #print(finalactivation[:,i])

        
        #print(np.shape(finalactivation[:,i]))
        #print(type(finalactivation[:,i]))
        fat = np.transpose(finalactivation)
        #print(fat)
        #print(actxyt)
        #print(i)
        #print(np.transpose(actxy))
        fat[i,:] = np.transpose(actxy)
        #print(fat)
        #print('finalact')
        finalactivation = np.transpose(fat)
        #print(finalactivation)
        # problem with broadcasting
        #finalactivation[:,i] = actxy[:,:]
        
        for L in range(len(network)-1): # for each input, loop through each layer        
            delW[L] = delW[L] + gradW[L]
            delB[L] = delB[L] + gradB[L]
    
    # After computing the gradients, update weights and bias.
    for Layer in range(len(network)-1):
        newW = network[Layer].weights - alpha*((1/m)*delW[Layer] + weightdecay*delW[Layer]) # calculate new weight using weight decay term
        newB = network[Layer].bias - alpha*((1/m)*delB[Layer])
        # calculate new bias without weight decay term
        
        # update weights and bias
        network[Layer].weights = newW
        network[Layer].bias = newB
        
    # calculate the cost for this round of gradient descent       
    totalcost = None
    return totalcost, finalactivation


       
def autoencoder():
    #net = 9 #for testing
    
    # Make three layers for the network
    networkShape = np.array([3,2,3]) # number of nodes in each layer of the network
    inputsarray = np.array([1,0,0]) # input layer with 8 nodes
    inputs = Layer(networkShape[0], inputsarray, networkShape[1])
    hidden = Layer(networkShape[1], inputs, networkShape[2]) # hidden layer with 3 nodes, takes in inputs layer
    outputs = Layer(networkShape[2], hidden) # output layer with 8 nodes, takes in hidden layer
    
    # Make the network as an array of layers
    network = np.array([inputs, hidden, outputs])
    
    
    # make a new network compatible with multiple inputs
    ###for input 8
    networkShape = np.array([8,3,8])
    inputs = Layer(networkShape[0], None, networkShape[1])
    hidden = Layer(networkShape[1], inputs, networkShape[2]) # hidden layer with 3 nodes, takes in inputs layer
    outputs = Layer(networkShape[2], hidden) 
    newnet = np.array([inputs, hidden, outputs])
    
    #print(inputs.weights)
    
#    gradW, gradB = backpropagation(newnet, np.array([[1],[0],[0]]),np.array([[1],[0],[0]]))
#    print(gradB)
#    print(network[0])
#    print(network[0].weights)
#    print(network[1])
#    print(network[1].weights)
#    print(network[2])

    
    # Use the idenity matrix as an input to test the autoencoder.
    identityinput3 = np.identity(3)
    identityinput8 = np.identity(8)
    #print(newnet[1].bias)
    x = np.array([[1],[0],[0],[0],[0],[0],[0],[0]])
#    testcost, final = gradientdescent(newnet, x, x)
 
    # repeat steps of gradient descent
    for i in range(10000):
        testcost, finalactivation = gradientdescent(newnet, identityinput8, identityinput8)
 #       testcost, final = gradientdescent(network, x, x)
        #print(final)
#    print(newnet[0].weights)
    
    # need to initialize network w/o one input#################3
    #cost = gradientdescent(network, identityinput3, identityinput3)
    
   
#    activations = [None]*(len(newnet)) # store activations, each item in list is a column vector of activations for that layer
#        
#        # check output
#    for layer in range(len(newnet)):
#        if layer == 0: # for the input layer
#            activations[layer] = x
#            #print(x)
#        else: # for all other layers after input
#            activations[layer] = network[layer].activation(newnet[layer-1],activations[layer-1])
#    
#    a = activations[-1]
#    print(a)
#    print(finalactivation)
    
    return finalactivation

#final = autoencoder()
  
if __name__=="__main__":
    final = autoencoder()