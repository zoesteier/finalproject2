# Final Project

[![Build
Status](https://travis-ci.org/zoesteier/finalproject2.svg?branch=master)](https://travis-ci.org/zoesteier/finalproject2)

Skeleton for neural network project.

## assignment
1. Implement an 8x3x8 autoencoder.
2. Implement a neural network to identify RAP1 transcription factor binding sites.
#. Apply the RAP1 network to test data to predict binding sites.



## structure

Autoencoder contained in neuralnet.py.
RAP1 neural network contained in rap1net.py.
Output sequences and scores from test data are contained in TestOutput.txt.


## usage

To use the package, first run

```
conda install --yes --file requirements.txt
```

to install all the dependencies in `requirements.txt`. Then the package's
main function (located in `hw2skeleton/__main__.py`) can be run as
follows

```
python -m hw2skeleton -P data test.txt
```

## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.


## contributors

Original design by Scott Pegg. Refactored and updated by Tamas Nagy.
