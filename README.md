# Final Project

[![Build
Status](https://travis-ci.org/zoesteier/finalproject2.svg?branch=master)](https://travis-ci.org/zoesteier/finalproject2)

Skeleton for neural network project.

## assignment
1. Implement Smith-Waterman algorithm.
2. Optimize a scoring matrix using a simulated annealing algorithm.


## structure

All functions are contained within SmithWaterman.py.


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
