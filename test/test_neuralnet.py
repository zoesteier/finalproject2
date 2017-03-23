from finalproject2 import neuralnet
import numpy as np


def test_autoencoder():
    finalactivation = neuralnet.autoencoder()
    # For the 8x3x8 autoencoder, test to make sure the output is the 8x8 identity matrix.
    # Raise an assertion error if each element in finalactivation has a difference of > 0.1 from the identity matrix.
    np.testing.assert_array_almost_equal(finalactivation, np.identity(8), decimal=1)