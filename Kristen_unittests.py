import unittest
import numpy as np
import core_functions as cf
import wrapped_functions as wf
import wrappers as w
import pickle

source_image = 'pathname'
#use image C2 for testing, convert to numpy array, need both 3D and 2D arrays for testing
base_3D_image = cf.tiff_stack_2_np_arr(source_image)
base_2D_image = cf.max_projection(base_3D_image)

class Kristen_pipeline_tester(unittest.TestCase):
    def test_gamma_stabilize(self):
        #not the best one to start testing with since it has 4 conditions
        pass
    def test_robust_binarize(self):
        # important function to test
        base_image = cf.max_projection()
        pass



if __name__ == "__main__":
    unittest.main()