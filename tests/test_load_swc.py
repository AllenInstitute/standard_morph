import os
import unittest
import pandas as pd
from standard_morph.Standardizer import Standardizer  

class TestHortaLoadSwc(unittest.TestCase):

    def setUp(self):
        """Set up a dummy neuron morphology DataFrame."""
        test_dir = os.path.dirname(__file__)
        file_path = os.path.join(test_dir, "test_horta_swc.swc")
        
        self.worker_df = Standardizer(path_to_swc=file_path)._swc_df

    def test_hort_load(self):
        """Test case where orphan nodes exist."""
        loaded_xyz = self.worker_df[['x', 'y', 'z']]
        
        expected_coords = [
            [50.0, 50.0, 50.0],
            [60.0, 60.0, 60.0],
            [70.0, 70.0, 70.0],
        ]
        expected_df = pd.DataFrame(expected_coords, columns=['x', 'y', 'z'])

        pd.testing.assert_frame_equal(loaded_xyz, expected_df)

if __name__ == "__main__":
    unittest.main()
