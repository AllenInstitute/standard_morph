import unittest
import pandas as pd
from standard_morph.Standardizer import Standardizer

class TestStandardizer(unittest.TestCase):

    def setUp(self):
        """Set up test cases for validating inputs in Standardizer class."""
        self.valid_swc = "valid_path_to_swc.swc"
        self.invalid_swc = 12345  # Invalid path type (int instead of string)
        self.empty_swc = None
        
        # Valid DataFrame to simulate correct input morphology
        self.valid_df = pd.DataFrame({
            "node_id": [1, 2, 3],
            "compartment": [1, 3, 3],
            "x": [0, 10, 20],
            "y": [0, 0, 0],
            "z": [0, 0, 0],
            "r": [5, 2, 2],
            "parent": [-1, 1, 2]
        })
        
        # Invalid DataFrame with missing required columns
        self.invalid_df = pd.DataFrame({
            "node_id": [1, 2],
            "compartment": [1, 3],
            "x": [0, 10],
            "y": [0, 0],
            "r": [5, 2],  # Missing 'z' and 'parent' columns
        }) #.set_index("node_id")

        # Valid DataFrame that would be seen from aind_morphology_utils 
        self.valid_rename_df = pd.DataFrame.from_dict(
            {
                0: {"node_id":1, "struct_type": 1, "x": 0.0, "y": 0.0, "z": 0.0, "radius": 0.5, "parent_id": -1},
                1: {"node_id":2, "struct_type": 3, "x": 1.0, "y": 0.0, "z": 0.0, "radius": 0.3, "parent_id": 1},
                2: {"node_id":3, "struct_type": 3, "x": 1.0, "y": 1.0, "z": 0.0, "radius": 0.4, "parent_id": 1},
                3: {"node_id":4, "struct_type": 2, "x": 0.0, "y": 1.0, "z": 0.0, "radius": 0.6, "parent_id": 1},
            },
            orient="index",
        )


        self.valid_filename_format = "AIND"
        self.invalid_filename_format = "INVALID"
        self.missing_filename_format = "None"
        
    def test_valid_path_to_swc(self):
        """Test that a valid path_to_swc does not trigger errors."""
        with self.assertRaises(FileNotFoundError):
            standardizer = Standardizer(path_to_swc=self.valid_swc, input_morphology_df=None)
    
    def test_invalid_path_to_swc(self):
        """Test that invalid path_to_swc raises a ValueError."""
        with self.assertRaises(ValueError):
            standardizer = Standardizer(path_to_swc=self.invalid_swc, input_morphology_df=None)
    
    def test_empty_path_and_dataframe(self):
        """Test that both path_to_swc and input_morphology_df being empty raises an error."""
        with self.assertRaises(ValueError):
            standardizer = Standardizer(path_to_swc=self.empty_swc, input_morphology_df=None)



    def test_valid_dataframe(self):
        """Test that a valid dataframe passes without errors."""
        standardizer = Standardizer(path_to_swc=None, input_morphology_df=self.valid_df)
        self.assertEqual(standardizer.input_morphology_df.shape, self.valid_df.shape)

    def test_valid_rename_dataframe(self):
        """Test that a valid aind_morphology_utils dataframe passes without errors."""
        standardizer = Standardizer(path_to_swc=None, input_morphology_df=self.valid_rename_df)
        self.assertEqual(standardizer.input_morphology_df.shape, self.valid_rename_df.shape)

    def test_invalid_dataframe_columns(self):
        """Test that missing columns in the dataframe raises an error."""
        with self.assertRaises(ValueError):
            standardizer = Standardizer(path_to_swc=None, input_morphology_df=self.invalid_df)

    def test_rename_dataframe_columns(self):
        """Test that the dataframe columns are correctly renamed."""
        standardizer = Standardizer(path_to_swc=None, input_morphology_df=self.valid_df)
        expected_columns = ['node_id', 'compartment', 'x', 'y', 'z', 'r', 'parent']
        self.assertListEqual(list(standardizer.input_morphology_df.columns), expected_columns)

    def test_input_df_wrong_type(self):
        """Test that passing non-DataFrame as input_morphology_df raises ValueError."""
        with self.assertRaises(ValueError):
            Standardizer(path_to_swc=None, input_morphology_df={"x": [0, 1], "y": [1, 2]})

    def test_empty_dataframe_with_valid_columns(self):
        empty_df = pd.DataFrame(columns=self.valid_df.columns)
        with self.assertRaises(ValueError):
            Standardizer(path_to_swc=None, input_morphology_df=empty_df)


    def test_no_input(self):
        """Test that passing in nothing rasies an error."""
        with self.assertRaises(ValueError):
            standardizer = Standardizer(path_to_swc=None, input_morphology_df=None)


        
    def test_invalid_filename_format(self):
        """Test that an invalid filename format raises a ValueError."""
        with self.assertRaises(ValueError):
            standardizer = Standardizer(path_to_swc=None, input_morphology_df=self.valid_df, valid_filename_format=self.invalid_filename_format)

    def test_missing_filename_format(self):
        """Test that missing filename format does not trigger errors."""
        standardizer = Standardizer(path_to_swc=None, input_morphology_df=self.valid_df, valid_filename_format=self.missing_filename_format)
        self.assertEqual(standardizer.valid_filename_format, self.missing_filename_format)

if __name__ == '__main__':
    unittest.main()
