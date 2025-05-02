import unittest
import pandas as pd
from standard_morph.tools import axon_origination_qc


class TestAxonOriginationQC(unittest.TestCase):

    def setUp(self):
        """Set up sample neuron morphology data."""
        self.sample_data = pd.DataFrame({
            'node_id': [1, 2, 3, 4, 5],
            'compartment': [1, 3, 2, 2, 2],  # 1: Soma, 3: Basal dendrite, 2: Axon
            'parent': [-1, 1, 1, 3, 3],
            'x': [0, 5, 10, 15, 20],
            'y': [0, 5, 10, 15, 20],
            'z': [0, 5, 10, 15, 20],
            'parent_node_type': [None, 1, 1, 2, 2]  # Parent node types
        }) #.set_index('node_id')

    def test_valid_axon_origination(self):
        """Test case where axon originates correctly from soma or basal dendrite."""
        result = axon_origination_qc(self.sample_data)[0]
        self.assertIsNone(result['nodes_with_error'], "Axon origination should be valid.")

    def test_invalid_multiple_axon_origins(self):
        """Test case where axon originates from multiple locations."""
        data = {
            'node_id': [1, 2, 3, 4, 5, 6],
            'compartment': [1, 3, 2, 2, 2, 2],  # Multiple axon origins
            'parent': [-1, 1, 1, 3, 1, 5],
            'x': [0, 5, 10, 15, 20, 25],
            'y': [0, 5, 10, 15, 20, 25],
            'z': [0, 5, 10, 15, 20, 25],
            'parent_node_type': [None, 1, 1, 3, 1, 2]
        }
        df = pd.DataFrame(data) #.set_index('node_id')
        result = axon_origination_qc(df)[0]
        
        self.assertIsNotNone(result['nodes_with_error'], "Should detect multiple axon origins.")
        self.assertGreater(len(result['nodes_with_error']), 1, "More than one invalid axon origin should be found.")

    def test_invalid_axon_origin_type(self):
        """Test case where axon originates from an invalid compartment."""
        data = {
            'node_id': [1, 2, 3, 4],
            'compartment': [1, 4, 2, 2],  # Axon originating from an invalid type (4)
            'parent': [-1, 1, 2, 2],
            'x': [0, 5, 10, 15],
            'y': [0, 5, 10, 15],
            'z': [0, 5, 10, 15],
            'parent_node_type': [None, 1, 4, 2]  # Invalid parent node type
        }
        df = pd.DataFrame(data) #.set_index('node_id')
        result = axon_origination_qc(df)[0]
        
        self.assertIsNotNone(result['nodes_with_error'], "Should detect invalid axon origin type.")
        self.assertEqual(result['nodes_with_error'], [(3,10,10,10)], "Only node 3 should be flagged as an invalid origin.")


if __name__ == '__main__':
    unittest.main()
