import unittest
from pathlib import Path
import pandas as pd
from standard_morph._archived.Standardizer import Standardizer
from standard_morph._archived.tools import axon_origination_qc, axon_origin_distance_qc


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


class TestAxonOriginDistanceQC(unittest.TestCase):

    def setUp(self):
        """Load SWC fixture through Standardizer to compute parent relationships."""
        swc_path = Path(__file__).parent / "swcs" / "N024-648434-CONSENSUS.swc"
        self.standardizer = Standardizer(path_to_swc=str(swc_path))
        self.morph_df = self.standardizer.morph_df

    def test_fixture_axon_origin_distance_facts(self):
        """Test expected axon origin node, parent, and first edge distance."""
        axon_origins = self.morph_df[
            (self.morph_df["compartment"] == 2)
            & (self.morph_df["parent_node_type"] != 2)
        ]

        self.assertEqual(len(axon_origins), 1)
        axon_origin = axon_origins.iloc[0]
        self.assertEqual(axon_origin["node_id"], 93)
        self.assertEqual(axon_origin["parent"], 68)
        self.assertEqual(axon_origin["parent_node_type"], 3)
        self.assertEqual(axon_origin["parent_distance"], 5)

    def test_axon_origin_distance_passes_default_threshold(self):
        """Test axon origin distance passes default threshold."""
        result = axon_origin_distance_qc(self.morph_df)[0]

        self.assertEqual(result["test"], "AxonOriginDistance")
        self.assertIsNone(result["nodes_with_error"])

    def test_axon_origin_distance_fails_custom_threshold(self):
        """Test axon origin distance fails when threshold is below fixture edge length."""
        result = axon_origin_distance_qc(
            self.morph_df, axon_origin_distance_threshold=4
        )[0]

        self.assertEqual(
            result["nodes_with_error"], [(93, 30060.0, 9835.0, 12166.0)]
        )

    def test_standardizer_validate_includes_default_axon_origin_distance_check(
        self,
    ):
        """Test Standardizer validates axon origin distance with default threshold."""
        self.standardizer.write_all_tests_to_report = True
        self.standardizer.validate()

        axon_origin_distance_tests = [
            test
            for test in self.standardizer.StandardizationReport["tests"]
            if test["test"] == "AxonOriginDistance"
        ]

        self.assertEqual(len(axon_origin_distance_tests), 1)
        self.assertIsNone(axon_origin_distance_tests[0]["nodes_with_error"])


if __name__ == "__main__":
    unittest.main()
