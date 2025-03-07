import unittest
import pandas as pd
from standard_morph.tools import orphan_node_check  

class TestOrphanNodeCheck(unittest.TestCase):

    def setUp(self):
        """Set up a dummy neuron morphology DataFrame."""
        data = [
            {"node_id": 1, "compartment": 1, "parent": -1, "parnet_node_type": None},  # Root (soma)
            {"node_id": 2, "compartment": 3, "parent": 1, "parnet_node_type": 1},  # Valid child
            {"node_id": 3, "compartment": 3, "parent": 1, "parnet_node_type": 1},  # Valid child
            {"node_id": 4, "compartment": 3, "parent": 2, "parnet_node_type": 3},  # Valid branch
            {"node_id": 5, "compartment": 3, "parent": 10, "parnet_node_type": 3},  # Orphaned (invalid parent)
            {"node_id": 6, "compartment": 3, "parent": 2, "parnet_node_type": None},  # Orphaned (missing parent type)
        ]

        self.df = pd.DataFrame(data).set_index("node_id")

    def test_no_orphan_nodes(self):
        """Test case where all nodes have valid parents."""
        df_valid = self.df.drop([5, 6])  # Remove orphaned nodes
        result = orphan_node_check(df_valid)
        self.assertIsNone(result["node_ids_with_error"], "Expected no orphan nodes.")

    def test_orphan_nodes_detected(self):
        """Test case where orphan nodes exist."""
        result = orphan_node_check(self.df)
        self.assertIsNotNone(result["node_ids_with_error"], "Expected orphan nodes but none were detected.")
        self.assertIn(5, result["node_ids_with_error"], "Node 5 should be detected as an orphan.")
        self.assertIn(6, result["node_ids_with_error"], "Node 6 should be detected as an orphan.")

    def test_orphan_nodes_fixed(self):
        """Test case where orphaned nodes are corrected."""
        self.df.at[5, "parent"] = 2  # Assign a valid parent
        self.df.at[6, "parnet_node_type"] = 1  # Assign a valid parent type
        result = orphan_node_check(self.df)
        self.assertIsNone(result["node_ids_with_error"], "Expected no orphan nodes after correction.")

if __name__ == "__main__":
    unittest.main()