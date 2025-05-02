import unittest
import pandas as pd
from standard_morph.tools import orphan_node_check  

class TestOrphanNodeCheck(unittest.TestCase):

    def setUp(self):
        """Set up a dummy neuron morphology DataFrame."""
        data = [
            {"node_id": 1, "compartment": 1, "parent": -1, "parent_node_type": None, "x": 1, "y": 1, "z": 1},  # Root (soma)
            {"node_id": 2, "compartment": 3, "parent": 1,  "parent_node_type": 1,    "x": 2, "y": 2, "z": 2},  # Valid child
            {"node_id": 3, "compartment": 3, "parent": 1,  "parent_node_type": 1,    "x": 3, "y": 3, "z": 3},  # Valid child
            {"node_id": 4, "compartment": 3, "parent": 2,  "parent_node_type": 3,    "x": 4, "y": 4, "z": 4},  # Valid branch
            {"node_id": 5, "compartment": 3, "parent": 10, "parent_node_type": None, "x": 5, "y": 5, "z": 5},  # Orphaned
            {"node_id": 6, "compartment": 3, "parent": 20, "parent_node_type": None, "x": 6, "y": 6, "z": 6},  # Orphaned
        ]

        self.df = pd.DataFrame(data).set_index("node_id")
        self.df['node_id'] = self.df.index

    def test_orphan_nodes_detected(self):
        """Test case where orphan nodes exist."""
        result = orphan_node_check(self.df)[0]
        print(result)
        self.assertIsNotNone(result["nodes_with_error"], "Expected orphan nodes but none were detected.")
        self.assertIn((5,5,5,5), result["nodes_with_error"], "Node 5 should be detected as an orphan.")
        self.assertIn( (6,6,6,6), result["nodes_with_error"], "Node 6 should be detected as an orphan.")

    def test_no_orphan_nodes(self):
        """Test case where all nodes have valid parents."""
        df_valid = self.df.copy().drop([5, 6]) 
        result = orphan_node_check(df_valid)[0]
        self.assertIsNone(result["nodes_with_error"], "Expected no orphan nodes.")

if __name__ == "__main__":
    unittest.main()