import unittest
import pandas as pd
from collections import defaultdict, deque
from itertools import permutations
from standard_morph.tools import check_cycles_and_topological_sort  


class TestCheckCyclesAndTopologicalSort(unittest.TestCase):

    def setUp(self):
        """Set up dummy dataframes and child dictionaries for testing."""
        # Valid acyclic graph with multiple DFS possibilities
        self.df_valid = pd.DataFrame([
            {"node_id": 3, "parent": -1},  
            {"node_id": 2, "parent": 3},
            {"node_id": 1, "parent": 3},
            {"node_id": 4, "parent": 1},
            {"node_id": 5, "parent": 2},
        ])
        self.child_dict_valid = defaultdict(list, {
            3: [2, 1],  # The order in which 2 and 1 are visited can vary
            2: [5],
            1: [4],
            4: [],
            5: []
        })
        
        self.valid_dfs_orders = [
            {3:1, 2:2, 5:3, 1:4, 4:5},
            {3:1, 1:2, 4:3, 2:4, 5:5}
        ]
        
        # Graph with a cycle (invalid)
        self.df_cyclic = pd.DataFrame([
            {"node_id": 1, "parent": -1},  # Root
            {"node_id": 2, "parent": 1},
            {"node_id": 3, "parent": 2},
            {"node_id": 4, "parent": 3},
            {"node_id": 5, "parent": 4},
            {"node_id": 2, "parent": 5},  # Cycle (node 2 pointing back)
        ])
        self.child_dict_cyclic = defaultdict(list, {
            1: [2],
            2: [3],
            3: [4],
            4: [5],
            5: [2]  # Cycle here
        })

    def test_valid_graph(self):
        """Test that a valid acyclic graph returns a correct DFS-based topological sorting."""
        has_cycle, node_mapping = check_cycles_and_topological_sort(self.df_valid, self.child_dict_valid)
        
        self.assertEqual(cycle_report['node_ids_with_error'], None,  "Expected no cycles in the valid graph.")
        self.assertEqual(set(node_mapping.keys()), set(self.df_valid["node_id"]), "Node IDs should match.")
        self.assertEqual(len(node_mapping), len(self.df_valid), "All nodes should be assigned a new label.")

        # Get all valid DFS orderings
        valid_dfs_orders = self.valid_dfs_orders

        self.assertIn(node_mapping, valid_dfs_orders, "New node labels should match a valid DFS order.")

    def test_cyclic_graph(self):
        """Test that a cyclic graph is detected correctly."""
        cycle_report, node_mapping = check_cycles_and_topological_sort(self.df_cyclic, self.child_dict_cyclic)

        self.assertEqual(cycle_report['node_ids_with_error'], [1], "Expected a cycle to be detected.")
        self.assertEqual(node_mapping, {}, "Expected an empty mapping when a cycle is detected.")


if __name__ == "__main__":
    unittest.main()
