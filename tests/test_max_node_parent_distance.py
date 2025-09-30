import unittest
import pandas as pd
from standard_morph.tools import distance_to_parent_node_check  


class TestQcChecks(unittest.TestCase):
    """
    Tests for quality control checks related to node degree and parent distance.
    """

    def _create_mock_dataframe(self) -> pd.DataFrame:
        """Helper method to create a mock DataFrame for testing QC functions."""
        data = {
            'node_id': [1, 2, 3, 4, 5, 6, 7],
            'parent_id': [-1, 1, 1, 2, 2, 3, 3],
            'compartment': [1, 3, 3, 3, 3, 3, 3],  # 1=Soma, 3=Dendrite
            'x': [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
            'y': [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            'z': [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            # Default values for children/distance (will be manually adjusted below)
            'number_of_children': [0] * 7,
            'parent_distance': [0.0, 50.0, 99.9, 100.1, 150.0, 100.0, 99.0]
        }
        df = pd.DataFrame(data)

        # Manually adjust children counts for specific test cases
        df.loc[df['node_id'] == 1, 'number_of_children'] = 6  # Soma node, degree 6 (ignored by QC)
        df.loc[df['node_id'] == 2, 'number_of_children'] = 5  # Node 2: 5 children (FAIL if max_degree=4)
        df.loc[df['node_id'] == 3, 'number_of_children'] = 4  # Node 3: 4 children (PASS if max_degree=4)

        return df


    def test_distance_to_parent_node_check_with_error(self):
        """
        Tests that nodes with parent_distance strictly greater than 100.0 
        (Nodes 4 and 5) are correctly identified.
        """
        df = self._create_mock_dataframe()
        max_distance = 100.0
        results = distance_to_parent_node_check(df, max_distance_to_parent_node=max_distance)

        # Expected error nodes: Node 4 (100.1) and Node 5 (150.0). Node 6 (100.0) passes.
        expected_errors = [
            (4, 30.0, 3.0, 3.0),
            (5, 40.0, 4.0, 4.0)
        ]

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0]['test'], 'MaxNodeDistanceFromParent')
        self.assertEqual(results[0]['description'], f'Nodes that are further than {max_distance} from their parent.')
        # Assert the found errors match the expected list
        self.assertEqual(results[0]['nodes_with_error'], expected_errors)

    def test_distance_to_parent_node_check_no_error(self):
        """
        Tests the passing case by setting max_distance high enough (200.0) 
        to ensure no errors are found.
        """
        df = self._create_mock_dataframe()
        max_distance = 200.0
        results = distance_to_parent_node_check(df, max_distance_to_parent_node=max_distance)

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0]['test'], 'MaxNodeDistanceFromParent')
        self.assertEqual(results[0]['description'], f'Nodes that are further than {max_distance} from their parent.')
        self.assertIsNone(results[0]['nodes_with_error'])

if __name__ == "__main__":
    unittest.main()