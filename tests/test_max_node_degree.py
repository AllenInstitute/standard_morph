import unittest
import pandas as pd
from standard_morph._archived.tools import node_degree_check  


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
    
    def setUp(self):
        self._create_mock_dataframe()


    def test_node_degree_check_with_error(self):
        """
        Tests that a node exceeding the max_node_degree (default 4) is correctly identified.
        Node 2 has 5 children and is not the soma, so it should fail.
        """
        df = self._create_mock_dataframe()
        max_degree = 4
        results = node_degree_check(df, max_node_degree=max_degree)

        # Expected error node details for node_id=2
        expected_error = (2, 10.0, 1.0, 1.0)
        
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0]['test'], 'MaxNodeDegree')
        self.assertEqual(results[0]['description'], f'Nodes with more than {max_degree} children.')
        self.assertEqual(results[0]['nodes_with_error'], [expected_error])

    def test_node_degree_check_no_error(self):
        """
        Tests the passing case by setting max_node_degree high enough (6) 
        to ensure no errors are found.
        """
        df = self._create_mock_dataframe()
        max_degree = 5
        results = node_degree_check(df, max_node_degree=max_degree)

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0]['test'], 'MaxNodeDegree')
        self.assertEqual(results[0]['description'], f'Nodes with more than {max_degree} children.')
        self.assertIsNone(results[0]['nodes_with_error'])

if __name__ == "__main__":
    unittest.main()