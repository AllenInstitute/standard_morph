import unittest
import warnings
import pandas as pd
from standard_morph.tools import soma_and_soma_children_qc

class TestSomaQC(unittest.TestCase):

    def setUp(self):
        """Set up test cases with different SWC structures."""
        self.valid_swc = pd.DataFrame([
            {"node_id": 1, "compartment": 1, "x": 0, "y": 0, "z": 0, "r": 5, "parent": -1, "number_of_children": 1},
            {"node_id": 2, "compartment": 3, "x": 10, "y": 0, "z": 0, "r": 2, "parent": 1, "number_of_children": 1},
            {"node_id": 3, "compartment": 3, "x": 20, "y": 0, "z": 0, "r": 2, "parent": 2, "number_of_children": 0},
        ])#.set_index("node_id")

        self.multiple_somas = pd.DataFrame([
            {"node_id": 1, "compartment": 1, "x": 0, "y": 0, "z": 0, "r": 5, "parent": -1, "number_of_children": 1},
            {"node_id": 2, "compartment": 1, "x": 0, "y": 0, "z": 0, "r": 5, "parent": -1, "number_of_children": 0},
            {"node_id": 3, "compartment": 3, "x": 10, "y": 0, "z": 0, "r": 2, "parent": 1, "number_of_children": 1},
            {"node_id": 4, "compartment": 3, "x": 20, "y": 0, "z": 0, "r": 2, "parent": 3, "number_of_children": 0},
        ])#.set_index("node_id")

        self.no_soma = pd.DataFrame([
            {"node_id": 1, "compartment": 3, "x": 0, "y": 0, "z": 0, "r": 5, "parent": -1, "number_of_children": 1},
            {"node_id": 2, "compartment": 3, "x": 0, "y": 0, "z": 0, "r": 5, "parent": -1, "number_of_children": 0},
            {"node_id": 3, "compartment": 3, "x": 10, "y": 0, "z": 0, "r": 2, "parent": 1, "number_of_children": 1},
            {"node_id": 4, "compartment": 3, "x": 20, "y": 0, "z": 0, "r": 2, "parent": 3, "number_of_children": 0},
        ])#.set_index("node_id")

        self.soma_children_branching = pd.DataFrame([
            {"node_id": 1, "compartment": 1, "x": 0, "y": 0, "z": 0, "r": 5, "parent": -1, "number_of_children": 1},
            {"node_id": 2, "compartment": 3, "x": 10, "y": 0, "z": 0, "r": 2, "parent": 1, "number_of_children": 2},
            {"node_id": 3, "compartment": 3, "x": 20, "y": 0, "z": 0, "r": 2, "parent": 2, "number_of_children": 1},
            {"node_id": 4, "compartment": 3, "x": 30, "y": 0, "z": 0, "r": 2, "parent": 2, "number_of_children": 0},
        ])#.set_index("node_id")

        self.soma_children_exceed_distance = pd.DataFrame([
            {"node_id": 1, "compartment": 1, "x": 0, "y": 0, "z": 0, "r": 5, "parent": -1, "number_of_children": 1},
            {"node_id": 2, "compartment": 3, "x": 60, "y": 0, "z": 0, "r": 2, "parent": 1, "number_of_children": 1},  # Exceeds threshold
            {"node_id": 3, "compartment": 3, "x": 20, "y": 0, "z": 0, "r": 2, "parent": 2, "number_of_children": 0},
        ])#.set_index("node_id")

    def test_valid_swc(self):
        """Test that a valid SWC structure does not trigger errors."""
        errors = soma_and_soma_children_qc(self.valid_swc)
        for error in errors:
            self.assertIsNone(error["nodes_with_error"])
        
    def test_multiple_somas(self):
        """Test that multiple soma nodes are correctly detected."""
        with self.assertWarns(UserWarning) as cm:
            errors = soma_and_soma_children_qc(self.multiple_somas)

        self.assertIn("2 somas detected", str(cm.warning))
        soma_error = next(err for err in errors if err["test"] == "NumberOfSomas")
        self.assertIsNotNone(soma_error["nodes_with_error"])
        self.assertEqual(soma_error["nodes_with_error"], [(1, 0, 0, 0), (2, 0, 0, 0)])
        
    def test_no_somas(self):
        """Test that no soma nodes error is correctly detected."""
        with self.assertWarns(UserWarning) as cm:
            errors = soma_and_soma_children_qc(self.no_soma)
        self.assertIn("0 somas detected", str(cm.warning))
        soma_error = next(err for err in errors if err["test"] == "NumberOfSomas")
        self.assertIsNotNone(soma_error["nodes_with_error"])
        self.assertEqual(soma_error["nodes_with_error"], [(1,0,0,0)] )

    def test_soma_children_branching(self):
        """Test that improperly branching soma children are detected."""
        errors = soma_and_soma_children_qc(self.soma_children_branching)
        branching_error = next(err for err in errors if err["test"] == "SomaChildrenFurcation")
        self.assertIsNotNone(branching_error["nodes_with_error"])
        self.assertEqual(branching_error["nodes_with_error"], [(2, 10, 0, 0)])

    def test_soma_children_exceed_distance(self):
        """Test that soma children exceeding the distance threshold are detected."""
        errors = soma_and_soma_children_qc(self.soma_children_exceed_distance, soma_child_distance_threshold=50)
        distance_error = next(err for err in errors if err["test"] == "SomaChildrenDistance")
        self.assertIsNotNone(distance_error["nodes_with_error"])
        self.assertEqual(distance_error["nodes_with_error"] , [(2, 60, 0, 0)] )


if __name__ == '__main__':
    unittest.main()
