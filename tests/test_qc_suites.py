import unittest

from standard_morph.suites import resolve_suite, available_suites
from standard_morph.registry import REGISTRY
from standard_morph.models.qc_context import Space

CCF_ONLY = {"soma_inside_ccf_mesh", "nodes_outside_ccf_mesh"}


class TestSuiteReconciliation(unittest.TestCase):
    def setUp(self):
        self.pre = resolve_suite("default_pre_registration_tests")
        self.post = resolve_suite("default_post_registration_tests")

    def test_only_difference_is_the_ccf_metrics(self):
        # The rule: every check runs in both suites except the two CCF-only ones.
        self.assertEqual(set(self.post) - set(self.pre), CCF_ONLY)
        self.assertEqual(set(self.pre) - set(self.post), set())

    def test_pre_has_no_ccf_metrics(self):
        self.assertEqual(set(self.pre) & CCF_ONLY, set())

    def test_pre_is_a_prefix_of_post(self):
        self.assertEqual(self.post[: len(self.pre)], self.pre)

    def test_all_suite_metrics_are_registered(self):
        for name in set(self.pre) | set(self.post):
            self.assertIn(name, REGISTRY.names(), name)

    def test_ccf_metrics_are_the_only_ccf_restricted_ones_in_suites(self):
        # Any suite metric that is NOT compatible with image space must be one of
        # the CCF-only pair -- catches a future CCF metric being added to pre.
        for name in set(self.pre) | set(self.post):
            spaces = REGISTRY.get(name).applicability.spaces
            if Space.IMAGE_SPACE not in spaces:
                self.assertIn(name, CCF_ONLY, name)

    def test_available_suites(self):
        self.assertEqual(
            available_suites(),
            ["default_post_registration_tests", "default_pre_registration_tests"],
        )


if __name__ == "__main__":
    unittest.main()
