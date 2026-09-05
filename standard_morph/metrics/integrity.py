"""Input-integrity metrics (Phase 1).

These metrics read the raw SWC *table* (``swc_df``) -- an unvetted DataFrame,
not yet a graph -- and check that it is well-formed enough to build a faithful
``PreparedMorphology`` from. They are the reason a malformed file produces a
clean QC report ("duplicate node ids: [2]") instead of a crash: the engine runs
them *before* building the morphology, and each declares -- via
``blocks_on_failure`` (a ``BlockScope``) -- what its failure invalidates:

* ``required_columns``  (BUILD)    -- the canonical SWC columns are present
  (incl. radius); without them no arrays can be built at all, so a failure
  skips *every* morphology metric.
* ``non_empty``         (BUILD)    -- there is at least one node.
* ``castable_columns``  (BUILD)    -- every value casts to its column's numeric
  type (float coords, int ids); catches nulls and non-numeric/garbage values
  that would otherwise crash or silently corrupt the build.
* ``unique_node_ids``   (TOPOLOGY) -- node ids are unique, so the id -> index
  map is well-defined. A duplicate id corrupts the *topology* but not the
  coordinates, so a failure skips only the metrics that need topology;
  coordinate/attribute-only metrics still run.
* ``valid_parent_references`` (TOPOLOGY) -- every non-root ``parent`` id exists
  as a ``node_id`` (no dangling edges / orphaned subtrees).
* ``acyclic``           (TOPOLOGY) -- the parent chain has no cycle (incl. a
  self-parent), so tree walks (``segments``) cannot loop.
* ``parent_before_child`` (none) -- opt-in ordering report: each parent row
  precedes its children (a topological order). Non-blocking.

The first four plus ``valid_parent_references`` and ``acyclic`` always run (they
are ``PreparedMorphology``'s preconditions for a trustworthy tree);
``parent_before_child`` is opt-in.

Together these cover exactly ``PreparedMorphology.from_dataframe``'s hard
preconditions, so once ``required_columns`` and ``non_empty`` pass, building the
morphology is guaranteed safe.
"""
import time

import numpy as np
import pandas as pd

from standard_morph.metrics.base import Metric, Applicability, EvaluationPhase, BlockScope
from standard_morph.models.qc_context import ALL_COORDINATE_SPACES, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register
from standard_morph.preparation import SWC_COLUMN_NAMES

#: Columns that must hold integers (the rest of SWC_COLUMN_NAMES are floats).
INT_COLUMNS = ("node_id", "compartment", "parent")

#: Applicability shared by all integrity metrics: a malformed table is malformed
#: regardless of coordinate space, morphology kind, or resources.
_UNIVERSAL = Applicability(
    spaces=ALL_COORDINATE_SPACES,
    morphology_kinds=ALL_MORPHOLOGY_KINDS,
    required_resources=frozenset(),
)


class RequiredColumnsMetric(Metric):
    """Every canonical SWC column is present in the table.

    The framework works on the seven canonical columns
    ``node_id, compartment, x, y, z, r, parent`` (see ``SWC_COLUMN_NAMES``).
    ``PreparedMorphology`` builds one numpy array per column, so a missing column
    means an array cannot be built and *nothing* downstream can run -- hence this
    is a ``BUILD``-scope check: on failure the engine skips the entire morphology
    phase. Radius (``r``) is included, so a file without a radius column fails
    here rather than silently defaulting.

    Result
    ------
    * ``value`` / ``value_label`` -- ``n_missing_columns`` (0 means pass).
    * ``measurements`` -- ``required_columns`` / ``present_columns`` /
      ``missing_columns`` for a precise diff.
    """

    name = "required_columns"
    display_name = "Required SWC columns present"
    metric_number = None
    evaluation_phase = EvaluationPhase.INPUT_INTEGRITY
    blocks_on_failure = BlockScope.BUILD  # no columns -> no arrays -> nothing runs
    applicability = _UNIVERSAL

    def evaluate(self, swc_df, context, policy):
        """Fail if any canonical SWC column is absent from ``swc_df``."""
        t0 = time.perf_counter()
        result = MetricResult(name=self.name, status="pass")

        missing = [c for c in SWC_COLUMN_NAMES if c not in swc_df.columns]
        result.value = len(missing)
        result.value_label = "n_missing_columns"
        result.measurements = {
            "required_columns": list(SWC_COLUMN_NAMES),
            "present_columns": [str(c) for c in swc_df.columns],
            "missing_columns": missing,
        }
        result.counts = {"n_missing_columns": len(missing)}

        if missing:
            result.status = self.violation_severity.value
            result.message = f"missing required SWC column(s): {missing}"
        else:
            result.message = "All required SWC columns are present."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


class NonEmptyMetric(Metric):
    """The table contains at least one node.

    An empty table yields zero-length arrays -- there is no morphology to check
    at all -- so this is a ``BUILD``-scope check: on failure the morphology phase
    is skipped. (An SWC file that is all comments/headers, or a DataFrame that
    was filtered down to nothing, lands here.)

    Result
    ------
    * ``value`` / ``value_label`` -- ``n_rows`` (0 means fail).
    """

    name = "non_empty"
    display_name = "SWC table is non-empty"
    metric_number = None
    evaluation_phase = EvaluationPhase.INPUT_INTEGRITY
    blocks_on_failure = BlockScope.BUILD  # no rows -> no arrays -> nothing runs
    applicability = _UNIVERSAL

    def evaluate(self, swc_df, context, policy):
        """Fail if ``swc_df`` has no rows."""
        t0 = time.perf_counter()
        result = MetricResult(name=self.name, status="pass")

        n_rows = int(len(swc_df))
        result.value = n_rows
        result.value_label = "n_rows"
        result.measurements = {"n_rows": n_rows}
        result.counts = {"n_rows": n_rows}

        if n_rows == 0:
            result.status = self.violation_severity.value
            result.message = "SWC table has no rows."
        else:
            result.message = f"SWC table has {n_rows} node(s)."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


class UniqueNodeIdsMetric(Metric):
    """No two rows share a ``node_id``.

    The whole ``PreparedMorphology`` is indexed by mapping each node id to a row
    index; every ``parent`` pointer is resolved through that map. A duplicate id
    makes the map ambiguous (the last row silently wins), so parent resolution --
    and therefore the entire tree *topology* (``parent`` / ``children`` /
    ``roots`` / ``segments``) -- becomes untrustworthy. The per-node data
    (``xyz``, ``compartment``, ``r``) is *not* affected, though.

    Hence this is a ``TOPOLOGY``-scope check: on failure the engine skips only
    morphology metrics with ``requires_topology=True``; coordinate/attribute-only
    metrics (e.g. ``nodes_outside_ccf_mesh``) still run. The morphology can still
    be built -- it is just not safe to trust its structure.

    Result
    ------
    * ``value`` / ``value_label`` -- ``n_duplicate_ids`` (0 means pass).
    * ``flagged_node_ids`` -- the duplicated ids themselves (coordinates are
      ambiguous when an id is duplicated, so no coordinates are reported).
    """

    name = "unique_node_ids"
    display_name = "Node ids are unique"
    metric_number = None
    evaluation_phase = EvaluationPhase.INPUT_INTEGRITY
    blocks_on_failure = BlockScope.TOPOLOGY  # ids -> index map; corrupts topology only
    applicability = _UNIVERSAL

    def evaluate(self, swc_df, context, policy):
        """Fail if any ``node_id`` appears on more than one row."""
        t0 = time.perf_counter()
        result = MetricResult(name=self.name, status="pass")

        # Defensive: this normally runs only after required_columns has passed
        # (the engine orders and short-circuits phase 1), but guard so it never
        # crashes if invoked standalone or out of order.
        if "node_id" not in swc_df.columns:
            result.status = "error"
            result.message = "cannot check node id uniqueness: 'node_id' column is missing"
            result.runtime_ms = (time.perf_counter() - t0) * 1000
            return result

        counts = swc_df["node_id"].value_counts()
        duplicate_ids = [int(i) for i in counts[counts > 1].index.tolist()]

        result.value = len(duplicate_ids)
        result.value_label = "n_duplicate_ids"
        result.measurements = {"duplicate_ids": duplicate_ids, "n_duplicate_ids": len(duplicate_ids)}
        result.counts = {"n_duplicate_ids": len(duplicate_ids)}
        result.flagged_node_ids = duplicate_ids

        if duplicate_ids:
            result.status = self.violation_severity.value
            result.message = f"duplicate node_id(s): {duplicate_ids}"
        else:
            result.message = "All node ids are unique."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


class CastableColumnsMetric(Metric):
    """Every value casts cleanly to its column's required numeric type.

    The morphology is built by casting each column to a numpy array:
    ``x, y, z, r`` to float and ``node_id, compartment, parent`` to int. Without
    this check those casts fail in two bad ways:

    * a non-numeric string (e.g. a stray ``"NA"``) *crashes* the build with a
      ``ValueError``; and
    * a null / ``NaN`` in an int column silently becomes a garbage integer, while
      a null coordinate silently propagates ``NaN`` into every geometry metric.

    Validating the casts up front turns both failure modes into a clean report.
    It is ``BUILD`` scope: if a value cannot be cast, the arrays cannot be
    trusted, so the whole morphology phase is skipped. This subsumes a dedicated
    "no null values" check -- a null value fails to cast.

    A value is flagged when it is null, non-numeric, non-finite (``inf``), or --
    for the integer columns -- numeric but non-integral (e.g. ``1.5`` as a
    ``node_id``).

    Result
    ------
    * ``value`` / ``value_label`` -- ``n_uncastable_values`` (0 means pass).
    * ``measurements['problems_by_column']`` -- ``{column: count}`` of offending
      values, per column, so a reconstructor can see exactly where to look.
    """

    name = "castable_columns"
    display_name = "Column values cast to required types"
    metric_number = None
    evaluation_phase = EvaluationPhase.INPUT_INTEGRITY
    blocks_on_failure = BlockScope.BUILD  # uncastable values -> unusable arrays
    applicability = _UNIVERSAL

    def evaluate(self, swc_df, context, policy):
        """Fail if any value cannot be cast to its column's required type."""
        t0 = time.perf_counter()
        result = MetricResult(name=self.name, status="pass")

        problems = {}
        for col in SWC_COLUMN_NAMES:
            if col not in swc_df.columns:
                continue  # a missing column is required_columns' concern, not ours
            numeric = pd.to_numeric(swc_df[col], errors="coerce")
            # Null / non-numeric OR non-finite. `inf` casts to a legal float, so
            # to_numeric keeps it and isna() misses it -- but it poisons geometry
            # exactly like NaN, so isinf is flagged too.
            uncastable = numeric.isna() | np.isinf(numeric)
            if col in INT_COLUMNS:
                # numeric but non-integral (e.g. 1.5 as a node id) is also invalid
                uncastable = uncastable | ((~uncastable) & (numeric % 1 != 0))
            n_bad = int(uncastable.sum())
            if n_bad:
                problems[col] = n_bad

        total_bad = int(sum(problems.values()))
        result.value = total_bad
        result.value_label = "n_uncastable_values"
        result.measurements = {
            "problems_by_column": problems,
            "n_uncastable_values": total_bad,
        }
        result.counts = {"n_uncastable_values": total_bad}

        if problems:
            result.status = self.violation_severity.value
            result.message = (
                f"{total_bad} value(s) cannot be cast to their required numeric "
                f"type: {problems}"
            )
        else:
            result.message = "All column values cast cleanly to their required types."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


def _find_cyclic_nodes(id_to_parent):
    """Return the set of node ids lying on a cycle in the parent chain.

    ``id_to_parent`` maps node id -> parent id. A chain ends cleanly at ``-1`` (a
    root) or at a parent id absent from the map (a dangling reference -- not a
    cycle). Iterative (no recursion), so it is safe on very large morphologies.
    Each node is coloured once, so it runs in O(N).
    """
    SAFE, ONPATH, CYC = 1, 2, 3
    state = {}
    cyclic = set()
    for start in id_to_parent:
        if state.get(start) in (SAFE, CYC):
            continue
        path = []
        pos = {}
        node = start
        while True:
            st = state.get(node)
            if st in (SAFE, CYC):
                break  # merged into already-resolved territory
            if st == ONPATH:  # revisited a node on the current path -> a cycle
                i = pos[node]
                for n in path[i:]:
                    state[n] = CYC
                    cyclic.add(n)
                for n in path[:i]:
                    state[n] = SAFE
                path = []  # already resolved above
                break
            if node == -1 or node not in id_to_parent:
                break  # clean end: a root or a dangling reference
            state[node] = ONPATH
            pos[node] = len(path)
            path.append(node)
            node = id_to_parent[node]
        for n in path:  # nodes that reached a clean/safe end are safe
            if state.get(n) == ONPATH:
                state[n] = SAFE
    return cyclic


class ValidParentReferencesMetric(Metric):
    """Every non-root ``parent`` id exists as a ``node_id`` in the table.

    A parent id that is neither ``-1`` nor any node's id is a *dangling edge*:
    the node it names is missing, so that subtree is orphaned. ``from_dataframe``
    tolerates it (the node becomes an ``orphan`` / an extra component), but it is
    a genuine topology defect, so it is flagged here directly rather than only
    surfacing downstream as an unexpected connected component.

    ``TOPOLOGY`` scope: the coordinates are intact, but the tree wiring is not.
    """

    name = "valid_parent_references"
    display_name = "Parent references exist"
    metric_number = None
    evaluation_phase = EvaluationPhase.INPUT_INTEGRITY
    blocks_on_failure = BlockScope.TOPOLOGY
    applicability = _UNIVERSAL

    def evaluate(self, swc_df, context, policy):
        t0 = time.perf_counter()
        result = MetricResult(name=self.name, status="pass")
        if "node_id" not in swc_df.columns or "parent" not in swc_df.columns:
            result.status = "error"
            result.message = "cannot check parent references: 'node_id'/'parent' column missing"
            result.runtime_ms = (time.perf_counter() - t0) * 1000
            return result

        nid = pd.to_numeric(swc_df["node_id"], errors="coerce")
        par = pd.to_numeric(swc_df["parent"], errors="coerce")
        valid_targets = {int(v) for v in nid.dropna().tolist()}

        nid_list, par_list = nid.tolist(), par.tolist()
        dangling_rows = [
            i for i, p in enumerate(par_list)
            if pd.notna(p) and int(p) != -1 and int(p) not in valid_targets
        ]

        result.value = len(dangling_rows)
        result.value_label = "n_dangling_parents"
        result.flagged_node_ids = [int(nid_list[i]) for i in dangling_rows if pd.notna(nid_list[i])]
        result.measurements = {"n_dangling_parents": len(dangling_rows)}
        result.counts = {"n_dangling_parents": len(dangling_rows)}
        if dangling_rows:
            result.status = self.violation_severity.value
            bad = sorted({int(par_list[i]) for i in dangling_rows})
            result.message = (
                f"{len(dangling_rows)} node(s) reference a parent id absent from the "
                f"file (dangling edge); missing parent id(s): {bad[:10]}"
            )
        else:
            result.message = "All parent ids reference an existing node."
        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


class AcyclicMetric(Metric):
    """The parent chain contains no cycle.

    Each node has one ``parent`` pointer, so the ids form a functional graph;
    following parents must terminate at a root (``-1``) or a dangling reference.
    A cycle (including a self-parent, ``parent == node_id``) makes tree walks such
    as ``PreparedMorphology.segments`` loop, so it must be caught before any
    topology metric runs -- and because those walks are lazy and gated on
    ``requires_topology``, flagging the cycle here means they are never invoked.

    ``TOPOLOGY`` scope; detection is an iterative O(N) parent walk.
    """

    name = "acyclic"
    display_name = "Parent chain is acyclic"
    metric_number = None
    evaluation_phase = EvaluationPhase.INPUT_INTEGRITY
    blocks_on_failure = BlockScope.TOPOLOGY
    applicability = _UNIVERSAL

    def evaluate(self, swc_df, context, policy):
        t0 = time.perf_counter()
        result = MetricResult(name=self.name, status="pass")
        if "node_id" not in swc_df.columns or "parent" not in swc_df.columns:
            result.status = "error"
            result.message = "cannot check for cycles: 'node_id'/'parent' column missing"
            result.runtime_ms = (time.perf_counter() - t0) * 1000
            return result

        nid = pd.to_numeric(swc_df["node_id"], errors="coerce")
        par = pd.to_numeric(swc_df["parent"], errors="coerce")
        id_to_parent = {
            int(v): (int(p) if pd.notna(p) else -1)
            for v, p in zip(nid.tolist(), par.tolist()) if pd.notna(v)
        }

        cyclic = _find_cyclic_nodes(id_to_parent)

        result.value = len(cyclic)
        result.value_label = "n_cyclic_nodes"
        result.flagged_node_ids = sorted(cyclic)
        result.measurements = {"n_cyclic_nodes": len(cyclic)}
        result.counts = {"n_cyclic_nodes": len(cyclic)}
        if cyclic:
            result.status = self.violation_severity.value
            result.message = (
                f"parent chain contains a cycle involving {len(cyclic)} node(s): "
                f"{sorted(cyclic)[:10]}"
            )
        else:
            result.message = "Parent chain is acyclic."
        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


class ParentBeforeChildMetric(Metric):
    """Each node's parent appears before it in the file (a topological order).

    SWC files are conventionally written so a parent row precedes its children.
    It is *not* required for correctness here -- the morphology is indexed by id,
    not row order -- so this is an opt-in, non-blocking report. It checks
    parent-before-child (any topological order), which is deliberately weaker than
    a strict BFS/DFS order; the "soma is the first row" check lives in
    ``single_root_node``.
    """

    name = "parent_before_child"
    display_name = "Parent appears before child"
    metric_number = None
    evaluation_phase = EvaluationPhase.INPUT_INTEGRITY
    blocks_on_failure = None  # a file-ordering convention, not a precondition
    applicability = _UNIVERSAL

    def evaluate(self, swc_df, context, policy):
        t0 = time.perf_counter()
        result = MetricResult(name=self.name, status="pass")
        if "node_id" not in swc_df.columns or "parent" not in swc_df.columns:
            result.status = "error"
            result.message = "cannot check ordering: 'node_id'/'parent' column missing"
            result.runtime_ms = (time.perf_counter() - t0) * 1000
            return result

        nid = pd.to_numeric(swc_df["node_id"], errors="coerce").tolist()
        par = pd.to_numeric(swc_df["parent"], errors="coerce").tolist()
        first_row = {}
        for i, v in enumerate(nid):
            if pd.notna(v) and int(v) not in first_row:
                first_row[int(v)] = i

        out_of_order = []
        for i, p in enumerate(par):
            if pd.isna(p) or int(p) == -1:
                continue
            j = first_row.get(int(p))
            if j is not None and j >= i:  # parent at/after its child in the file
                out_of_order.append(i)

        result.value = len(out_of_order)
        result.value_label = "n_out_of_order"
        result.flagged_node_ids = [int(nid[i]) for i in out_of_order if pd.notna(nid[i])]
        result.measurements = {"n_out_of_order": len(out_of_order)}
        result.counts = {"n_out_of_order": len(out_of_order)}
        if out_of_order:
            result.status = self.violation_severity.value
            result.message = (
                f"{len(out_of_order)} node(s) appear before their parent in the file "
                "(not in parent-before-child order)."
            )
        else:
            result.message = "Every parent appears before its children."
        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(RequiredColumnsMetric())
register(NonEmptyMetric())
register(CastableColumnsMetric())
register(UniqueNodeIdsMetric())
register(ValidParentReferencesMetric())
register(AcyclicMetric())
register(ParentBeforeChildMetric())
