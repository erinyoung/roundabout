import pytest
import pandas as pd
from pathlib import Path
from roundabout.run_grouping_by_similarity import define_groups_by_similarity


def test_define_groups_by_similarity(tmp_path):
    """Tests if the core ANI network clustering groups identical plasmids."""

    # 1. Create a mock local ANI matrix (Square format)
    # A and B are 99% identical. C is an outlier (80%).
    data = {
        "sampleA_contig1": [100.0, 99.0, 80.0],
        "sampleB_contig1": [99.0, 100.0, 82.0],
        "sampleC_contig1": [80.0, 82.0, 100.0],
    }
    idx = ["sampleA_contig1", "sampleB_contig1", "sampleC_contig1"]
    local_matrix_df = pd.DataFrame(data, index=idx, columns=idx)

    # 2. Create a mock raw Skani DataFrame (Long format)
    raw_skani_data = {
        "query_name": ["sampleA_contig1", "sampleB_contig1", "sampleC_contig1"],
        "refseq_hit": ["sampleB_contig1", "sampleA_contig1", "sampleA_contig1"],
        "ani": [99.0, 99.0, 80.0],
        "align_fraction_query": [100.0, 100.0, 100.0],
        "align_fraction_ref": [100.0, 100.0, 100.0],
    }
    raw_skani_df = pd.DataFrame(raw_skani_data)

    # 3. Run the grouping logic
    groups = define_groups_by_similarity(
        raw_skani_df=raw_skani_df,
        local_matrix_df=local_matrix_df,
        min_ani=95.0,
        min_ani_align_fraction_query=90.0,
        min_ani_align_fraction_ref=90.0,
        outdir=tmp_path,
    )

    # 4. Assertions: Ensure A and B formed a group, and C was excluded
    assert len(groups) == 1, "Should create exactly one group for A and B"
    assert "sampleA_contig1" in groups[0]
    assert "sampleB_contig1" in groups[0]
    assert "sampleC_contig1" not in groups[0], "Outlier C should not be in the group"

    # Check that the JSON summary was saved
    assert (tmp_path / "skani_results" / "similarity_groups.json").exists()
