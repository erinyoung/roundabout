import pytest
from pathlib import Path
from roundabout.run_grouping_by_amr import define_groups_by_amr
from roundabout.run_grouping_by_pf import define_groups_by_plasmidfinder


def test_define_groups_by_amr(tmp_path):
    """Tests if sequences with identical AMR profiles are grouped correctly."""
    # Dummy AMR dictionary
    amr_dict = {
        "isolate_A": ["blaNDM-5", "sul1"],
        "isolate_B": ["blaNDM-5", "sul1"],  # Should group with A
        "isolate_C": ["tet(A)"],  # Should be in its own group
        "isolate_D": [],  # Should be skipped (no genes)
    }

    # Run the grouping function (using tmp_path so the JSON saves safely)
    groups = define_groups_by_amr(amr_dict, outdir=tmp_path)

    # Assertions
    assert len(groups) == 2, "Should create exactly two groups"
    assert ["isolate_A", "isolate_B"] in groups, "A and B must group together"
    assert ["isolate_C"] in groups, "C must be isolated"
    assert not any("isolate_D" in g for g in groups), "D should be excluded"


def test_define_groups_by_amr_with_target(tmp_path):
    """Tests if the grouping properly filters by a target gene string."""
    amr_dict = {
        "isolate_A": ["blaNDM-5", "sul1"],
        "isolate_B": ["blaNDM-1", "tet(A)"],
        "isolate_C": ["tet(A)", "sul1"],
    }

    # Run targeting only "ndm"
    groups = define_groups_by_amr(amr_dict, outdir=tmp_path, target_gene="ndm")

    # Assertions
    assert len(groups) == 2, "A and B have NDM, C does not."
    assert ["isolate_A"] in groups
    assert ["isolate_B"] in groups
    # C is excluded completely because it lacks "ndm"
    assert not any("isolate_C" in g for g in groups)


def test_define_groups_by_plasmidfinder(tmp_path):
    """Tests if sequences with identical replicons are grouped correctly."""
    pf_dict = {
        "seq1": ["IncFII", "IncFIA"],
        "seq2": ["IncFII", "IncFIA"],  # Groups with seq1
        "seq3": ["IncM1"],  # Isolated
    }

    groups = define_groups_by_plasmidfinder(pf_dict, outdir=tmp_path)

    assert len(groups) == 2
    assert ["seq1", "seq2"] in groups
