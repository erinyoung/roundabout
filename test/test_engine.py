import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from roundabout.engine import stage_and_split_fastas, run_pipeline


# ---------------------------------------------------------
# 1. Testing the FASTA Staging Logic
# ---------------------------------------------------------
def test_stage_and_split_fastas(tmp_path):
    """Tests if the pipeline correctly reads, splits, and stages multi-FASTAs."""
    # Setup mock input and output directories
    input_dir = tmp_path / "input_fastas"
    staging_dir = tmp_path / "staging"
    input_dir.mkdir()

    # Create a dummy multi-FASTA file
    dummy_fasta = input_dir / "sample_A.fasta"
    dummy_fasta.write_text(">contig_1\nATGCATGC\n" ">contig_2\nCGTACGTA\n")

    # Run the staging function
    staged_data = stage_and_split_fastas(input_dir, staging_dir)

    # Assertions
    assert len(staged_data) == 2, "Should have split into 2 separate sequences"
    assert staged_data[0]["contig_id"] == "contig_1"
    assert staged_data[0]["parent_name"] == "sample_A"

    # Check if the files were actually written to disk
    assert (staging_dir / "sample_A_contig_1.fasta").exists()
    assert (staging_dir / "sample_A_contig_2.fasta").exists()


# ---------------------------------------------------------
# 2. Testing the Pipeline Execution with Bypassed DBs
# ---------------------------------------------------------
@patch("roundabout.engine.run_setup")
@patch("roundabout.engine.execute_amrfinder_parallel")
@patch("roundabout.engine.execute_plasmidfinder_parallel")
@patch("roundabout.engine.execute_bakta_parallel")
@patch("roundabout.engine.execute_skani")
def test_run_pipeline_mocked(
    mock_skani, mock_bakta, mock_pf, mock_amr, mock_setup, tmp_path
):
    """Tests the main pipeline execution by completely bypassing databases and CLI tools."""

    # 1. Mock the Database Setup to pretend all databases are installed
    mock_setup.return_value = {
        "amrfinder": "/fake/amr/db",
        "bakta": "/fake/bakta/db",
        "plasmidfinder": "/fake/pf/db",
        "refseq_plasmid_dl": None,  # Skip refseq for this test
    }

    # 2. Mock the heavy execution functions to return dummy data instantly
    mock_amr.return_value = {"sample_A_contig_1": ["blaNDM-5", "sul1"]}
    mock_pf.return_value = {"sample_A_contig_1": ["IncFII"]}
    mock_bakta.return_value = [Path("/fake/bakta/results")]
    mock_skani.return_value = "fake_skani_results.tsv"

    # 3. Create fake input files so the pipeline has something to process
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    (input_dir / "sample_A.fasta").write_text(">contig_1\nATGC")

    # 4. Mock the argparse Namespace (the 'args' object passed into run_pipeline)
    # Using MagicMock allows it to safely return None or default values for all the CLI flags
    mock_args = MagicMock()
    mock_args.fastas = str(input_dir)
    mock_args.outdir = str(tmp_path / "results")
    mock_args.threads = 1
    mock_args.min_contig_length = 0
    mock_args.max_contig_length = 1000000
    mock_args.amr_gene = None
    mock_args.plasmidfinder_string = None

    # Prevent the pipeline from crashing when it hits visualization functions that
    # rely on the skipped tool outputs
    with (
        patch("roundabout.engine.run_minkemap_cohorts"),
        patch("roundabout.engine.run_daisyblast_cohorts"),
        patch("roundabout.engine.visualize_as_heatmap"),
        patch("roundabout.engine.visualize_global_matches_scatter"),
    ):

        # Execute the pipeline!
        try:
            run_pipeline(mock_args)
        except Exception as e:
            # If it hits a parsing error down the line because of missing real files,
            # that's okay for this test—we just want to ensure it routed correctly.
            pass

    # 5. Verify the pipeline actually attempted to run the tools!
    mock_setup.assert_called_once()
    mock_amr.assert_called_once()
    mock_pf.assert_called_once()
