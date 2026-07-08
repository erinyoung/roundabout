import pytest
from pathlib import Path
import pandas as pd
from unittest.mock import patch, MagicMock

# Update this import path to match your actual package file name
from roundabout.run_pygenomeviz import (
    run_pygenomeviz_blast,
    run_pygenomeviz_mummer,
    run_pygenomeviz_pmauve,
    run_pygenomeviz_mmseqs,
)

# ---------------------------------------------------------
# Testing Refactored PyGenomeViz Modules
# ---------------------------------------------------------


@patch("roundabout.run_pygenomeviz.Blast")
@patch("roundabout.run_pygenomeviz.Genbank")
@patch("roundabout.run_pygenomeviz.initial_gv_canvas")
@patch("roundabout.run_pygenomeviz.render_links_and_save")
@patch("roundabout.run_pygenomeviz.AlignCoord")
def test_run_pygenomeviz_blast_gbff(
    mock_aligncoord,
    mock_render,
    mock_initial_canvas,
    mock_genbank,
    mock_blast,
    tmp_path,
):
    """Tests the protein BLAST module workflow with GBFF input tracking."""
    dummy_gbks = [Path("sample1.gbff"), Path("sample2.gbff")]
    out_path = tmp_path / "test_blast_gbff.png"
    pgv_opts = {"length_thr": 100, "identity_thr": 30.0}

    # Setup Mocks
    mock_gbk_obj = MagicMock()
    mock_genbank.return_value = mock_gbk_obj
    mock_canvas_instance = MagicMock()
    mock_initial_canvas.return_value = mock_canvas_instance

    mock_blast_instance = MagicMock()
    mock_blast.return_value = mock_blast_instance
    mock_blast_instance.run.return_value = ["fake_coord_1", "fake_coord_2"]
    mock_aligncoord.filter.return_value = ["filtered_coord_1"]

    # Execute
    run_pygenomeviz_blast(
        input_paths=dummy_gbks, out_path=out_path, file_type="gbff", pgv_opts=pgv_opts
    )

    # Assertions
    assert mock_genbank.call_count == 2
    mock_initial_canvas.assert_called_once_with(pgv_opts)
    # Ensure BLAST handles protein mode for GBFF files
    mock_blast.assert_called_once_with([mock_gbk_obj, mock_gbk_obj], seqtype="protein")
    mock_aligncoord.filter.assert_called_once_with(
        ["fake_coord_1", "fake_coord_2"], length_thr=100, identity_thr=30.0
    )
    mock_render.assert_called_once_with(
        mock_canvas_instance, ["filtered_coord_1"], out_path, pgv_opts, method="blast"
    )


@patch("roundabout.run_pygenomeviz.ProgressiveMauve")
@patch("roundabout.run_pygenomeviz.Genbank")
@patch("roundabout.run_pygenomeviz.initial_gv_canvas")
@patch("roundabout.run_pygenomeviz.render_links_and_save")
@patch("roundabout.run_pygenomeviz.AlignCoord")
def test_run_pygenomeviz_pmauve_gbff(
    mock_aligncoord,
    mock_render,
    mock_initial_canvas,
    mock_genbank,
    mock_pmauve,
    tmp_path,
):
    """Tests that progressiveMauve maps out tracking metrics properly."""
    dummy_gbks = [Path("sample1.gbff"), Path("sample2.gbff")]
    out_path = tmp_path / "test_pmauve.png"
    pgv_opts = {"length_thr": 500, "refid": 0}

    # Setup Mocks
    mock_gbk_obj = MagicMock()
    mock_gbk_obj.get_seqid2size.return_value = {"seq1": 1000}
    mock_genbank.return_value = mock_gbk_obj

    mock_canvas_instance = MagicMock()
    mock_initial_canvas.return_value = mock_canvas_instance

    mock_pmauve_instance = MagicMock()
    mock_pmauve.return_value = mock_pmauve_instance
    mock_pmauve_instance.run.return_value = ["fake_coord_1"]
    mock_aligncoord.filter.return_value = ["filtered_coord_1"]

    # Execute
    run_pygenomeviz_pmauve(
        input_paths=dummy_gbks, out_path=out_path, file_type="gbff", pgv_opts=pgv_opts
    )

    # Assertions
    mock_pmauve.assert_called_once_with([mock_gbk_obj, mock_gbk_obj], refid=0)
    # Critical Check: Verify progressiveMauve does NOT receive identity_thr during filter tracking
    mock_aligncoord.filter.assert_called_once_with(["fake_coord_1"], length_thr=500)
    mock_render.assert_called_once_with(
        mock_canvas_instance, ["filtered_coord_1"], out_path, pgv_opts, method="pmauve"
    )


@patch("roundabout.run_pygenomeviz.MMseqs")
@patch("roundabout.run_pygenomeviz.Genbank")
@patch("roundabout.run_pygenomeviz.initial_gv_canvas")
@patch("roundabout.run_pygenomeviz.render_links_and_save")
@patch("roundabout.run_pygenomeviz.AlignCoord")
def test_run_pygenomeviz_mmseqs_guard_rail(
    mock_aligncoord,
    mock_render,
    mock_initial_canvas,
    mock_genbank,
    mock_mmseqs,
    tmp_path,
):
    """Tests that MMseqs correctly triggers an early exit if file_type is 'fasta'."""
    dummy_fastas = [Path("sample1.fasta"), Path("sample2.fasta")]
    out_path = tmp_path / "test_mmseqs_fail.png"

    # Execute passing fasta instead of gbff
    run_pygenomeviz_mmseqs(
        input_paths=dummy_fastas, out_path=out_path, file_type="fasta", pgv_opts={}
    )

    # Assertions: Verify no engine initialization occurred
    mock_mmseqs.assert_not_called()
    mock_initial_canvas.assert_not_called()
    mock_render.assert_not_called()
