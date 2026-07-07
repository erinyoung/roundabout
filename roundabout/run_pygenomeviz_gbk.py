import logging
from pathlib import Path
import pandas as pd
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank
from pygenomeviz.align import Blast, MUMmer, ProgressiveMauve, MMseqs, AlignCoord

# --- INTERNAL HELPERS ---

def _normalize_and_sort_gbks(gbk_paths: list[Path], local_matrix_df: pd.DataFrame) -> list[Genbank]:
    """
    Parses GenBank files, standardizes internal record IDs to filename stems 
    to prevent tool linkage errors, and sorts them by the ANI matrix.
    """
    raw_gbks = [Genbank(p) for p in gbk_paths]
    
    for gbk in raw_gbks:
        filename_stem = Path(gbk.file_path).stem if hasattr(gbk, "file_path") else gbk.name
        gbk._name = filename_stem
        if hasattr(gbk, "records"):
            for record in gbk.records:
                record.id = filename_stem
                for feature in record.features:
                    feature.ref = filename_stem

    if len(raw_gbks) < 3:
        return raw_gbks

    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform
    import numpy as np

    def clean_id(val):
        return Path(str(val)).stem.split('.')[0].strip()

    valid_gbks = []
    valid_ids = []
    for gbk in raw_gbks:
        cid = clean_id(gbk.name)
        if cid in local_matrix_df.index:
            valid_gbks.append(gbk)
            valid_ids.append(cid)

    if len(valid_gbks) < 3:
        return raw_gbks

    sub_matrix = local_matrix_df.loc[valid_ids, valid_ids].copy().fillna(80.0)
    distance_matrix = 100.0 - sub_matrix.values
    distance_matrix = (distance_matrix + distance_matrix.T) / 2.0
    np.fill_diagonal(distance_matrix, 0.0)

    Z = linkage(squareform(distance_matrix), method='average') 
    return [valid_gbks[i] for i in leaves_list(Z)]

def _build_gv_canvas(sorted_gbks: list[Genbank]) -> GenomeViz:
    """
    Initializes the GenomeViz canvas and paints the green CDS arrows.
    """
    gv = GenomeViz(track_align_type="center")
    gv.set_scale_bar()

    for gbk in sorted_gbks:
        track = gv.add_feature_track(gbk.name, gbk.get_seqid2size())
        for seqid, features in gbk.get_seqid2features("CDS").items():
            segment = track.get_segment(seqid)
            segment.add_features(features, plotstyle="bigarrow", fc="limegreen", lw=0.5)
            
    return gv

def _render_links_and_save(gv: GenomeViz, align_coords: list, out_path: Path, has_identity: bool):
    """
    Filters alignment strings, draws gray/red linkages, and saves the image.
    """
    if align_coords:
        color, inverted_color = "grey", "red"
        if has_identity and any(ac.identity for ac in align_coords):
            min_ident = int(min([ac.identity for ac in align_coords if ac.identity]))
            for ac in align_coords:
                gv.add_link(ac.query_link, ac.ref_link, color=color, inverted_color=inverted_color, v=ac.identity, vmin=min_ident)
            gv.set_colorbar([color, inverted_color], vmin=min_ident)
        else:
            for ac in align_coords:
                gv.add_link(ac.query_link, ac.ref_link, color=color, inverted_color=inverted_color)
    else:
        logging.warning("No alignment links passed filtering thresholds.")
        
    gv.savefig(out_path, dpi=300)

# --- DEDICATED API FUNCTIONS FOR ENGINE.PY ---

def run_pygenomeviz_gbk_blast(gbk_paths: list[Path], local_matrix_df: pd.DataFrame, out_path: Path, min_length: int = 100, min_identity: float = 30.0):
    sorted_gbks = _normalize_and_sort_gbks(gbk_paths, local_matrix_df)
    gv = _build_gv_canvas(sorted_gbks)
    try:
        align_coords = Blast(sorted_gbks, seqtype="protein").run()
        align_coords = AlignCoord.filter(align_coords, length_thr=min_length, identity_thr=min_identity)
        _render_links_and_save(gv, align_coords, out_path, has_identity=True)
    except Exception as e:
        logging.error(f"GenBank BLAST alignment failed: {e}")

def run_pygenomeviz_gbk_mummer(gbk_paths: list[Path], local_matrix_df: pd.DataFrame, out_path: Path, min_length: int = 100, min_identity: float = 30.0):
    sorted_gbks = _normalize_and_sort_gbks(gbk_paths, local_matrix_df)
    gv = _build_gv_canvas(sorted_gbks)
    try:
        align_coords = MUMmer(sorted_gbks).run()
        align_coords = AlignCoord.filter(align_coords, length_thr=min_length, identity_thr=min_identity)
        _render_links_and_save(gv, align_coords, out_path, has_identity=True)
    except Exception as e:
        logging.error(f"GenBank MUMmer alignment failed: {e}")

def run_pygenomeviz_gbk_pmauve(gbk_paths: list[Path], local_matrix_df: pd.DataFrame, out_path: Path, min_length: int = 100, min_identity: float = 30.0):
    sorted_gbks = _normalize_and_sort_gbks(gbk_paths, local_matrix_df)
    gv = _build_gv_canvas(sorted_gbks)
    try:
        align_coords = ProgressiveMauve(sorted_gbks).run()
        align_coords = AlignCoord.filter(align_coords, length_thr=min_length)
        _render_links_and_save(gv, align_coords, out_path, has_identity=False)
    except Exception as e:
        logging.error(f"GenBank progressiveMauve alignment failed: {e}")

def run_pygenomeviz_gbk_mmseqs(gbk_paths: list[Path], local_matrix_df: pd.DataFrame, out_path: Path, min_length: int = 100, min_identity: float = 30.0):
    sorted_gbks = _normalize_and_sort_gbks(gbk_paths, local_matrix_df)
    gv = _build_gv_canvas(sorted_gbks)
    try:
        align_coords = MMseqs(sorted_gbks).run()
        align_coords = AlignCoord.filter(align_coords, length_thr=min_length, identity_thr=min_identity)
        _render_links_and_save(gv, align_coords, out_path, has_identity=True)
    except Exception as e:
        logging.error(f"GenBank MMseqs alignment failed: {e}")