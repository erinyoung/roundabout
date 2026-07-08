import pytest
import json
from pathlib import Path
from roundabout.database import check_bakta_db, check_plasmidfinder_db


def test_check_bakta_db_valid(tmp_path):
    """Tests if the Bakta checker accepts a properly formatted database folder."""
    db_dir = tmp_path / "db-light"
    db_dir.mkdir()

    # Create the required fake files
    (db_dir / "bakta.db").touch()
    with open(db_dir / "version.json", "w") as f:
        json.dump({"type": "light", "major": 5, "minor": 1}, f)

    # Should return the path since it is valid
    assert check_bakta_db(str(db_dir)) == str(db_dir)


def test_check_bakta_db_invalid(tmp_path):
    """Tests if the Bakta checker rejects a corrupted/incomplete database."""
    db_dir = tmp_path / "db-light"
    db_dir.mkdir()

    # We create bakta.db, but intentionally FORGET version.json
    (db_dir / "bakta.db").touch()

    # Should return None because it's missing files
    assert check_bakta_db(str(db_dir)) is None


def test_check_plasmidfinder_db_valid(tmp_path):
    """Tests if the PlasmidFinder checker looks for KMA indexes."""
    db_dir = tmp_path / "pf_db"
    db_dir.mkdir()

    # Create required files
    (db_dir / "VERSION").write_text("2.1")
    (db_dir / "config").touch()

    # Create the KMA index footprint
    (db_dir / "enterobacteriaceae.seq.b").touch()

    assert check_plasmidfinder_db(str(db_dir)) == str(db_dir)


def test_check_plasmidfinder_db_unindexed(tmp_path):
    """Tests if it rejects a downloaded but un-indexed PF database."""
    db_dir = tmp_path / "pf_db"
    db_dir.mkdir()

    (db_dir / "VERSION").write_text("2.1")
    (db_dir / "config").touch()

    # NO .seq.b files exist (user forgot to run INSTALL.py)
    assert check_plasmidfinder_db(str(db_dir)) is None
