import pytest
from roundabout.cli import extract_version_from_text


def test_extract_version_from_text_standard():
    """Tests if the parser can grab standard semantic versions."""
    text = "blastn: 2.12.0+\nSome other text"
    version = extract_version_from_text(text, "blastn")
    assert version == "2.12.0+"


def test_extract_version_from_text_v_prefix():
    """Tests if the parser can grab 'v' prefixed versions."""
    text = "skani v0.2.1"
    version = extract_version_from_text(text, "skani")
    assert version == "v0.2.1"


def test_extract_version_from_text_unknown():
    """Tests the fallback for text with no numbers."""
    text = "command not found"
    version = extract_version_from_text(text, "fake_tool")
    assert version == "Unknown Version"


def test_extract_version_from_text_ignores_tool_name():
    """Tests that the parser doesn't just return the echoed tool name."""
    text = "minimap2\n2.26-r1175"
    version = extract_version_from_text(text, "minimap2")
    assert version == "2.26-r1175"
