import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


@pytest.fixture
def fixture_content_root() -> Path:
    """A tiny synthetic skill+task tree, isolated from the real, human-authored
    civic-editorial-review/ content."""
    return FIXTURES_DIR / "content_root"


@pytest.fixture
def review_evidence_item_input() -> dict:
    """A realistic-but-fictional inputs dict matching the real review_evidence_item task's shape
    (an evidence_item, a list of pending revisions, and optional source_full_text)."""
    return json.loads((FIXTURES_DIR / "review_evidence_item_input.json").read_text())
