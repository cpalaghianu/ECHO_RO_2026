"""Make the local package importable when scripts are run directly (python scripts/01_...py)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
