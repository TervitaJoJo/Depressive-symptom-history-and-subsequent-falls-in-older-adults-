"""Keep all generated data, caches and output outside the submission package."""
from pathlib import Path
PACKAGE_ROOT = Path(__file__).resolve().parents[2]
def external_path(value):
    path=Path(value).expanduser().resolve()
    if path==PACKAGE_ROOT or PACKAGE_ROOT in path.parents:
        raise ValueError("Use a private working directory outside the submission package.")
    return path
