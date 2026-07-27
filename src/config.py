"""Load config.yaml and expose it as a plain dict."""
from pathlib import Path
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]


def load_config(path: str | Path = None) -> dict:
    path = Path(path) if path else REPO_ROOT / "config.yaml"
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    cfg["_repo_root"] = str(REPO_ROOT)
    return cfg


if __name__ == "__main__":
    import json
    print(json.dumps(load_config(), indent=2, default=str))
