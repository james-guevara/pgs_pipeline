"""Configuration loader for PGS pipeline."""

from pathlib import Path
from typing import Optional, Union

# Try Python 3.11+ tomllib, fall back to tomli
try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib
    except ImportError:
        tomllib = None


DEFAULT_CONFIG = {
    "paths": {
        "rsid_maps": "rsid_maps",
        "score_list": "sbayesrc_sumstats_filepaths.txt",
    },
    "directories": {
        "filter": "01a_filter",
        "annotate": "01b_annotate",
        "maf": "01c_maf",
        "merge": "01d_merge",
        "concat": "02_concat",
        "missingness": "03_missingness",
        "summary": "04_summary",
        "scores": "scores",
        "ancestry_maf": "06a_maf",
        "ancestry_prune": "06b_prune",
        "ancestry_relatedness": "06c_relatedness",
        "ancestry_pca": "06d_pca",
    },
    "defaults": {
        "threads": 16,
        "memory": 16000,
    },
    "output": {
        "format": "pgen",
    },
    "filters": {
        "maf": 0.01,
        "geno": 0.05,
        "mac": 10,
        "sample_miss": 0.05,
        "variant_miss": 0.05,
    },
    "ancestry": {
        "maf": 0.05,
        "king_cutoff": 0.0884,
        "num_pcs": 10,
        "ld_window": 200,
        "ld_step": 50,
        "ld_r2": 0.1,
    },
}


def load_config(config_path: Optional[Union[Path, str]] = None) -> dict:
    """Load configuration from TOML file, with defaults as fallback.

    Args:
        config_path: Path to config.toml. If None, looks for config.toml in
                     the same directory as this module, then current directory.

    Returns:
        Configuration dictionary with all settings.
    """
    config = DEFAULT_CONFIG.copy()

    if config_path is None:
        # Look for config.toml in module directory, then cwd
        module_dir = Path(__file__).parent
        candidates = [module_dir / "config.toml", Path("config.toml")]
        for candidate in candidates:
            if candidate.exists():
                config_path = candidate
                break

    if config_path and Path(config_path).exists():
        if tomllib is None:
            print("Warning: tomllib/tomli not available, using defaults")
            return config

        with open(config_path, "rb") as f:
            file_config = tomllib.load(f)

        # Deep merge file_config into config
        for section, values in file_config.items():
            if section in config and isinstance(config[section], dict):
                config[section].update(values)
            else:
                config[section] = values

    return config


# Convenience function to get a nested config value
def get(config: dict, *keys, default=None):
    """Get a nested config value safely.

    Example: get(config, "filters", "maf") -> 0.01
    """
    value = config
    for key in keys:
        if isinstance(value, dict) and key in value:
            value = value[key]
        else:
            return default
    return value
