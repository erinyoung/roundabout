import importlib.metadata

try:
    # This must match the exact 'name' string defined in your pyproject.toml
    __version__ = importlib.metadata.version("roundabout")
except importlib.metadata.PackageNotFoundError:
    # Fallback for when the package is run locally without being installed
    __version__ = "0.1.0-dev"
