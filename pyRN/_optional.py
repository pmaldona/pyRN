"""Helpers for optional runtime dependencies."""

from importlib import import_module


def require_dependency(module_name, package_name=None, extra=None):
    try:
        return import_module(module_name)
    except ImportError as exc:
        package = package_name or module_name.split(".")[0]
        if extra is None:
            message = f"{package} is required for this feature."
        else:
            message = (
                f"{package} is required for this feature. "
                f"Install it with: python -m pip install pyRN[{extra}]"
            )
        raise ImportError(message) from exc


def require_matplotlib_pyplot():
    return require_dependency("matplotlib.pyplot", "matplotlib", "viz")


def require_pypoman():
    return require_dependency("pypoman", "pypoman", "geometry")


def require_pyvis_network():
    return require_dependency("pyvis.network", "pyvis", "viz").Network


def require_roadrunner():
    return require_dependency("roadrunner", "libroadrunner", "roadrunner")
