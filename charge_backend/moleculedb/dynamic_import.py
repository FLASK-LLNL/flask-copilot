import importlib.util
import sys
from types import ModuleType


def import_from_path(module_name: str, file_path: str) -> ModuleType | None:
    """
    Dynamically load a module from file path

    :param module_name: The name of the module to use
    :param file_path: The path to the Python file containing the module
    :return: Imported module, or None if import failed.
    """
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    if spec is None:
        return None
    module = importlib.util.module_from_spec(spec)

    # Add to sys.modules so it can be found by other imports
    sys.modules[module_name] = module

    # Execute the module
    assert spec.loader is not None
    spec.loader.exec_module(module)

    return module
