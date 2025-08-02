# utils/config_loader.py
import os
import sys
import yaml

# Add the project root to Python path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

def load_config(config_path=os.path.join(project_root, 'config.yaml')):
    if not os.path.exists(config_path):
        raise FileNotFoundError(
            f"Configuration file '{config_path}' not found. "
            f"Please create one based on '{os.path.join(project_root, 'config.example.yaml')}'."
        )
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    return config

def get_tool_path(tool_name, config=None):
    if config is None:
        config = load_config()
    path = config.get("tools", {}).get(tool_name)
    if not path or not os.path.exists(path):
        raise FileNotFoundError(
            f"Path for '{tool_name}' not set or does not exist in config.yaml.\n"
            f"Set it under 'tools.{tool_name}'"
        )
    return path

print(get_tool_path('ovitosif'))