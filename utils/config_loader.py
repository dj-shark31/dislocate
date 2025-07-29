# utils/config_loader.py

import os
import yaml

def load_config(config_path="config.yaml"):
    if not os.path.exists(config_path):
        raise FileNotFoundError(
            f"Configuration file '{config_path}' not found. "
            f"Please create one based on 'config.example.yaml'."
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
