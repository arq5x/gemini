"""Configuration YAML files for Gemini.

Provide Gemini configuration files in two alternative locations:

- Global:    /usr/local/share/gemini/gemini-config.yaml
- User only: $HOME/.gemini/gemini-config.yaml

Prefer the former if you have system level permissions for installation
since it will work for all system users.
"""
import os

import yaml

CONFIG_FILE = "gemini-config.yaml"
CONFIG_DIRS = [os.path.join(os.environ["HOME"], ".gemini"),
               "/usr/local/share/gemini"]

def _get_config_file(dirs=None):
    dnames = CONFIG_DIRS if dirs is None else dirs + CONFIG_DIRS
    for dname in dnames:
        fname = os.path.join(dname, CONFIG_FILE)
        if os.path.exists(fname):
            return fname
    raise ValueError("Gemini configuration file {0} not found in {1}".format(
            CONFIG_FILE, dnames))

def read_gemini_config(dirs=None, allow_missing=False):
    try:
        fname = _get_config_file(dirs)
    except ValueError:
        if allow_missing:
            return {}
        else:
            raise
    with open(fname) as in_handle:
        return yaml.load(in_handle)

def _find_best_config_file(dirs=None):
    dnames = CONFIG_DIRS if dirs is None else dirs + CONFIG_DIRS
    dnames.reverse()
    for dname in dnames:
        if os.access(dname, os.W_OK) or os.access(os.path.dirname(dname), os.W_OK):
            return os.path.join(dname, CONFIG_FILE)
    raise ValueError("Gemini configuration: "
                     "Could not find writeable directory: {0}".format(
            dnames))

def write_gemini_config(new_config, dirs=None):
    try:
        fname = _get_config_file(dirs)
    except ValueError:
        fname = _find_best_config_file(dirs)
    if not os.path.exists(os.path.dirname(fname)):
        os.makedirs(os.path.dirname(fname))
    with open(fname, "w") as out_handle:
        yaml.dump(new_config, out_handle, allow_unicode=False,
                  default_flow_style=False)
