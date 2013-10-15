"""Configuration YAML files for Gemini.

Provide Gemini configuration files in alternative locations:

- Installer based: gemini-virtualenv/../data or gemini-virtualenv/../gemini/data
- Global:    /usr/local/share/gemini/gemini-config.yaml
- User only: $HOME/.gemini/gemini-config.yaml

Prefer installer based or global if you have system level permissions for
installation since it will work for all system users.
"""
import os
import yaml

CONFIG_FILE = "gemini-config.yaml"

def get_config_dirs():
    virtualenv_loc = __file__.find("gemini-virtualenv")
    anaconda_loc = __file__.find("anaconda")
    if anaconda_loc >= 0:
        base = __file__[:anaconda_loc]
        dirs = [os.path.join(base), os.path.join(base, "gemini")]
    elif virtualenv_loc >= 0:
        base = __file__[:virtualenv_loc]
        dirs = [os.path.join(base), os.path.join(base, "gemini")]
    else:
        dirs = []
    dirs.append("/usr/local/share/gemini")
    dirs.append(os.path.join(os.environ["HOME"], ".gemini"))
    return dirs

def _get_config_file(dirs=None):
    dirs = [] if dirs is None else dirs
    dnames = dirs + get_config_dirs()
    for dname in dnames:
        fname = os.path.join(dname, CONFIG_FILE)
        if os.path.exists(fname):
            return fname
    raise ValueError("GEMINI configuration file {0} not found in {1}.\n"
                     "Please ensure the GEMINI data is installed using the install-data.py script\n"
                     "http://gemini.readthedocs.org/en/latest/content/installation.html"
                     .format(CONFIG_FILE, dnames))

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
    dirs = [] if dirs is None else dirs
    dnames = dirs + get_config_dirs()
    for dname in dnames:
        if os.access(dname, os.W_OK) or \
                os.access(os.path.dirname(dname), os.W_OK):
            return os.path.join(dname, CONFIG_FILE)

    raise ValueError("Gemini configuration: "
                     "Could not find writeable directory: {0}".format(dnames))

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
