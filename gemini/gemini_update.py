"""Perform in-place updates of gemini and databases when installed into virtualenv.
"""
import os
import subprocess
import sys

import gemini.config

def release(parser, args):
    """Update gemini to the latest release, along with associated data files.
    """
    url = "https://raw.github.com/arq5x/gemini/master/requirements.txt"
    # update locally isolated python
    pip_bin = os.path.join(os.path.dirname(sys.executable), "pip")
    activate_bin = os.path.join(os.path.dirname(sys.executable), "activate")
    conda_bin = os.path.join(os.path.dirname(sys.executable), "conda")
    if os.path.exists(conda_bin):
        pkgs = ["cython", "distribute", "ipython", "nose", "numpy",
                "pip", "pycrypto", "pyparsing", "pysam", "pyyaml",
                "pyzmq", "pandas", "scipy"]
        subprocess.check_call([conda_bin, "install", "--yes"] + pkgs)
    elif os.path.exists(activate_bin):
        subprocess.check_call([pip_bin, "install", "--upgrade", "distribute"])
    else:
        raise NotImplementedError("Can only upgrade gemini installed in anaconda or virtualenv")
    # update libraries
    subprocess.check_call([pip_bin, "install", "-r", url])
    # update datafiles
    config = gemini.config.read_gemini_config()
    install_script = os.path.join(os.path.dirname(__file__), "install-data.py")
    subprocess.check_call([sys.executable, install_script, config["annotation_dir"]])
    print "Gemini upgraded to latest version"
    # update tests
    test_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(pip_bin))),
                            "gemini")
    if os.path.exists(test_dir) and os.path.exists(os.path.join(test_dir, "master-test.sh")):
        os.chdir(test_dir)
        subprocess.check_call(["git", "pull", "origin", "master"])
        print "Run test suite with: cd %s && bash master-test.sh" % test_dir
