"""Perform in-place updates of gemini and databases when installed into virtualenv.
"""
import os
import shutil
import subprocess
import sys

import gemini.config

def release(parser, args):
    """Update gemini to the latest release, along with associated data files.
    """
    url = "https://raw.github.com/arq5x/gemini/master/requirements.txt"
    repo = "https://github.com/arq5x/gemini.git"
    # update locally isolated python
    pip_bin = os.path.join(os.path.dirname(sys.executable), "pip")
    ei_bin = os.path.join(os.path.dirname(sys.executable), "easy_install")
    activate_bin = os.path.join(os.path.dirname(sys.executable), "activate")
    conda_bin = os.path.join(os.path.dirname(sys.executable), "conda")
    # Work around issue with distribute where asks for 'distribute==0.0'
    try:
        subprocess.check_call([ei_bin, "--upgrade", "distribute"])
    except subprocess.CalledProcessError:
        subprocess.check_call([pip_bin, "install", "--upgrade", "distribute"])
    if os.path.exists(conda_bin):
        pkgs = ["cython", "ipython", "jinja2", "nose", "numpy",
                "pip", "pycrypto", "pyparsing", "pysam", "pyyaml",
                "pyzmq", "pandas", "scipy"]
        subprocess.check_call([conda_bin, "install", "--yes"] + pkgs)
    elif os.path.exists(activate_bin):
        pass
    else:
        raise NotImplementedError("Can only upgrade gemini installed in anaconda or virtualenv")
    # update libraries
    subprocess.check_call([pip_bin, "install", "-r", url])
    if args.devel:
        print("Installing latest GEMINI development version")
        subprocess.check_call([pip_bin, "install", "--upgrade", "--no-deps",
                               "git+%s" % repo])
    # update datafiles
    config = gemini.config.read_gemini_config()
    subprocess.check_call([sys.executable, _get_install_script(), config["annotation_dir"]])
    print "Gemini upgraded to latest version"
    # update tests
    test_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(pip_bin))),
                            "gemini")
    _update_testbase(test_dir, repo)
    print "Run test suite with: cd %s && bash master-test.sh" % test_dir

def _get_install_script():
    try:
        import pkg_resources
        return pkg_resources.resource_filename(__name__, "install-data.py")
    except ImportError:
        return os.path.join(os.path.dirname(__file__), "install-data.py")

def _update_testbase(repo_dir, repo):
    cur_dir = os.getcwd()
    needs_git = True
    if os.path.exists(repo_dir):
        os.chdir(repo_dir)
        try:
            subprocess.check_call(["git", "pull", "origin", "master"])
            needs_git = False
        except:
            os.chdir(cur_dir)
            shutil.rmtree(repo_dir)
    if needs_git:
        os.chdir(os.path.split(repo_dir)[0])
        subprocess.check_call(["git", "clone", repo])
    os.chdir(cur_dir)
