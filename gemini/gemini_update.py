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
    base = os.path.dirname(os.path.realpath(sys.executable))
    gemini_cmd = os.path.join(base, "gemini")
    pip_bin = os.path.join(base, "pip")
    ei_bin = os.path.join(base, "easy_install")
    activate_bin = os.path.join(base, "activate")
    conda_bin = os.path.join(base, "conda")
    print conda_bin
    if not args.dataonly:
        # Work around issue with distribute where asks for 'distribute==0.0'
        # try:
        #     subprocess.check_call([ei_bin, "--upgrade", "distribute"])
        # except subprocess.CalledProcessError:
        #     try:
        #         subprocess.check_call([pip_bin, "install", "--upgrade", "distribute"])
        #     except subprocess.CalledProcessError:
        #         pass
        if os.path.exists(conda_bin):
            pkgs = ["bx-python", "conda", "cython", "ipython", "jinja2", "nose", "numpy",
                    "pip", "pycrypto", "pyparsing", "pysam", "pyyaml",
                    "pyzmq", "pandas", "scipy"]
            channels = ["-c", "https://conda.binstar.org/bcbio"]
            subprocess.check_call([conda_bin, "install", "--yes", "numpy"])
            subprocess.check_call([conda_bin, "install", "--yes"] + channels + pkgs)
        elif os.path.exists(activate_bin):
            pass
        else:
            raise NotImplementedError("Can only upgrade gemini installed in anaconda or virtualenv")
        # allow downloads excluded in recent pip (1.5 or greater) versions
        try:
            p = subprocess.Popen([pip_bin, "--version"], stdout=subprocess.PIPE)
            pip_version = p.communicate()[0].split()[1]
        except:
            pip_version = ""
        pip_compat = []
        if pip_version >= "1.5":
            for req in ["python-graph-core", "python-graph-dot"]:
                pip_compat += ["--allow-external", req, "--allow-unverified", req]
        # update libraries
        subprocess.check_call([pip_bin, "install"] + pip_compat + ["-r", url])
        if args.devel:
            print("Installing latest GEMINI development version")
            subprocess.check_call([pip_bin, "install", "--upgrade", "--no-deps",
                                   "git+%s" % repo])
        print "Gemini upgraded to latest version"
    # update datafiles
    config = gemini.config.read_gemini_config( args = args )
    extra_args = ["--extra=%s" % x for x in args.extra]
    subprocess.check_call([sys.executable, _get_install_script(), config["annotation_dir"]] + extra_args)
    print "Gemini data files updated"
    # update tests
    if not args.dataonly:
        test_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(pip_bin))),
                                "gemini")
        if not os.path.exists(test_dir) or os.path.isdir(test_dir):
            _update_testbase(test_dir, repo, gemini_cmd)
            print "Run test suite with: cd %s && bash master-test.sh" % test_dir

def _get_install_script():
    try:
        import pkg_resources
        return pkg_resources.resource_filename(__name__, "install-data.py")
    except ImportError:
        return os.path.join(os.path.dirname(__file__), "install-data.py")

def _update_testbase(repo_dir, repo, gemini_cmd):
    cur_dir = os.getcwd()
    needs_git = True
    if os.path.exists(repo_dir):
        os.chdir(repo_dir)
        try:
            subprocess.check_call(["git", "pull", "origin", "master", "--tags"])
            needs_git = False
        except:
            os.chdir(cur_dir)
            shutil.rmtree(repo_dir)
    if needs_git:
        os.chdir(os.path.split(repo_dir)[0])
        subprocess.check_call(["git", "clone", repo])
    os.chdir(repo_dir)
    _update_testdir_revision(gemini_cmd)
    os.chdir(cur_dir)

def _update_testdir_revision(gemini_cmd):
    """Update test directory to be in sync with a tagged installed version or development.
    """
    try:
        p = subprocess.Popen([gemini_cmd, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        gversion = p.communicate()[0].split()[1]
    except:
        gversion = ""
    tag = ""
    if gversion:
        try:
            p = subprocess.Popen("git tag -l | grep %s" % gversion, stdout=subprocess.PIPE, shell=True)
            tag = p.communicate()[0].strip()
        except:
            tag = ""
    if tag:
        subprocess.check_call(["git", "checkout", "tags/%s" % tag])
        pass
    else:
        subprocess.check_call(["git", "reset", "--hard", "HEAD"])
