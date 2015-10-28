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
    fab_cmd = os.path.join(base, "fab")
    activate_bin = os.path.join(base, "activate")
    conda_bin = os.path.join(base, "conda")
    if not args.dataonly:
        if os.path.exists(conda_bin):
            clean_env_variables()
            pkgs = ["bcolz", "bx-python", "conda", "cyordereddict", "cython", "grabix", "ipyparallel",
                    "ipython-cluster-helper",
                    "jinja2", "nose", "numexpr", "numpy", "openssl", "pip", "pybedtools",
                    "pycrypto", "pyparsing", "python-graph-core", "python-graph-dot",
                    "pysam", "pyyaml", "pyzmq", "pandas", "scipy"]
            channels = ["-c", "bcbio", "-c", "bioconda"]
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
        # Set PIP SSL certificate to installed conda certificate to avoid SSL errors
        if os.path.exists(conda_bin):
            anaconda_dir = os.path.dirname(os.path.dirname(conda_bin))
            cert_file = os.path.join(anaconda_dir, "ssl", "cert.pem")
            if os.path.exists(cert_file):
                os.environ["PIP_CERT"] = cert_file
        subprocess.check_call([pip_bin, "install"] + pip_compat + ["-r", url])
        if args.devel:
            print("Installing latest GEMINI development version")
            subprocess.check_call([pip_bin, "install", "--upgrade", "--no-deps",
                                   "git+%s" % repo])
        print "Gemini upgraded to latest version"
    if args.tooldir:
        print "Upgrading associated tools..."
        if os.path.exists(conda_bin):
            anaconda_dir = os.path.dirname(os.path.dirname(conda_bin))
            link_tools(args.tooldir, anaconda_dir)
    # update datafiles
    config = gemini.config.read_gemini_config(args=args)
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

def clean_env_variables():
    """Adjust environmental variables which can cause conflicts with installed anaconda python.
    """
    for k in ["PYTHONPATH", "PYTHONHOME"]:
        os.environ.pop(k, None)
    # https://docs.python.org/2/using/cmdline.html#envvar-PYTHONNOUSERSITE
    os.environ["PYTHONNOUSERSITE"] = "1"

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
        p = os.path.split(repo_dir)[0]
        print "cloning %s to %s" % (repo, p)
        os.chdir(p)
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

def link_tools(tooldir, anaconda_dir):
    """Link tools installed via conda into the tool directory.
    """
    bins = ["grabix"]
    for binfile in bins:
        orig_file = os.path.join(anaconda_dir, "bin", binfile)
        final_file = os.path.join(tooldir, "bin", binfile)
        if os.path.exists(orig_file):
            _do_link(orig_file, final_file)

def _do_link(orig_file, final_file):
    """Perform a soft link of the original file into the final location.

    We need the symlink to point to /anaconda/bin directory, not the real location
    in the pkgs directory so conda can resolve LD_LIBRARY_PATH and the interpreters.
    """
    needs_link = True
    # working symlink, check if already in the right place or remove it
    if os.path.exists(final_file):
        if (os.path.realpath(final_file) == os.path.realpath(orig_file) and
              orig_file == os.path.normpath(os.path.join(os.path.dirname(final_file), os.readlink(final_file)))):
            needs_link = False
        else:
            os.remove(final_file)
    # broken symlink
    elif os.path.lexists(final_file):
        os.unlink(final_file)
    if needs_link:
        os.symlink(os.path.relpath(orig_file, os.path.dirname(final_file)), final_file)
