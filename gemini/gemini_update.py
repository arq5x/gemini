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
    repo = "https://github.com/arq5x/gemini.git"
    # update locally isolated python
    base = os.path.dirname(os.path.realpath(sys.executable))
    gemini_cmd = os.path.join(base, "gemini")
    pip_bin = os.path.join(base, "pip")
    conda_bin = os.path.join(base, "conda")
    if not args.dataonly:
        if not os.path.exists(conda_bin):
            raise NotImplementedError("Can only upgrade gemini installed with anaconda")
        subprocess.check_call([conda_bin, "install", "-y", "-c", "bioconda", "gemini", "pip"])
        if args.devel:
            subprocess.check_call([pip_bin, "install", "--upgrade", "--no-deps",
                                   "git+%s" % repo])
        print "Gemini upgraded to latest version"
    if args.tooldir:
        print "Upgrading associated tools..."
        if os.path.exists(conda_bin):
            anaconda_dir = os.path.dirname(os.path.dirname(conda_bin))
            link_tools(args.tooldir, anaconda_dir)
    # update datafiles
    config = gemini.config.read_gemini_config(args=args, allow_missing=True)
    gemini.config.write_gemini_config(config)
    if args.install_data:
        extra_args = ["--extra=%s" % x for x in args.extra]
        subprocess.check_call([sys.executable, _get_install_script(), config["annotation_dir"]] + extra_args)
        print "Gemini data files updated"
    # update tests
    if not args.dataonly:
        test_dir = os.path.join(os.path.dirname(os.path.dirname(base)), "github_gemini")
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
        print "cloning %s to %s" % (repo, repo_dir)
        subprocess.check_call(["git", "clone", repo, repo_dir])
    os.chdir(repo_dir)
    try:
        _update_testdir_revision(gemini_cmd)
    except Exception, msg:
        print("Unable to update to revision, skipping: %s" % msg)
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
            p = subprocess.Popen("git tag -l | grep %s$" % gversion, stdout=subprocess.PIPE, shell=True)
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
    bin_dir = os.path.join(tooldir, "bin")
    if not os.path.exists(bin_dir):
        os.makedirs(bin_dir)
    bins = [("", "grabix"), ("", "gemini"),
            ("gemini_", "python"), ("gemini_", "conda"), ("gemini_", "pip")]
    for prefix, binfile in bins:
        orig_file = os.path.join(anaconda_dir, "bin", binfile)
        final_file = os.path.join(bin_dir, "%s%s" % (prefix, binfile))
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
