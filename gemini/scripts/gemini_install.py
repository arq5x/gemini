#!/usr/bin/env python
"""Installer for gemini: a lightweight db framework for disease and population genetics.

https://github.com/arq5x/gemini

Handles installation of:

- Required third party software
- Required Python libraries
- Gemini application
- Associated data files

Requires: Python 2.7 (or 2.6 and argparse), git, and compilers (gcc, g++)

Run gemini_install.py -h for usage.
"""
from __future__ import print_function
import argparse
import platform
import os
import shutil
import subprocess
import sys
import urllib2
import urllib

remotes = {"requirements_conda": "",
           "versioned_installations": "https://raw.githubusercontent.com/arq5x/gemini/master/versioning/",
           "anaconda": "http://repo.continuum.io/miniconda/Miniconda-latest-%s-x86%s.sh"}

remotes_dev = remotes.copy()
remotes_dev.update({
    "requirements_conda": "https://raw.githubusercontent.com/arq5x/gemini/dev/versioning/unstable/requirements_conda.txt",
})

remotes_bp = remotes_dev
remotes_bp.update({
    "requirements_conda": "https://raw.githubusercontent.com/brentp/gemini/dev/versioning/unstable/requirements_conda.txt",
})

def main(args, remotes=remotes):
    check_dependencies()
    clean_env_variables()
    work_dir = os.path.join(os.getcwd(), "tmpgemini_install")
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    os.chdir(work_dir)

    if args.gemini_version in ("unstable", "bp"):
        if args.gemini_version == "unstable":
            remotes = remotes_dev
        else:
            remotes = remotes_bp

        requirements_conda = remotes['requirements_conda']
        urllib.urlretrieve(requirements_conda, filename='_conda_dev.txt')

        # quick hack to support testing installs:
        if args.gemini_version == "bp":
            for f in '_conda_dev.txt':
                contents = open(f).read().replace('arq5x', 'brentp')
                with open(f, 'w') as fh:
                    fh.write(contents)

        remotes.update({'requirements_conda': '_conda_dev.txt'})

    elif args.gemini_version != 'latest':
        requirements_conda = os.path.join(remotes['versioned_installations'],
                                          args.gemini_version, 'requirements_conda.txt')
        try:
            urllib2.urlopen(requirements_conda)
        except:
            sys.exit('Gemini version %s could not be found. Try the latest version.' % args.gemini_version)
        remotes.update({'requirements_conda': requirements_conda})

    print("Installing isolated base python installation")
    make_dirs(args)
    anaconda = install_anaconda_python(args, remotes)
    print("Installing base gemini package...")
    gemini = install_conda_pkgs(anaconda, remotes, args)
    install_rest(gemini, args)
    print("Finished: gemini, tools and data installed")
    if args.tooldir:
        print(" Tools installed in:\n  %s" % args.tooldir)
        print(" NOTE: be sure to add %s/bin to your PATH." % args.tooldir)
    if args.install_data:
        print(" Data installed in:\n  %s" % args.datadir)
        print(" NOTE: Install data files for GERP_bp & CADD_scores (not installed by default).\n ")
    print(" Run tests with:\n  cd %s && bash master-test.sh" %
          (os.path.join(os.path.dirname(anaconda["dir"]), "github_gemini")))

    shutil.rmtree(work_dir)

def install_conda_pkgs(anaconda, remotes, args):
    if args.gemini_version != 'latest':
        req_file = '_conda-requirements-%s.txt' % args.gemini_version
        urllib.urlretrieve(remotes["requirements_conda"], filename=req_file)
        pkgs = ["--file", req_file]
    else:
        req_file = None
        pkgs = ["gemini"]
    channels = ["-c", "bioconda"]
    print(" ".join([anaconda["conda"], "install", "--yes"] + channels + pkgs))
    subprocess.check_call([anaconda["conda"], "install", "--yes"] + channels + pkgs)
    if req_file and os.path.exists(req_file):
        os.remove(req_file)
    return os.path.join(anaconda["dir"], "bin", "gemini")

def install_anaconda_python(args, remotes):
    """Provide isolated installation of Anaconda python.
    http://docs.continuum.io/anaconda/index.html
    """
    anaconda_dir = os.path.join(args.datadir, "anaconda")
    bindir = os.path.join(anaconda_dir, "bin")
    conda = os.path.join(bindir, "conda")
    if platform.mac_ver()[0]:
        distribution = "macosx"
    else:
        distribution = "linux"
    if platform.architecture()[0] == "32bit":
        arch = ""
    else:
        arch = "_64"
    if not os.path.exists(anaconda_dir) or not os.path.exists(conda):
        if os.path.exists(anaconda_dir):
            shutil.rmtree(anaconda_dir)
        url = remotes["anaconda"] % ("MacOSX" if distribution == "macosx" else "Linux", arch)
        if not os.path.exists(os.path.basename(url)):
            subprocess.check_call(["wget", "--continue", url])
        subprocess.check_call("bash %s -b -p %s" %
                              (os.path.basename(url), anaconda_dir), shell=True)
    return {"conda": conda,
            "dir": anaconda_dir}

def install_rest(gemini, args):
    """Install biological data and tests used by gemini.
    """
    if os.path.exists(os.path.join(args.datadir, "data")):
        annotation_dir = os.path.join(args.datadir, "data")
    else:
        annotation_dir = os.path.join(args.datadir, "gemini_data")
    cmd = [gemini, "--annotation-dir", annotation_dir, "update", "--dataonly"]
    if not args.install_data:
        cmd += ["--nodata"]
    if args.tooldir:
        cmd += ["--tooldir", args.tooldir]
    print(" ".join(cmd))
    subprocess.check_call(cmd)

def make_dirs(args):
    for dname in [args.datadir, args.tooldir]:
        if dname and not os.path.exists(dname):
            os.makedirs(dname)

def clean_env_variables():
    """Adjust environmental variables which can cause conflicts with installed anaconda python.
    """
    for k in ["PYTHONPATH", "PYTHONHOME"]:
        os.environ.pop(k, None)
    # https://docs.python.org/2/using/cmdline.html#envvar-PYTHONNOUSERSITE
    os.environ["PYTHONNOUSERSITE"] = "1"

def check_dependencies():
    """Ensure required tools for installation are present.
    """
    print("Checking required dependencies...")
    for cmd, url in [("wget", "http://www.gnu.org/software/wget/")]:
        try:
            retcode = subprocess.call([cmd, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except OSError:
            retcode = 127
        if retcode == 127:
            raise OSError("gemini requires %s (%s)" % (cmd, url))
        else:
            print(" %s found" % cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automated installer for gemini framework.")
    parser.add_argument("tooldir", help="Directory to install 3rd party software tools",
                        type=os.path.abspath)
    parser.add_argument("datadir", help="Directory to install anaconda python and gemini data files",
                        type=os.path.abspath)
    parser.add_argument("--nodata", help="Do not install data dependencies",
                        dest="install_data", action="store_false", default=True)
    parser.add_argument("--gemini-version", dest="gemini_version", default="latest",
                        help="Install one specific gemini version with a fixed dependency chain.")
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        main(parser.parse_args())
