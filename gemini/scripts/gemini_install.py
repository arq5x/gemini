#!/usr/bin/env python
"""Installer for gemini: a lightweight db framework for disease and population genetics.

https://github.com/arq5x/gemini

Handles installation of:

- Required third party software
- Required Python libraries
- Gemini application
- Associated data files

Requires: Python 2.7, git, and compilers (gcc, g++)

Run gemini_install.py -h for usage.
"""
import argparse
import platform
import os
import shutil
import subprocess
import sys
import urllib2

remotes = {"requirements_pip":
           "https://raw.github.com/arq5x/gemini/master/requirements.txt",
            "requirements_conda":
           "",
            "versioned_installations":
            "https://raw.githubusercontent.com/arq5x/gemini/master/versioning/",
           "cloudbiolinux":
           "https://github.com/chapmanb/cloudbiolinux.git",
           "gemini":
           "https://github.com/arq5x/gemini.git",
           "anaconda":
           "http://repo.continuum.io/miniconda/Miniconda-3.5.5-%s-x86_64.sh"}

def main(args):
    check_dependencies()
    work_dir = os.path.join(os.getcwd(), "tmpgemini_install")
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    os.chdir(work_dir)

    if args.gemini_version != 'latest':
        requirements_pip = os.path.join(remotes['versioned_installations'],
                                        args.gemini_version, 'requirements_pip.txt')
        requirements_conda = os.path.join(remotes['versioned_installations'],
                                          args.gemini_version, 'requirements_conda.txt')
        try:
            urllib2.urlopen(requirements_pip)
        except:
            sys.exit('Gemini version %s could not be found. Try the latest version.' % args.gemini_version)
        remotes.update({'requirements_pip': requirements_pip, 'requirements_conda': requirements_conda})

    print "Installing isolated base python installation"
    make_dirs(args)
    anaconda = install_anaconda_python(args, remotes)
    print "Installing gemini..."
    install_conda_pkgs(anaconda, remotes, args)
    gemini = install_gemini(anaconda, remotes, args.datadir, args.tooldir, args.sudo)
    if args.install_tools:
        cbl = get_cloudbiolinux(remotes["cloudbiolinux"])
        fabricrc = write_fabricrc(cbl["fabricrc"], args.tooldir, args.datadir, args.sudo)
        print "Installing associated tools..."
        install_tools(gemini["fab"], cbl["tool_fabfile"], fabricrc)
    os.chdir(work_dir)
    install_data(gemini["python"], gemini["data_script"], args)
    os.chdir(work_dir)
    test_script = install_testbase(args.datadir, remotes["gemini"], gemini)
    print "Finished: gemini, tools and data installed"
    print " Tools installed in:\n  %s" % args.tooldir
    print " Data installed in:\n  %s" % args.datadir
    print " Run tests with:\n  cd %s && bash %s" % (os.path.dirname(test_script),
                                                    os.path.basename(test_script))
    print " NOTE: be sure to add %s/bin to your PATH." % args.tooldir
    print " NOTE: Install data files for GERP_bp & CADD_scores (not installed by default).\n "

    shutil.rmtree(work_dir)

def install_gemini(anaconda, remotes, datadir, tooldir, use_sudo):
    """Install gemini plus python dependencies inside isolated Anaconda environment.
    """
    # Work around issue with distribute where asks for 'distribute==0.0'
    # try:
    #     subprocess.check_call([anaconda["easy_install"], "--upgrade", "distribute"])
    # except subprocess.CalledProcessError:
    #     try:
    #         subprocess.check_call([anaconda["pip"], "install", "--upgrade", "distribute"])
    #     except subprocess.CalledProcessError:
    #         pass
    # Ensure latest version of fabric for running CloudBioLinux
    subprocess.check_call([anaconda["pip"], "install", "fabric>=1.7.0"])
    # allow downloads excluded in recent pip (1.5 or greater) versions
    try:
        p = subprocess.Popen([anaconda["pip"], "--version"], stdout=subprocess.PIPE)
        pip_version = p.communicate()[0].split()[1]
    except:
        pip_version = ""
    pip_compat = []
    if pip_version >= "1.5":
        for req in ["python-graph-core", "python-graph-dot"]:
            pip_compat += ["--allow-external", req, "--allow-unverified", req]
    subprocess.check_call([anaconda["pip"], "install"] + pip_compat + ["-r", remotes["requirements_pip"]])
    python_bin = os.path.join(anaconda["dir"], "bin", "python")
    _cleanup_problem_files(anaconda["dir"])
    _add_missing_inits(python_bin)
    for final_name, ve_name in [("gemini", "gemini"), ("gemini_python", "python"),
                                ("gemini_pip", "pip")]:
        final_script = os.path.join(tooldir, "bin", final_name)
        ve_script = os.path.join(anaconda["dir"], "bin", ve_name)
        sudo_cmd = ["sudo"] if use_sudo else []
        if os.path.lexists(final_script):
            subprocess.check_call(sudo_cmd + ["rm", "-f", final_script])
        else:
            subprocess.check_call(sudo_cmd + ["mkdir", "-p", os.path.dirname(final_script)])
        cmd = ["ln", "-s", ve_script, final_script]
        subprocess.check_call(sudo_cmd + cmd)
    library_loc = subprocess.check_output("%s -c 'import gemini; print gemini.__file__'" % python_bin,
                                          shell=True)
    return {"fab": os.path.join(anaconda["dir"], "bin", "fab"),
            "data_script": os.path.join(os.path.dirname(library_loc.strip()), "install-data.py"),
            "python": python_bin,
            "cmd": os.path.join(anaconda["dir"], "bin", "gemini")}

def install_conda_pkgs(anaconda, remotes, args):
    if args.gemini_version != 'latest':
        pkgs = ["--file", remotes['requirements_conda']]
    else:
        pkgs = ["bx-python", "conda", "cython", "ipython", "jinja2", "nose", "numpy",
            "pip", "pycrypto", "pyparsing", "pysam", "pyyaml",
            "pyzmq", "pandas", "scipy"]
    channels = ["-c", "https://conda.binstar.org/bcbio"]
    subprocess.check_call([anaconda["conda"], "install", "--yes"] + channels + pkgs)

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
    if not os.path.exists(anaconda_dir) or not os.path.exists(conda):
        if os.path.exists(anaconda_dir):
            shutil.rmtree(anaconda_dir)
        url = remotes["anaconda"] % ("MacOSX" if distribution == "macosx" else "Linux")
        if not os.path.exists(os.path.basename(url)):
            subprocess.check_call(["wget", url])
        subprocess.check_call("bash %s -b -p %s" %
                              (os.path.basename(url), anaconda_dir), shell=True)
    return {"conda": conda,
            "pip": os.path.join(bindir, "pip"),
            "easy_install": os.path.join(bindir, "easy_install"),
            "dir": anaconda_dir}

def _add_missing_inits(python_bin):
    """pip/setuptools strips __init__.py files with namespace declarations.
    I have no idea why, but this adds them back.
    """
    library_loc = subprocess.check_output("%s -c 'import pygraph.classes.graph; "
                                          "print pygraph.classes.graph.__file__'" % python_bin,
                                          shell=True)
    pygraph_init = os.path.normpath(os.path.join(os.path.dirname(library_loc.strip()), os.pardir,
                                                 "__init__.py"))
    if not os.path.exists(pygraph_init):
        with open(pygraph_init, "w") as out_handle:
            out_handle.write("__import__('pkg_resources').declare_namespace(__name__)\n")

def _cleanup_problem_files(venv_dir):
    """Remove problem bottle items in PATH which conflict with site-packages
    """
    for cmd in ["bottle.py", "bottle.pyc"]:
        bin_cmd = os.path.join(venv_dir, "bin", cmd)
        if os.path.exists(bin_cmd):
            os.remove(bin_cmd)

def install_tools(fab_cmd, fabfile, fabricrc):
    """Install 3rd party tools used by Gemini using a custom CloudBioLinux flavor.
    """
    tools = ["tabix", "grabix", "samtools", "bedtools"]
    flavor_dir = os.path.join(os.getcwd(), "gemini-flavor")
    if not os.path.exists(flavor_dir):
        os.makedirs(flavor_dir)
    with open(os.path.join(flavor_dir, "main.yaml"), "w") as out_handle:
        out_handle.write("packages:\n")
        out_handle.write("  - bio_nextgen\n")
        out_handle.write("libraries:\n")
    with open(os.path.join(flavor_dir, "custom.yaml"), "w") as out_handle:
        out_handle.write("bio_nextgen:\n")
        for tool in tools:
            out_handle.write("  - %s\n" % tool)
    cmd = [fab_cmd, "-f", fabfile, "-H", "localhost", "-c", fabricrc,
           "install_biolinux:target=custom,flavor=%s" % flavor_dir]
    subprocess.check_call(cmd)

def install_data(python_cmd, data_script, args):
    """Install biological data used by gemini.
    """
    data_dir = os.path.join(args.datadir, "gemini_data") if args.sharedpy else args.datadir
    cmd = [python_cmd, data_script, data_dir]
    if args.install_data:
        print "Installing gemini data..."
    else:
        cmd.append("--nodata")
    subprocess.check_call(cmd)

def install_testbase(datadir, repo, gemini):
    """Clone or update gemini code so we have the latest test suite.
    """
    gemini_dir = os.path.join(datadir, "gemini")
    cur_dir = os.getcwd()
    needs_git = True
    if os.path.exists(gemini_dir):
        os.chdir(gemini_dir)
        try:
            subprocess.check_call(["git", "pull", "origin", "master", "--tags"])
            needs_git = False
        except:
            os.chdir(cur_dir)
            shutil.rmtree(gemini_dir)
    if needs_git:
        os.chdir(os.path.split(gemini_dir)[0])
        subprocess.check_call(["git", "clone", repo])
    os.chdir(gemini_dir)
    _update_testdir_revision(gemini["cmd"])
    os.chdir(cur_dir)
    return os.path.join(gemini_dir, "master-test.sh")

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

def write_fabricrc(base_file, tooldir, datadir, use_sudo):
    out_file = os.path.join(os.getcwd(), os.path.basename(base_file))
    with open(base_file) as in_handle:
        with open(out_file, "w") as out_handle:
            for line in in_handle:
                if line.startswith("system_install"):
                    line = "system_install = %s\n" % tooldir
                elif line.startswith("local_install"):
                    line = "local_install = %s/install\n" % tooldir
                elif line.startswith("data_files"):
                    line = "data_files = %s\n" % datadir
                elif line.startswith("use_sudo"):
                    line = "use_sudo = %s\n" % use_sudo
                elif line.startswith("edition"):
                    line = "edition = minimal\n"
                elif line.startswith("#galaxy_home"):
                    line = "galaxy_home = %s\n" % os.path.join(datadir, "galaxy")
                out_handle.write(line)
    return out_file

def make_dirs(args):
    sudo_cmd = ["sudo"] if args.sudo else []
    for dname in [args.datadir, args.tooldir]:
        if not os.path.exists(dname):
            subprocess.check_call(sudo_cmd + ["mkdir", "-p", dname])
            username = subprocess.check_output("echo $USER", shell=True).strip()
            subprocess.check_call(sudo_cmd + ["chown", username, dname])

def get_cloudbiolinux(repo):
    base_dir = os.path.join(os.getcwd(), "cloudbiolinux")
    if not os.path.exists(base_dir):
        subprocess.check_call(["git", "clone", repo])
    return {"fabricrc": os.path.join(base_dir, "config", "fabricrc.txt"),
            "tool_fabfile": os.path.join(base_dir, "fabfile.py")}

def check_dependencies():
    """Ensure required tools for installation are present.
    """
    print "Checking required dependencies..."
    for cmd, url in [("git", "http://git-scm.com/"),
                     ("wget", "http://www.gnu.org/software/wget/"),
                     ("curl", "http://curl.haxx.se/")]:
        try:
            retcode = subprocess.call([cmd, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except OSError:
            retcode = 127
        if retcode == 127:
            raise OSError("gemini requires %s (%s)" % (cmd, url))
        else:
            print " %s found" % cmd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automated installer for gemini framework.")
    parser.add_argument("tooldir", help="Directory to install 3rd party software tools",
                        type=os.path.abspath)
    parser.add_argument("datadir", help="Directory to install gemini data files",
                        type=os.path.abspath)
    parser.add_argument("--gemini-version", dest="gemini_version", default="latest",
                        help="Install one specific gemini version with a fixed dependency chain.")
    parser.add_argument("--nosudo", help="Specify we cannot use sudo for commands",
                        dest="sudo", action="store_false", default=True)
    parser.add_argument("--notools", help="Do not install tool dependencies",
                        dest="install_tools", action="store_false", default=True)
    parser.add_argument("--nodata", help="Do not install data dependencies",
                        dest="install_data", action="store_false", default=True)
    parser.add_argument("--sharedpy", help=("Indicate we share an Anaconda Python directory with "
                                            "another project. Creates unique gemini data directory."),
                        action="store_true", default=False)
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        main(parser.parse_args())
