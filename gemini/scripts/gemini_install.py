#!/usr/bin/env python
"""Installer for gemini: a lightweight db framework for disease and population genetics.

https://github.com/arq5x/gemini

Handles installation of:

- Required third party software
- Required Python libraries
- Gemini application
- Associated data files

Requires: Python 2.7, git

Run gemini_install.py -h for usage.
"""
import argparse
import os
import shutil
import subprocess
import sys

remotes = {"requirements":
           "https://raw.github.com/arq5x/gemini/master/requirements.txt",
           "cloudbiolinux":
           "https://github.com/chapmanb/cloudbiolinux.git",
           "virtualenv":
           "https://raw.github.com/pypa/virtualenv/master/virtualenv.py",
           "gemini":
           "https://github.com/arq5x/gemini.git"}

def main(args):
    check_dependencies()
    work_dir = os.path.join(os.getcwd(), "tmpgemini_install")
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    os.chdir(work_dir)
    print "Installing gemini..."
    make_dirs(args)
    gemini = install_gemini(remotes, args.datadir, args.tooldir, args.sudo)
    cbl = get_cloudbiolinux(remotes["cloudbiolinux"])
    fabricrc = write_fabricrc(cbl["fabricrc"], args.tooldir, args.datadir,
                              "ubuntu", args.sudo)
    if args.install_tools:
        print "Installing associated tools..."
        install_tools(gemini["fab"], cbl["tool_fabfile"], fabricrc)
    if args.install_data:
        print "Installing gemini data..."
        install_data(gemini["python"], gemini["data_script"], args.datadir)
    test_script = install_testbase(args.datadir, remotes["gemini"])
    print "Finished: gemini, tools and data installed"
    print " Tools installed in:\n  %s" % args.tooldir
    print " Data installed in:\n  %s" % args.datadir
    print " Run tests with:\n  cd %s && bash %s" % (os.path.dirname(test_script),
                                                    os.path.basename(test_script))
    shutil.rmtree(work_dir)

def install_gemini(remotes, datadir, tooldir, use_sudo):
    """Install a virtualenv containing gemini plus dependencies.
    """
    virtualenv_dir = os.path.join(datadir, "gemini-virtualenv")
    if not os.path.exists(virtualenv_dir):
        subprocess.check_call(["wget", remotes["virtualenv"]])
        subprocess.check_call([sys.executable, "virtualenv.py", "--system-site-packages",
                               "--distribute", virtualenv_dir])
        os.remove("virtualenv.py")
    pip_cmd = os.path.join(virtualenv_dir, "bin", "pip")
    ez_cmd = os.path.join(virtualenv_dir, "bin", "easy_install")
    # work around issue with latest version of pip on MacOSX: https://github.com/pypa/pip/issues/829
    subprocess.check_call([ez_cmd, "pip==1.2.1"])
    subprocess.check_call([pip_cmd, "install", "--upgrade", "fabric"])
    subprocess.check_call([pip_cmd, "install", "--upgrade", "distribute"])
    subprocess.check_call([pip_cmd, "install", "--upgrade", "cython"])
    subprocess.check_call([pip_cmd, "install", "--upgrade", "pyyaml"])
    subprocess.check_call([pip_cmd, "install", "-r", remotes["requirements"]])
    for script in ["gemini"]:
        final_script = os.path.join(tooldir, "bin", script)
        ve_script = os.path.join(virtualenv_dir, "bin", script)
        if not os.path.exists(final_script):
            sudo_cmd = ["sudo"] if use_sudo else []
            subprocess.check_call(sudo_cmd + ["mkdir", "-p", os.path.dirname(final_script)])
            cmd = ["ln", "-s", ve_script, final_script]
            subprocess.check_call(sudo_cmd + cmd)
    python_bin = os.path.join(virtualenv_dir, "bin", "python")
    library_loc = subprocess.check_output("%s -c 'import gemini; print gemini.__file__'" % python_bin,
                                          shell=True)
    return {"fab": os.path.join(virtualenv_dir, "bin", "fab"),
            "data_script": os.path.join(os.path.dirname(library_loc.strip()), "install-data.py"),
            "python": python_bin}

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

def install_data(python_cmd, data_script, datadir):
    """Install biological data used by gemini.
    """
    subprocess.check_call([python_cmd, data_script, datadir])

def install_testbase(datadir, repo):
    """Clone or update gemini code so we have the latest test suite.
    """
    gemini_dir = os.path.join(datadir, "gemini")
    cur_dir = os.getcwd()
    if not os.path.exists(gemini_dir):
        os.chdir(os.path.split(gemini_dir)[0])
        subprocess.check_call(["git", "clone", repo])
    else:
        os.chdir(gemini_dir)
        subprocess.check_call(["git", "pull", "origin", "master"])
    os.chdir(cur_dir)
    return os.path.join(gemini_dir, "master-test.sh")

def write_fabricrc(base_file, tooldir, datadir, distribution, use_sudo):
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
                elif line.startswith("distribution"):
                    line = "distribution = %s\n" % distribution
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
    print "Checking required dependencies"
    try:
        subprocess.check_call(["git", "--version"])
    except OSError:
        raise OSError("bcbio-nextgen installer requires Git (http://git-scm.com/)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automated installer for gemini framework.")
    parser.add_argument("tooldir", help="Directory to install 3rd party software tools",
                        type=os.path.abspath)
    parser.add_argument("datadir", help="Directory to install gemini data files",
                        type=os.path.abspath)
    parser.add_argument("--nosudo", help="Specify we cannot use sudo for commands",
                        dest="sudo", action="store_false", default=True)
    parser.add_argument("--notools", help="Do not install tool dependencies",
                        dest="install_tools", action="store_false", default=True)
    parser.add_argument("--nodata", help="Do not install data dependencies",
                        dest="install_data", action="store_false", default=True)
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        main(parser.parse_args())
