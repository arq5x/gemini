#!/usr/bin/env python
"""
Store the path for GEMINI data-dir
"""

from gemini.config import read_gemini_config

config = read_gemini_config()
anno_dirname = config["annotation_dir"]

print anno_dirname



