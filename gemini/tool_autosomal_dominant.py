#!/usr/bin/env python
import os
import sys
from gemini_inheritance_model_utils import GeminiInheritanceModelFactory

def run(parser, args):
    if os.path.exists(args.db):
        auto_dominant_factory = \
            GeminiInheritanceModelFactory(args, model="auto_dom")
        auto_dominant_factory.get_candidates()

