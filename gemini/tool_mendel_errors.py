#!/usr/bin/env python
import os
from gemini_inheritance_model_utils import GeminiInheritanceModelFactory

def run(parser, args):
    if os.path.exists(args.db):
        factory = GeminiInheritanceModelFactory(args, model="mendel_violations")
        factory.get_candidates()

