#!/usr/bin/env python
import os
import sys
from gemini_inheritance_model_utils import GeminiInheritanceModelFactory


def run(parser, args):
    if os.path.exists(args.db):
        de_novo_factory = \
            GeminiInheritanceModelFactory(args, model="de_novo")
        de_novo_factory.get_candidates()

