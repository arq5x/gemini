#!/usr/bin/env python

import os
import sys

# gemini imports
import GeminiQuery


def query(parser, args):

    if (args.db is None):
        parser.print_help()

    if os.path.exists(args.db):

        gq = GeminiQuery.GeminiQuery(args.db)
        gq.run(args.query, args.gt_filter)

        if args.use_header:
            print gq.header

        for row in gq:
            print row

if __name__ == "__main__":
    main()
