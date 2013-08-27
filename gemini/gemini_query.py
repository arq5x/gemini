#!/usr/bin/env python

import os

# gemini imports
import GeminiQuery


def query(parser, args):

    if (args.db is None):
        parser.print_help()

    if os.path.exists(args.db):

        gq = GeminiQuery.GeminiQuery(args.db)
        gq.run(args.query, args.gt_filter, 
               args.show_variant_samples,
               args.sample_delim)

        if args.use_header and not args.use_json:
            print gq.header

        for row in gq:
            if not args.use_json:
                print row
            else:
                print row.get_json()

if __name__ == "__main__":
    main()
