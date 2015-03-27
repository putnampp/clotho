#!/usr/bin/env python
#

import sys, json

from optparse import OptionParser

parser = OptionParser()

parser.add_option( "-i", "--input", dest="in_file", help="Input JSON file", default="")
parser.add_option( "-o", "--output", dest="out_file", help="Output JSON file", default="")

(options,args)=parser.parse_args()

if len(options.in_file) == 0:
    sys.exit(1)

with open( options.in_file, 'r') as rfile:
    data = json.load(rfile)

if len(options.out_file) == 0:
    json.dumps( data, sort_keys=True, indent=2)
else:
    with open( options.out_file, 'w') as lfile:
        res = json.dump( data, lfile, sort_keys=True, indent=2)
