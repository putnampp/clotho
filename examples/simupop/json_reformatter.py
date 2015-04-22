#   Copyright 2015 Patrick Putnam
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
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
