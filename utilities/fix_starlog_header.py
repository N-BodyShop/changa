#!/usr/bin/env python3
# Utility to fix broken starlog files output by ChaNGa prior to the fix for
# issue #294 (https://github.com/N-BodyShop/changa_uw/issues/294).
#
# Starlog files affected by this issue will have a header reporting 104
# byte records, but actually have records that are 96 bytes long.  This
# utility will rewrite your starlog file (and make a backup of the old one)
# The new file should "just work".

import struct
import sys
import os

def is_header_good(filename):
    starlog = open(filename, mode="rb")
    head = struct.unpack("!i", starlog.read(4))[0]
    print('starlog header gives a record size of', head)
    starlog.close()
    size = os.path.getsize(filename)
    if (size-4) % head != 0:
        print("File size does not match record size in header")
        return False
    else:
        return True

def fix_header(filename):
    size = os.path.getsize(filename)
    if (size-4) % 96 == 0:
        os.rename(filename, filename+".bak")
        starlog = open(filename+".bak", mode="rb")
        starlog.read(4) # Advance 4 bytes to discard header
        newstarlog = open(filename, mode="wb")
        print("Header should be 96 bytes")
        head = struct.pack("!i", 96)
        newstarlog.write(head)
        newstarlog.write(starlog.read())
        print("Wrote out fixed starlog @", filename)

    starlog.close()
    newstarlog.close()

if __name__ == "__main__":
    fname = sys.argv[1]
    nofix = is_header_good(fname)
    if nofix == False:
        fix_header(fname)
