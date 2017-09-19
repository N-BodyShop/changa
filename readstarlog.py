# Utility to peak at the starlog file.
# The file is written in "network order"
# Based on the size in the header, guesses are made about the number
# of variables in each record, and the first record is printed out
# for confirmation.

import struct
import sys

fmthead = "!i"
fmtevent      = "!qqdddddddddd"
fmteventOLD   = "!IIdddddddddd"
fmteventH2    = "!qqddddddddddd"
fmteventH2OLD = "!IIddddddddddd"

filelog = open(sys.argv[1], mode="rb")
head = filelog.read(4)
size = struct.unpack(fmthead, head)
size = size[0]  # a tuple is returned above

print 'starlog record size is ', size

if size == struct.calcsize(fmteventOLD) :
    print '32 bit iOrders format'
    fmt = fmteventOLD
elif size == struct.calcsize(fmtevent) :
    # test read of first few bytes
    tstbuf = filelog.read(16)
    first = struct.unpack('qq', tstbuf)
    if first[0] < first[1] or first[0] < 0 or first[1] < 0:
        print '32 bit iOrders with H2 or BH'
        fmt = fmteventH2OLD
        filelog.seek(4)
    else :
        print '64 bit iOrders format'
        fmt = fmtevent
elif size == struct.calcsize(fmteventH2) :
    print '64 bit iOrders format with H2 or BH'
    fmt = fmteventH2

rec1 = filelog.read(struct.calcsize(fmt))
print struct.unpack(fmt, rec1)
