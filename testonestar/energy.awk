# script to sum up all the feedbacks in onestar.out
#
# Run with "awk -f energy.awk < onestar.out"
#
BEGIN {
    total_fb0 = 0.0;
    total_fb1 = 0.0;
    total_fb2 = 0.0;
    }

# order is SNII, SNIa, Wind    
/feedback 0:/ { total_fb0 += $4; }
/feedback 1:/ { total_fb1 += $4; }
/feedback 2:/ { total_fb2 += $4; }

END {
    print "Total SNII Energy: ", total_fb0;
    print "Total SNIa Energy: ", total_fb1;
    print "Total Wind Energy: ", total_fb2;
    print "Total Energy: ", total_fb0 + total_fb1 + total_fb2;
}

