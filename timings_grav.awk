# Awk script to get time profile of multistep run
# This version is for gravity only runs
# Run this script with:
# awk -f timings_grav.awk < DIAG
# where DIAG is the diagnostic output from ChaNGa
#
BEGIN { rung = 0;
    maxrung = 30;
    for (tmp = 0; tmp < maxrung; tmp++) {
	counts[tmp] = 0;
	gravity[tmp] = 0.0;
	dd[tmp] = 0.0;
	loadb[tmp] = 0.0;
	build[tmp] = 0.0;
    }
	tgrav = tdd = tloadb = tbuild = 0.0;
	rung = 0;
	tsf = 0.0;
}
/Gravity Active/ {
    rung = $6
    }
/^total / { dd[rung] = dd[rung] + $2; }
/^Domain decomposition ... total/ { dd[rung] = dd[rung] + $5; }
/^took / { loadb[rung] = loadb[rung] + $2; }
/Building trees / { build[rung] = build[rung] + $5; }
/Calculating gravity / { counts[rung]++ ;
    if ($NF == 9) gravity[rung] = gravity[rung] + $8;
    else  gravity[rung] = gravity[rung] + $10; }
END {
    while(counts[maxrung-1] == 0) {
	maxrung--;	 
	}
    steps = counts[maxrung-1];
    
    print "Rung, counts, Gravity, DomainD, LoadB and TBuild";
    for (tmp = 0; tmp < maxrung; tmp++) {
	print tmp, counts[tmp], gravity[tmp], dd[tmp], loadb[tmp], build[tmp];
	tgrav = tgrav + gravity[tmp];
	tdd = tdd + dd[tmp];
	tloadb = tloadb + loadb[tmp];
	tbuild = tbuild + build[tmp];
    }
    print "Totals\nGrav, DomainD, LoadB and TBuild";
    print  tgrav, tdd, tloadb, tbuild;
    }
