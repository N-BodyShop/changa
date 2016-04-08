# Awk script to get time profile of multistep run
# This version is for "consph" runs
# Run this script with:
# awk -f timings_sph.awk < DIAG
# where DIAG is the diagnostic output from ChaNGa
#
BEGIN { rung = 0;
    maxrung = 30;
    first = 1;
    for (tmp = 0; tmp < maxrung; tmp++) {
	counts[tmp] = 0;
	gravity[tmp] = 0.0;
	dd[tmp] = 0.0;
	loadb[tmp] = 0.0;
	build[tmp] = 0.0;
	udot[tmp] = 0.0;
	adjust[tmp] = 0.0;
	drift[tmp] = 0.0;
	kick[tmp] = 0.0;
        cache[tmp] = 0.0;
    }
	tgrav = tudot = tdd = tloadb = tbuild = 0.0;
	tadjust = tdrift = tkick = tcache = 0.0;
	rung = 0;
	tsf = 0.0;
	tfb = 0.0;
	doSF = 0;
        tsfDD = 0.0;
        tsfLB = 0.0;
}
/Gravity Active/ {
    rung = $6;
    doSF = 0;
    }
/for star formation/ {
    doSF = 1;
    }
/^total / { if(!doSF) dd[rung] = dd[rung] + $2; else tsfDD = tsfDD + $2;}
/^Domain decomposition ... total/ { dd[rung] = dd[rung] + $5; }
/^Load balancer ... took / { loadb[rung] = loadb[rung] + $5; }
/^took / { if(!doSF) loadb[rung] = loadb[rung] + $2; else tsfLB = tsfLB + $2;}
/Building trees / { build[rung] = build[rung] + $5; }
/Calculating gravity and SPH / { if(first) { print "First Gravity ", $6; first = 0;}
                                 else { counts[rung]++ ; gravity[rung] = gravity[rung] + $6; }}
/^uDot/ { udot[rung] += $7 ; }
/^Star Formation/ { tsf += $5 ; }
/^Distribute/ { tfb += $9 ; }
/^Kick took/ { kick[rung]  += $3; }
/^Adjust took/ { adjust[rung]  += $3; }
/^Drift took/ { drift[rung]  += $3; }
/^Finish NodeCache took/ { cache[rung]  += $4; }
END {
    while(counts[maxrung-1] == 0) {
	maxrung--;	 
	}
    steps = counts[maxrung-1];
    
    print "Rung, counts, Gravity/SPH, uDot, DomainD, LoadB, TBuild, Adjust, Kick, Drift, Cache";
    for (tmp = 0; tmp < maxrung; tmp++) {
	print tmp, counts[tmp], gravity[tmp], udot[tmp], dd[tmp], loadb[tmp], build[tmp], adjust[tmp], kick[tmp], drift[tmp], cache[tmp];
	tgrav = tgrav + gravity[tmp];
 	tudot += udot[tmp];
	tdd = tdd + dd[tmp];
	tloadb = tloadb + loadb[tmp];
	tbuild = tbuild + build[tmp];
	tadjust += adjust[tmp];
        tkick += kick[tmp];
	tdrift += drift[tmp];
	tcache += cache[tmp];
    }
    print "Totals\nGrav, Star Form, FeedBack, udot, DomainD, LoadB, TBuild, Adjust, Kick, Drift, Cache";
    print  tgrav, tsf, tfb, tudot, tdd, tloadb, tbuild, tadjust, tkick, tdrift, tcache;
    print "StarFormDD ", tsfDD, " StarFormLB ", tsfLB;
    }
