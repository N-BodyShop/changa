/*
 * Module written for ChaNGa to calculate Lyman Werner feedback.
 * Added by Elaad Applebaum.
 */

#ifdef STOCH24
#define ARRLENGTH 24
#else
#define ARRLENGTH 12
#endif

#include <assert.h>
#include "lymanwerner.h"
#include <stdio.h>

/* Initialize data for LW table */
LWDATA *LymanWernerTableInit()
{
    LWDATA *lwd;
    lwd = (LWDATA *) malloc(sizeof(LWDATA));
    assert(lwd!=NULL);

    return lwd;
}

void lwInitData(LWDATA *lwd)
{
    int nrows;
    int ncols;
    FILE *fp;

    fp = fopen("lwtable.txt", "r");
    assert(fp!=NULL);
    fscanf(fp, "%d %d \n", &nrows, &ncols);
    lwd->lwLuminosity = (float **)malloc(nrows*sizeof(float**));
    for(int i=0; i<nrows; i++){
        lwd->lwLuminosity[i] = (float *)malloc(ncols*sizeof(float *));
    }
    

    for(int i=0; i<nrows; i++){
        for(int j=0; j<ncols; j++){
            fscanf(fp, "%f", &lwd->lwLuminosity[i][j]);
        }
    }

    fclose(fp);

}

void LymanWernerTableFinalize(LWDATA *lwd)
{
    if (lwd->lwLuminosity != NULL) free(lwd->lwLuminosity);
    free(lwd);
}

double calcLogStochLymanWerner(double dAgelog, SFEvent *sfEvent, LWDATA *lwd)
{
/* Calculating the LW flux of individual high mass stars, by using the closest values
 * (interpolation would largely be pointless, since the values are either flat or uncertain anyway)
 *
 * The masses, times, and corresponding intervals between them are based on lwtable.txt,
 * which has nonuniform masses. See the README in the data directory.*/
    double dAge = pow(10, dAgelog);

    double fluxtot = 0.0;

    for (int i=0; i<ARRLENGTH; i++){
        if (sfEvent->rgdHMStars[i]>0){
            int tindex;
            int mindex;
            if (dAge < 1.0E6) tindex = 0;
            else if (dAge >= 9.701E7) index = 49; // beyond ~10^7.5, all hmstars have constant LW flux
            else tindex = (int)round((dAge-1.0E4)/2.0E6);
            CkPrintf("%f, %d\n", dAge, tindex);
            
            double dMass = sfEvent->rgdHMStars[i];
            if (dMass >= 98.5) mindex = 61;
            else if (dMass < 39.5) mindex = (int)round(dMass-8);
            else mindex = (int)round((dMass+23.0)/2.0);
            CkPrintf("%f, %d\n", dMass, mindex);
            
            fluxtot += pow(10, lwd->lwLuminosity[mindex][tindex]);
        }
    }
    return log10(fluxtot);
}
