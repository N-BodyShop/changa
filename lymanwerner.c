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
#include <charm.h>

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
    int i, j;
    FILE *fp;

    fp = fopen("lwtable.txt", "r");
    assert(fp!=NULL);
    fscanf(fp, "%d %d \n", &(lwd->nrows), &(lwd->ncols));
    lwd->lwLuminosity = (float **)malloc(lwd->nrows*sizeof(float**));
    for(i=0; i<lwd->nrows; i++){
        lwd->lwLuminosity[i] = (float *)malloc(lwd->ncols*sizeof(float *));
    }
    

    for(i=0; i<lwd->nrows; i++){
        for(j=0; j<lwd->ncols; j++){
            fscanf(fp, "%f", &lwd->lwLuminosity[i][j]);
        }
    }

    fclose(fp);

}

void LymanWernerTableFinalize(LWDATA *lwd)
{
    int i;
    if (lwd->lwLuminosity != NULL) free(lwd->lwLuminosity);
    for(i=0; i<lwd->nrows; i++){
        if (lwd->lwLuminosity[i] != NULL) free(lwd->lwLuminosity[i]);
    }
    free(lwd);
}

/// @brief Calculate LW flux from a set of HM stars
/// @param dAgelog log10 age of stars in years
/// @param rgdHMStars array of HM star masses
/// @param lwd LW flux table
double calcLogStochLymanWerner(double dAgelog, double *rgdHMStars, LWDATA *lwd)
{
/* Calculating the LW flux of individual high mass stars, by using the closest values
 * (interpolation would largely be pointless, since the values are either flat or uncertain anyway)
 *
 * The masses, times, and corresponding intervals between them are based on lwtable.txt,
 * which has nonuniform masses. See the README in the data directory.*/
    double dAge = pow(10, dAgelog);
    int i;
    double fluxtot = 0.0;
    double minAge = 1.e6; // Minimum age (gives tindex = 0)
    double maxAge = 9.701e7; // Maximum age (gives tindex = lwd->ncols-1)
    double tOff = 1.0e4; // Offset to align index in LWDATA table
    double tStep = 2.0e6; // Step size in LWDATA table
    
    int tindex;
    if (dAge < minAge) tindex = 0;
    else if (dAge >= maxAge) tindex = lwd->ncols-1; // beyond ~10^7.5, all hmstars have constant LW flux
    else tindex = (int)round((dAge-tOff)/tStep);

    for (i=0; i<ARRLENGTH; i++){
        if (rgdHMStars[i]>=8){
            int mindex;
            
            double dMass = rgdHMStars[i];
            if (dMass >= 98.5) mindex = 61;
            else if (dMass < 39.5) mindex = (int)round(dMass-8);
            else mindex = (int)round((dMass+23.0)/2.0);
            CkAssert(mindex < lwd->nrows);
            CkAssert(tindex < lwd->ncols);
            
            fluxtot += pow(10, lwd->lwLuminosity[mindex][tindex]);
        }
    }
    return log10(fluxtot);
}
