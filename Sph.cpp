/*
 * Routines to implement SPH.
 * Main author: James Wadsley, as first implemented in GASOLINE.
 * See Wadsley, J.~W., Stadel, J., Quinn, T.\ 2004.\ Gasoline: a flexible,
 * parallel implementation of TreeSPH.\ New Astronomy 9, 137-158.
 */

#include "ParallelGravity.h"
#include "DataManager.h"
#include "smooth.h"
#include "Sph.h"
#include "SphUtils.h"
#include "physconst.h"

#ifndef MAXPATHLEN
#define MAXPATHLEN PATH_MAX
#endif

#include <float.h>

///
/// @brief initialize SPH quantities
///
/// Initial calculation of densities and internal energies, and cooling rates.
///
void
Main::initSph() 
{
    if(param.bDoGas) {
	ckout << "Calculating densities/divv ...";
	// The following smooths all GAS, and also marks neighbors of
	// actives, and those who have actives as neighbors
        // Starting is true
	DenDvDxSmoothParams pDen(TYPE_GAS, 0, param.csm, dTime, 0,
				 param.bConstantDiffusion, 1, bHaveAlpha,
                                 param.dConstAlphaMax);
	double startTime = CkWallTimer();
	double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
	treeProxy.startSmooth(&pDen, 1, param.nSmooth, dfBall2OverSoft2,
			      CkCallbackResumeThread());
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
        if(verbosity > 1 && !param.bConcurrentSph)
            memoryStatsCache();
	double dTuFac = param.dGasConst/(param.dConstGamma-1)
	    /param.dMeanMolWeight;
	double z = 1.0/csmTime2Exp(param.csm, dTime) - 1.0;
	if(param.bGasCooling) {
	    // Update cooling on the datamanager
	    dMProxy.CoolingSetTime(z, dTime, CkCallbackResumeThread());
	    if(!bIsRestarting)  // Energy is already OK from checkpoint.
		treeProxy.InitEnergy(dTuFac, z, dTime, CkCallbackResumeThread());
	    }
	if(verbosity) CkPrintf("Initializing SPH forces\n");
	nActiveSPH = nTotalSPH;
	doSph(0, 0);
	double duDelta[MAXRUNG+1];
	double dStartTime[MAXRUNG+1];
	for(int iRung = 0; iRung <= MAXRUNG; iRung++) {
	    duDelta[iRung] = 0.5e-7*param.dDelta;
	    dStartTime[iRung] = dTime;
	    }
	treeProxy.updateuDot(0, duDelta, dStartTime, param.bGasCooling, 0, 1,
            (param.dConstGamma-1), CkCallbackResumeThread());
	}
    }

// see below for definition.
bool arrayFileExists(const std::string filename, const int64_t count) ;

#include <sys/stat.h>

///
/// @brief Initialize cooling constants and integration data structures.
///
void Main::initCooling()
{
#ifndef COOLING_NONE
    dMProxy.initCooling(param.dGmPerCcUnit, param.dComovingGmPerCcUnit,
		    param.dErgPerGmUnit, param.dSecUnit, param.dKpcUnit,
		    param.CoolParam, CkCallbackResumeThread());
    
    /* Read in tables from files as necessary */
    int cntTable = 0;
    int nTableRows;
    int nTableColumns;
    char TableFileSuffix[20];
    
    for (;;) {
	CoolTableReadInfo(&param.CoolParam, cntTable, &nTableColumns,
			  TableFileSuffix);
	if (!nTableColumns) break;
    
	cntTable++;
	nTableRows = ReadASCII(TableFileSuffix, nTableColumns, NULL);
	if (nTableRows) {
	    CkAssert(sizeof(double)*nTableRows*nTableColumns <= CL_NMAXBYTETABLE );
	    double *dTableData = (double *)malloc(sizeof(double)*nTableRows*nTableColumns);
	    CkAssert( dTableData != NULL );
	    nTableRows = ReadASCII(TableFileSuffix, nTableColumns, dTableData);
      
	    dMProxy.dmCoolTableRead(dTableData,nTableRows*nTableColumns,
				  CkCallbackResumeThread());
	    free(dTableData);
	    }
	}
    treeProxy.initCoolingData(CkCallbackResumeThread());
    if(!bIsRestarting) {  // meaning not restarting from a checkpoint.
        struct stat s;
        int err = stat(basefilename.c_str(), &s);
        if(err != -1 && S_ISDIR(s.st_mode)) {
            // The file is a directory; assume NChilada
            int64_t nGas = 0;
            nGas = ncGetCount(basefilename + "/gas/coolontime");
            if(nGas == nTotalSPH) {
                CkPrintf("Reading coolontime\n");
                coolontimeOutputParams pCoolOnOut(basefilename, 6, 0.0);
                treeProxy.readFloatBinary(pCoolOnOut, param.bParaRead,
                                          CkCallbackResumeThread());
                }
            }
        else {
            if(arrayFileExists(basefilename + ".coolontime", nTotalParticles)) {
                CkPrintf("Reading coolontime\n");
                coolontimeOutputParams pCoolOnOut(basefilename, 0, 0.0);
                treeProxy.readTipsyArray(pCoolOnOut, CkCallbackResumeThread());
                }
            }
        }
#endif
    }

/**
 * Initialized Cooling Read-only data on the DataManager, which
 * doesn't migrate.
 */
void
DataManager::initCooling(double dGmPerCcUnit, double dComovingGmPerCcUnit,
		       double dErgPerGmUnit, double dSecUnit, double dKpcUnit,
		       COOLPARAM inParam, const CkCallback& cb)
{
#ifndef COOLING_NONE
    clInitConstants(Cool, dGmPerCcUnit, dComovingGmPerCcUnit, dErgPerGmUnit,
		    dSecUnit, dKpcUnit, inParam);
    
    CoolInitRatesTable(Cool,inParam);
#endif
    contribute(cb);
    }

/**
 * Per thread initialization
 */
void
TreePiece::initCoolingData(const CkCallback& cb)
{
#ifndef COOLING_NONE
    bGasCooling = 1;
    dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    CoolData = CoolDerivsInit(dm->Cool);
#endif
    contribute(cb);
    }

void
DataManager::dmCoolTableRead(double *dTableData, int nData, const CkCallback& cb)
{
#ifndef COOLING_NONE
    CoolTableRead(Cool, nData*sizeof(double), (void *) dTableData);
#endif
    contribute(cb);
    }

///
/// @brief function from PKDGRAV to read an ASCII table
///
/// @param extension Appended to outName to determine file name to
/// read.
/// @param nDataPerLine Number of columns in the table.
/// @param dDataOut pointer to array in which to store the table.
/// Note if dDataOut is NULL it just counts the number of valid input
/// lines.
///
int Main::ReadASCII(char *extension, int nDataPerLine, double *dDataOut)
{
	FILE *fp;
	int i,ret;
	char achIn[160];
	double *dData;

	if (dDataOut == NULL) 
	    dData = (double *)malloc(sizeof(double)*nDataPerLine);
	else
	    dData = dDataOut;
	
	CkAssert(nDataPerLine > 0 && nDataPerLine <= 10);
	char achFile[MAXPATHLEN];
	sprintf(achFile, "%s.%s", param.achOutName, extension);
	fp = fopen(achFile,"r");
	if (!fp) {
            CkPrintf("WARNING: Could not open .%s input file:%s\n",
                      extension,achFile);
	    return 0;
	    }

	i = 0;
	while (1) {
	    if (!fgets(achIn,160,fp)) goto Done;
	    switch (nDataPerLine) {
	    case 1:
		ret = sscanf(achIn,"%lf",dData); 
		break;
	    case 2:
		ret = sscanf(achIn,"%lf %lf",dData,dData+1); 
		break;
	    case 3:
		ret = sscanf(achIn,"%lf %lf %lf",dData,dData+1,dData+2); 
		break;
	    case 4:
		ret = sscanf(achIn,"%lf %lf %lf %lf",dData,dData+1,dData+2,dData+3); 
		break;
	    case 5:
		ret = sscanf(achIn,"%lf %lf %lf %lf %lf",dData,dData+1,dData+2,dData+3,dData+4); 
		break;
	    case 6:
		ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf",dData,dData+1,dData+2,dData+3,dData+4,dData+5); 
		break;
	    case 7:
		ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf",
			     dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6); 
		break;
	    case 8:
		ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf %lf",
			     dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6,dData+7); 
		break;
	    case 9:
		ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
			     dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6,dData+7,dData+8); 
		break;
	    case 10:
		ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			     dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6,dData+7,dData+8,dData+9); 
		break;
	    default:
		ret = EOF;
		CkAssert(0);
		}
	    if (ret != nDataPerLine) goto Done;
	    ++i;
	    if (dDataOut != NULL) dData += nDataPerLine;
	    }
 Done:
	fclose(fp);
	if (dDataOut != NULL && verbosity)
	    printf("Read %i lines from %s\n",i,achFile);
	if (dDataOut == NULL) free(dData);
	return i;
	}

/*
 * Update the cooling functions to the current time.
 * This is on the DataManager to avoid duplication of effort.
 */
void
DataManager::CoolingSetTime(double z, // redshift
			    double dTime, // Time
			    const CkCallback& cb)
{
#ifndef COOLING_NONE
    CoolSetTime( Cool, dTime, z  );
#endif

    contribute(cb);
    }

/**
 * @brief DataManager::SetStarCM saves the total mass and center of mass of the
 * star(s) to the COOL struct Cool, making them available to the cool particles
 * @param dCenterOfMass Array(length 4) which contains the star(s) center of
 * mass as the first 3 entries and the total star mass as the final entry
 * @param cb    Callback
 */
void DataManager::SetStarCM(double dCenterOfMass[4], const CkCallback& cb) {
#ifndef COOLING_NONE
#ifdef COOLING_PLANET
    CoolSetStarCM(Cool, dCenterOfMass);
#endif
#endif
    contribute(cb);
}

/**
 *  @brief utility for checking array files
 */
bool
arrayFileExists(const std::string filename, const int64_t count) 
{
    FILE *fp = CmiFopen(filename.c_str(), "r");
    if(fp != NULL) {
        // Check if its a binary file
        unsigned int iDum;
        XDR xdrs;
        xdrstdio_create(&xdrs, fp, XDR_DECODE);
        xdr_u_int(&xdrs,&iDum);
        xdr_destroy(&xdrs);
        if(iDum == count) { // Assume a valid binary array file
            fclose(fp);
            return true;
            }
        fseek(fp, 0, SEEK_SET);
        int64_t nIOrd;
        fscanf(fp, "%ld", &nIOrd);
        CkAssert(nIOrd == count); // Valid ASCII file.
        fclose(fp);
        return true;
    }
    return false;
}

/// @brief Set total metals based on Ox and Fe mass fractions
void
TreePiece::resetMetals(const CkCallback& cb)
{
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
        // Use total metals to Fe and O based on Asplund et al 2009
	if (p->isGas())
            p->fMetals() = 1.06*p->fMFracIron() + 2.09*p->fMFracOxygen();
	if (p->isStar())
            p->fStarMetals() = 1.06*p->fStarMFracIron()
                + 2.09*p->fStarMFracOxygen();
        }

    contribute(cb);
}

#include <sys/stat.h>
/**
 *  @brief Read in array files for complete gas information.
 */
void
Main::restartGas() 
{
    if(verbosity)
        CkPrintf("Restarting Gas Simulation with array files.\n");
    
  struct stat s;
  int err = stat(basefilename.c_str(), &s);
  if(err != -1 && S_ISDIR(s.st_mode)) {
      // The file is a directory; assume NChilada
      int64_t nGas = 0;
      int64_t nDark = 0;
      int64_t nStar = 0;
      if(nTotalSPH > 0)
          nGas = ncGetCount(basefilename + "/gas/iord");
      if(nTotalDark > 0)
          nDark = ncGetCount(basefilename + "/dark/iord");
      if(nTotalStar > 0)
          nStar = ncGetCount(basefilename + "/star/iord");
      if(nGas + nDark + nStar == nTotalParticles) {
          IOrderOutputParams pIOrdOut(basefilename, 6, 0.0);
          treeProxy.readFloatBinary(pIOrdOut, param.bParaRead,
                                      CkCallbackResumeThread());
          CkReductionMsg *msg;
          treeProxy.getMaxIOrds(CkCallbackResumeThread((void*&)msg));
          CmiInt8 *maxIOrds = (CmiInt8 *)msg->getData();
          nMaxOrderGas = maxIOrds[0];
          nMaxOrderDark = maxIOrds[1];
          nMaxOrder = maxIOrds[2];
          delete msg;
          }
      else
          CkError("WARNING: no iorder file, or wrong format for restart\n");
      if(nTotalStar > 0)
          nStar = ncGetCount(basefilename + "/star/igasorder");
      if(nStar == nTotalStar) {
          IGasOrderOutputParams pIOrdOut(basefilename, 6, 0.0);
          treeProxy.readFloatBinary(pIOrdOut, param.bParaRead,
                                      CkCallbackResumeThread());
          }
      else
          CkError("WARNING: no igasorder file, or wrong format for restart\n");
      if(param.bFeedback) {
          if(nTotalSPH > 0)
              nGas = ncGetCount(basefilename + "/gas/ESNRate");
          if(nTotalStar > 0)
              nStar = ncGetCount(basefilename + "/star/ESNRate");
          if(nGas + nStar == nTotalSPH + nTotalStar) {
              ESNRateOutputParams pESNROut(basefilename, 6, 0.0);
              treeProxy.readFloatBinary(pESNROut, param.bParaRead,
                                   CkCallbackResumeThread());
              }
          else
              CkError("WARNING: no ESNRate file, or wrong format for restart\n");
          if(nTotalSPH > 0)
              nGas = ncGetCount(basefilename + "/gas/OxMassFrac");
          if(nTotalStar > 0)
              nStar = ncGetCount(basefilename + "/star/OxMassFrac");
          if(nGas + nStar == nTotalSPH + nTotalStar) {
              OxOutputParams pOxOut(basefilename, 6, 0.0);
              treeProxy.readFloatBinary(pOxOut, param.bParaRead,
                                   CkCallbackResumeThread());
              }
          else
              CkError("WARNING: no OxMassFrac file, or wrong format for restart\n");
          if(nTotalSPH > 0)
              nGas = ncGetCount(basefilename + "/gas/FeMassFrac");
          if(nTotalStar > 0)
              nStar = ncGetCount(basefilename + "/star/FeMassFrac");
          if(nGas + nStar == nTotalSPH + nTotalStar) {
              FeOutputParams pFeOut(basefilename, 6, 0.0);
              treeProxy.readFloatBinary(pFeOut, param.bParaRead,
                                   CkCallbackResumeThread());
              }
          else
              CkError("WARNING: no FeMassFrac file, or wrong format for restart\n");
          treeProxy.resetMetals(CkCallbackResumeThread());

          if(nTotalStar > 0)
              nStar = ncGetCount(basefilename + "/star/massform");
          if(nStar == nTotalStar) {
              MFormOutputParams pMFOut(basefilename, 6, 0.0);
              treeProxy.readFloatBinary(pMFOut, param.bParaRead,
                                   CkCallbackResumeThread());
              }
          else
              CkError("WARNING: no massform file, or wrong format for restart\n");
          }
#ifdef CULLENALPHA
      if(nTotalSPH > 0) {
          nGas = ncGetCount(basefilename + "/gas/alpha");
          if(nGas == nTotalSPH) {
              AlphaOutputParams pAlphaOut(basefilename, 6, 0.0);
              treeProxy.readFloatBinary(pAlphaOut, param.bParaRead,
                                        CkCallbackResumeThread());
              bHaveAlpha = 1;
              }
          else
              CkError("WARNING: no alpha file, or wrong format for restart\n");
          }
#endif
#ifndef COOLING_NONE
      if(param.bGasCooling && nTotalSPH > 0) {
          bool bFoundCoolArray = false;
          // read ionization fractions
          nGas = ncGetCount(basefilename + "/gas/" + COOL_ARRAY0_EXT);
          if(nGas == nTotalSPH) {
              Cool0OutputParams pCool0Out(basefilename, 6, 0.0);
              treeProxy.readFloatBinary(pCool0Out, param.bParaRead,
                                   CkCallbackResumeThread());
              bFoundCoolArray = true;
              }
          else
              CkError("WARNING: no CoolArray0 file, or wrong format for restart\n");
          nGas = ncGetCount(basefilename + "/gas/" + COOL_ARRAY1_EXT);
          if(nGas == nTotalSPH) {
              Cool1OutputParams pCool1Out(basefilename, 6, 0.0);
              treeProxy.readFloatBinary(pCool1Out, param.bParaRead,
                                   CkCallbackResumeThread());
              bFoundCoolArray = true;
              }
          else
              CkError("WARNING: no CoolArray1 file, or wrong format for restart\n");
          nGas = ncGetCount(basefilename + "/gas/" + COOL_ARRAY2_EXT);
          if(nGas == nTotalSPH) {
              Cool2OutputParams pCool2Out(basefilename, 6, 0.0);
              treeProxy.readFloatBinary(pCool2Out, param.bParaRead,
                                   CkCallbackResumeThread());
              bFoundCoolArray = true;
              }
          else
              CkError("WARNING: no CoolArray2 file, or wrong format for restart\n");
          nGas = ncGetCount(basefilename + "/gas/" + COOL_ARRAY3_EXT);
          if(nGas == nTotalSPH) {
              Cool3OutputParams pCool3Out(basefilename, 6, 0.0);
              treeProxy.readFloatBinary(pCool3Out, param.bParaRead,
                                   CkCallbackResumeThread());
              bFoundCoolArray = true;
              }
          else
              CkError("WARNING: no CoolArray3 file, or wrong format for restart\n");
        double dTuFac = param.dGasConst/(param.dConstGamma-1)
                /param.dMeanMolWeight;
        if(bFoundCoolArray) {
            // reset thermal energy with ionization fractions
            treeProxy.RestartEnergy(dTuFac, CkCallbackResumeThread());
        }
        else {
            double z = 1.0/csmTime2Exp(param.csm, dTime) - 1.0;
            dMProxy.CoolingSetTime(z, dTime, CkCallbackResumeThread());
            treeProxy.InitEnergy(dTuFac, z, dTime, CkCallbackResumeThread());
            }
      }
#endif
  } else {                      
    // Assume TIPSY arrays
    // read iOrder
    if(arrayFileExists(basefilename + ".iord", nTotalParticles)) {
        CkReductionMsg *msg;
        IOrderOutputParams pIOrdOut(basefilename, 0, 0.0);
        treeProxy.readTipsyArray(pIOrdOut, CkCallbackResumeThread());

        treeProxy.getMaxIOrds(CkCallbackResumeThread((void*&)msg));
        CmiInt8 *maxIOrds = (CmiInt8 *)msg->getData();
        nMaxOrderGas = maxIOrds[0];
        nMaxOrderDark = maxIOrds[1];
        nMaxOrder = maxIOrds[2];
        delete msg;
        }
    else
        CkError("WARNING: no iOrder file for restart\n");
    // read iGasOrder
    if(arrayFileExists(basefilename + ".igasorder", nTotalParticles)) {
        IGasOrderOutputParams pIOrdOut(basefilename, 0, 0.0);
        treeProxy.readTipsyArray(pIOrdOut, CkCallbackResumeThread());
        }
    else {
        CkError("WARNING: no igasorder file for restart\n");
        }
    if(param.bFeedback) {
        if(arrayFileExists(basefilename + ".ESNRate", nTotalParticles)) {
            ESNRateOutputParams pESNROut(basefilename, 0, 0.0);
            treeProxy.readTipsyArray(pESNROut, CkCallbackResumeThread());
            }
        if(arrayFileExists(basefilename + ".OxMassFrac", nTotalParticles)) {
            OxOutputParams pOxOut(basefilename, 0, 0.0);
            treeProxy.readTipsyArray(pOxOut, CkCallbackResumeThread());
            }
        if(arrayFileExists(basefilename + ".FeMassFrac", nTotalParticles)) {
            FeOutputParams pFeOut(basefilename, 0, 0.0);
            treeProxy.readTipsyArray(pFeOut, CkCallbackResumeThread());
            }
        treeProxy.resetMetals(CkCallbackResumeThread());
        if(arrayFileExists(basefilename + ".massform", nTotalParticles)) {
            MFormOutputParams pMFOut(basefilename, 0, 0.0);
            treeProxy.readTipsyArray(pMFOut, CkCallbackResumeThread());
            }
        }
#ifdef CULLENALPHA
        if(arrayFileExists(basefilename + ".alpha", nTotalParticles)) {
            AlphaOutputParams pAlphaOut(basefilename, 0, 0.0);
            treeProxy.readTipsyArray(pAlphaOut, CkCallbackResumeThread());
            bHaveAlpha = 1;
            }
        else
            CkError("WARNING: no alpha file, or wrong format for restart\n");
#endif
#ifndef COOLING_NONE
    if(param.bGasCooling) {
        bool bFoundCoolArray = false;
        // read ionization fractions
        if(arrayFileExists(basefilename + "." + COOL_ARRAY0_EXT, nTotalParticles)) {
            Cool0OutputParams pCool0Out(basefilename, 0, 0.0);
            treeProxy.readTipsyArray(pCool0Out, CkCallbackResumeThread());
            bFoundCoolArray = true;
            }
        else {
            CkError("WARNING: no CoolArray0 file for restart\n");
            }
        if(arrayFileExists(basefilename + "." + COOL_ARRAY1_EXT, nTotalParticles)) {
                
            Cool1OutputParams pCool1Out(basefilename, 0, 0.0);
            treeProxy.readTipsyArray(pCool1Out, CkCallbackResumeThread());
            bFoundCoolArray = true;
            }
        else {
            CkError("WARNING: no CoolArray1 file for restart\n");
            }
        if(arrayFileExists(basefilename + "." + COOL_ARRAY2_EXT, nTotalParticles)) {
            Cool2OutputParams pCool2Out(basefilename, 0, 0.0);
            treeProxy.readTipsyArray(pCool2Out, CkCallbackResumeThread());
            bFoundCoolArray = true;
            }
        else {
            CkError("WARNING: no CoolArray2 file for restart\n");
            }
        if(arrayFileExists(basefilename + "." + COOL_ARRAY3_EXT, nTotalParticles)) {
            Cool3OutputParams pCool3Out(basefilename, 0, 0.0);
            treeProxy.readTipsyArray(pCool3Out, CkCallbackResumeThread());
            bFoundCoolArray = true;
            }
        else {
            CkError("WARNING: no CoolArray3 file for restart\n");
            }
        double dTuFac = param.dGasConst/(param.dConstGamma-1)
                /param.dMeanMolWeight;
        if(bFoundCoolArray) {
            // reset thermal energy with ionization fractions
            treeProxy.RestartEnergy(dTuFac, CkCallbackResumeThread());
        }
        else {
            double z = 1.0/csmTime2Exp(param.csm, dTime) - 1.0;
            dMProxy.CoolingSetTime(z, dTime, CkCallbackResumeThread());
            treeProxy.InitEnergy(dTuFac, z, dTime, CkCallbackResumeThread());
            }
        }
#endif
  }
}

/*
 * Initialize energy on restart
 */
void TreePiece::RestartEnergy(double dTuFac, // T to internal energy
                              const CkCallback& cb)
{
#ifndef COOLING_NONE
    COOL *cl;

    dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    cl = dm->Cool;
#endif

    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
	if (p->isGas()) {
#ifndef COOLING_NONE
#ifndef COOLING_GRACKLE
	    double T;
	    T = p->u() / dTuFac;
            PERBARYON Y;
#ifdef COOLING_METAL
            CoolPARTICLEtoPERBARYON(cl, &Y, &p->CoolParticle(), p->fMetals());
#elif COOLING_MOLECULARH
            CoolPARTICLEtoPERBARYON(cl, &Y, &p->CoolParticle(), p->fMetals());
#else
            CoolPARTICLEtoPERBARYON(cl, &Y, &p->CoolParticle());
#endif
            
	    p->u() = clThermalEnergy(Y.Total,T)*cl->diErgPerGmUnit;
#endif
#endif
	    p->uPred() = p->u();
	    }
	}
    contribute(cb);
    }

/**
 *  @brief Perform the SPH force calculation.
 *  @param activeRung Timestep rung (and above) on which to perform
 *  SPH
 *  @param bNeedDensity Does the density calculation need to be done?
 *  Defaults to 1
 */
void
Main::doSph(int activeRung, int bNeedDensity) 
{
  if(bNeedDensity) {
    double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
    if (param.bFastGas && nActiveSPH < nTotalSPH*param.dFracFastGas) {
	ckout << "Calculating densities/divv on Actives ...";
	// This also marks neighbors of actives
	DenDvDxSmoothParams pDen(TYPE_GAS, activeRung, param.csm, dTime, 1,
				 param.bConstantDiffusion, 0, 0,
                                 param.dConstAlphaMax);
	double startTime = CkWallTimer();
	treeProxy.startSmooth(&pDen, 1, param.nSmooth, dfBall2OverSoft2,
			      CkCallbackResumeThread());
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;

	ckout << "Marking Neighbors ...";
	// This marks particles with actives as neighbors
	MarkSmoothParams pMark(TYPE_GAS, activeRung);
	startTime = CkWallTimer();
	treeProxy.startMarkSmooth(&pMark, CkCallbackResumeThread());
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
	
	ckout << "Density of Neighbors ...";
	// This does neighbors (but not actives),  It also does no
	// additional marking
	DenDvDxNeighborSmParams pDenN(TYPE_GAS, activeRung, param.csm, dTime,
				      param.bConstantDiffusion,
                                      param.dConstAlphaMax);
	startTime = CkWallTimer();
	treeProxy.startSmooth(&pDenN, 1, param.nSmooth, dfBall2OverSoft2,
			      CkCallbackResumeThread());
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
	}
    else {
	ckout << "Calculating densities/divv ...";
	// The following smooths all GAS, and also marks neighbors of
	// actives, and those who have actives as neighbors.
	DenDvDxSmoothParams pDen(TYPE_GAS, activeRung, param.csm, dTime, 0,
				 param.bConstantDiffusion, 0, 0,
                                 param.dConstAlphaMax);
	double startTime = CkWallTimer();
	treeProxy.startSmooth(&pDen, 1, param.nSmooth, dfBall2OverSoft2,
			      CkCallbackResumeThread());
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;

	if(verbosity > 1 && !param.bConcurrentSph)
	    memoryStatsCache();
	}
      }
    treeProxy.sphViscosityLimiter(param.iViscosityLimiter, activeRung,
			CkCallbackResumeThread());
    double a = csmTime2Exp(param.csm,dTime);
    double dDtCourantFac = param.dEtaCourant*a*2.0/1.6;
    if(param.bGasCooling)
	treeProxy.getCoolingGasPressure(param.dConstGamma,
					param.dConstGamma-1, param.dThermalCondCoeffCode*a, param.dThermalCond2CoeffCode*a,
                    param.dThermalCondSatCoeff/a, param.dThermalCond2SatCoeff/a, 
            param.dEvapMinTemp,	dDtCourantFac,
            param.dResolveJeans/a,
            CkCallbackResumeThread());
    else
	treeProxy.getAdiabaticGasPressure(param.dConstGamma,
					param.dConstGamma-1, param.dThermalCondCoeffCode*a, param.dThermalCond2CoeffCode*a,
                    param.dThermalCondSatCoeff/a, param.dThermalCond2SatCoeff/a, 
                    param.dEvapMinTemp,	dDtCourantFac, CkCallbackResumeThread());

    ckout << "Calculating pressure gradients ...";
    PressureSmoothParams pPressure(TYPE_GAS, activeRung, param.csm, dTime,
                                   param.dConstAlpha, param.dConstBeta,
                                   param.dThermalDiffusionCoeff, param.dMetalDiffusionCoeff,
                                   param.dEtaCourant, param.dEtaDiffusion);
    double startTime = CkWallTimer();
    treeProxy.startReSmooth(&pPressure, CkCallbackResumeThread());
    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	  << endl;
    
    treeProxy.ballMax(activeRung, 1.0+param.ddHonHLimit,
		      CkCallbackResumeThread());
    }

/*
 * Initialize energy and ionization state for cooling particles
 */
void TreePiece::InitEnergy(double dTuFac, // T to internal energy
			   double z,	  // redshift
			   double dTime,
			   const CkCallback& cb)
{
#ifndef COOLING_NONE
    COOL *cl;

    dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    cl = dm->Cool;
#endif

    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
	if (TYPETest(p, TYPE_GAS) && p->rung >= activeRung) {
#ifndef COOLING_NONE
	    double T,E;
	    T = p->u() / dTuFac;
	    CoolInitEnergyAndParticleData(cl, &p->CoolParticle(), &E,
					  p->fDensity, T, p->fMetals() );
	    p->u() = E;
#endif
	    p->uPred() = p->u();
	    }
	}
    // Use shadow array to avoid reduction conflict
    smoothProxy[thisIndex].ckLocal()->contribute(cb);
    }

/**
 * @brief Update the cooling rate (uDot)
 *
 * @param activeRung (minimum) rung being updated
 * @param duDelta    array of timesteps of length MAXRUNG+1
 * @param dStartTime array of start times of length MAXRUNG+1
 * @param bCool      Whether cooling is on
 * @param bUpdateState Whether the ionization factions need updating
 * @param bAll	     Do all rungs below activeRung
 * @param cb	     Callback.
 */
void TreePiece::updateuDot(int activeRung,
			   double duDelta[MAXRUNG+1], // timesteps
			   double dStartTime[MAXRUNG+1],
			   int bCool, // select equation of state
			   int bUpdateState, // update ionization fractions
			   int bAll, // update all rungs below activeRung
               double gammam1,
			   const CkCallback& cb)
{
#ifndef COOLING_NONE
    double dt; // time in seconds
    double fDensity;
    double E;
    double ExternalHeating;
    
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
	if (TYPETest(p, TYPE_GAS)
	    && (p->rung == activeRung || (bAll && p->rung >= activeRung))) {
	    dt = CoolCodeTimeToSeconds(dm->Cool, duDelta[p->rung] );
        fDensity = p->fDensity;
        ExternalHeating = p->PdV() + p->fESNrate();
	    if ( bCool ) {
		COOLPARTICLE cp = p->CoolParticle();
		double r[3];  // For conversion to C
		p->position.array_form(r);
        CkAssert(p->u() < LIGHTSPEED*LIGHTSPEED/dm->Cool->dErgPerGmUnit);
        CkAssert(p->uPred() < LIGHTSPEED*LIGHTSPEED/dm->Cool->dErgPerGmUnit);
#ifdef SUPERBUBBLE
        double frac = p->massHot()/p->mass;
        double PoverRho = gammam1*(p->uHotPred()*frac+p->uPred()*(1-frac));
        double uMean = frac*p->uHotPred()+(1-frac)*p->uPred();
        CkAssert(uMean > 0.0);
        CkAssert(p->uHotPred() < LIGHTSPEED*LIGHTSPEED/dm->Cool->dErgPerGmUnit);
        CkAssert(p->uHot() < LIGHTSPEED*LIGHTSPEED/dm->Cool->dErgPerGmUnit);
        /*
         * If we have mass in the hot phase, we need to cool it appropriately.
         */
        if (p->massHot() > 0) { 
            ExternalHeating = p->PdV()*p->uHotPred()/uMean + p->fESNrate();
            if (p->uHot() > 0) {
                E = p->uHot();
                fDensity = p->fDensity*PoverRho/(gammam1*p->uHot());
                cp = p->CoolParticleHot();
#ifdef COOLING_MOLECULARH
                // Assume the cold phase is a shell surrounding the hot phase,
                // which is a sphere
                double columnL = pow(0.75*p->massHot()/(M_PI*fDensity), 0.333);
#ifdef COOLDEBUG
                dm->Cool->iOrder = p->iOrder; /*For debugging purposes */
#endif
                CoolIntegrateEnergyCode(dm->Cool, CoolData, &cp, &E,
                            ExternalHeating, fDensity,
                            p->fMetals(), r, dt, columnL);
#else /*COOLING_MOLECULARH*/
                CoolIntegrateEnergyCode(dm->Cool, CoolData, &cp, &E, ExternalHeating, fDensity,
                        p->fMetals(), r, dt);
#endif
                p->uHotDot() = (E- p->uHot())/duDelta[p->rung];
                if(bUpdateState) p->CoolParticleHot() = cp;
            }
            else /* If we just got feedback, only set up the uDot */
            {
                p->uHotDot() = ExternalHeating;
                p->cpHotInit() = 1;
            }
            ExternalHeating = p->PdV()*p->uPred()/uMean;
        }
        else { /* We have a single phase particle, treat it normally*/
            p->uHotDot() = 0;
            ExternalHeating = p->PdV() + p->fESNrate();
        }
        fDensity = p->fDensity*PoverRho/(gammam1*p->uPred());
        if (p->fDensityU() < p->fDensity) fDensity = p->fDensityU()*PoverRho/(gammam1*p->uPred());
        CkAssert(fDensity > 0);
        cp = p->CoolParticle();
#endif
		E = p->u();
#ifdef COOLING_BOLEY
		cp.mrho = pow(p->mass/p->fDensity, 1./3.);
#endif
		double dtUse = dt;
		
		if(dStartTime[p->rung] + 0.5*duDelta[p->rung]
		   < p->fTimeCoolIsOffUntil()) {
		    /* This flags cooling shutoff (e.g., from SNe) to
		       the cooling functions. */
		    dtUse = -dt;
		    p->uDot() = ExternalHeating;
		    }

#ifdef COOLING_MOLECULARH
		/*		cp.dLymanWerner = 52.0; for testing CC */
		double columnL = sqrt(0.25)*p->fBall;
#ifdef SUPERBUBBLE
        // Assume the cold phase is a shell surrounding the hot phase,
        // which is a sphere
        columnL = p->massHot()/(4*M_PI*columnL*columnL*fDensity);
#endif
#ifdef COOLDEBUG
		dm->Cool->iOrder = p->iOrder; /*For debugging purposes */
#endif
		CoolIntegrateEnergyCode(dm->Cool, CoolData, &cp, &E,
					ExternalHeating, fDensity,
					p->fMetals(), r, dtUse, columnL);
#else /*COOLING_MOLECULARH*/
		CoolIntegrateEnergyCode(dm->Cool, CoolData, &cp, &E,
					ExternalHeating, fDensity,
					p->fMetals(), r, dtUse);
#endif /*COOLING_MOLECULARH*/
		CkAssert(E > 0.0);
		if(dtUse > 0 || ExternalHeating*duDelta[p->rung] + p->u() < 0)
		    // linear interpolation over interval
		    p->uDot() = (E - p->u())/duDelta[p->rung];
		if (bUpdateState) p->CoolParticle() = cp;
		}
	    else { 
		p->uDot() = ExternalHeating;
		}
	    }
	}
#endif
    // Use shadow array to avoid reduction conflict
    smoothProxy[thisIndex].ckLocal()->contribute(cb);
    }

/* Set a maximum ball for inverse Nearest Neighbor searching */
void TreePiece::ballMax(int activeRung, double dhFac, const CkCallback& cb)
{
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	if (TYPETest(&myParticles[i], TYPE_GAS)) {
	    myParticles[i].fBallMax() = myParticles[i].fBall*dhFac;
	    }
	}
    // Use shadow array to avoid reduction conflict
    smoothProxy[thisIndex].ckLocal()->contribute(cb);
    }
    
int DenDvDxSmoothParams::isSmoothActive(GravityParticle *p) 
{
    if(bActiveOnly && p->rung < activeRung)
	return 0;		// not active

    return (TYPETest(p, iType));
    }

// Non-active neighbors of Actives
int DenDvDxNeighborSmParams::isSmoothActive(GravityParticle *p) 
{
    if(p->rung < activeRung && TYPETest(p, iType)
       && TYPETest(p, TYPE_NbrOfACTIVE))
	return 1;

    return 0;
    }

// Only do actives

int MarkSmoothParams::isSmoothActive(GravityParticle *p) 
{
    if(p->rung < activeRung)
	return 0;		// not active

    return (TYPETest(p, iType));
    }

/// A remote neighbor particle is active.
void MarkSmoothParams::combSmoothCache(GravityParticle *p1,
                                       ExternalSmoothParticle *p2)
{
    p1->iType |= p2->iType;
}

void DenDvDxSmoothParams::initSmoothParticle(GravityParticle *p)
{
    TYPEReset(p, TYPE_NbrOfACTIVE);
    }

void DenDvDxSmoothParams::initTreeParticle(GravityParticle *p)
{
    TYPEReset(p, TYPE_NbrOfACTIVE);
    }


void DenDvDxSmoothParams::postTreeParticle(GravityParticle *p)
{
#ifdef CULLENALPHA
    if(p->isGas())
	p->dvds_old() = p->dvdsOnSFull();
#endif
}


void DenDvDxSmoothParams::initSmoothCache(GravityParticle *p)
{
	}

void DenDvDxSmoothParams::combSmoothCache(GravityParticle *p1,
					  ExternalSmoothParticle *p2)
{
	p1->iType |= p2->iType;
	}

/* Gather only version */
void DenDvDxSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth,
				    pqSmoothNode *nnList)
{
  double ih2,ih, r2,rs,rs1,fDensity,fNorm,fNorm1,vFac;
	double dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz;
	double dvx,dvy,dvz,dx,dy,dz,trace,grx,gry,grz;
#ifdef SUPERBUBBLE
    double fDensityU = 0;
#endif
#ifdef CULLENALPHA
	double R_CD, R_CDN;     ///< R in CD limiter, and
                                ///  normalization for R.
	double maxVSignal;      ///< Maximum signal velocity
        R_CD = 0.0; R_CDN = 0;  maxVSignal = 0.0;
#endif
        double divvnorm = 0.0;
	GravityParticle *q;
	int i;
	unsigned int qiActive;

        ih2 = invH2(p); 
        ih = sqrt(ih2); 
	vFac = 1./(a*a); /* converts v to xdot */
	fNorm = M_1_PI*ih2*sqrt(ih2);
	fDensity = 0.0;
	dvxdx = 0; dvxdy = 0; dvxdz= 0;
	dvydx = 0; dvydy = 0; dvydz= 0;
	dvzdx = 0; dvzdy = 0; dvzdz= 0;
	grx = 0; gry = 0; grz= 0;

	qiActive = 0;
	for (i=0;i<nSmooth;++i) {
		double fDist2 = nnList[i].fKey;
		r2 = fDist2*ih2;
		q = nnList[i].p;
		if(q == NULL)
		    CkAbort("NULL neighbor in DenDvDxSmooth");
		if (p->rung >= activeRung)
		    TYPESet(q,TYPE_NbrOfACTIVE); /* important for SPH */
		if(q->rung >= activeRung)
		    qiActive = 1;
		rs = KERNEL(r2, nSmooth);
		fDensity += rs*q->mass;
#ifdef SUPERBUBBLE
		fDensityU += rs*q->mass*q->uPred();
#endif
		rs1 = DKERNEL(r2);
		rs1 *= q->mass;
		dx = nnList[i].dx.x; /* NB: dx = px - qx */
		dy = nnList[i].dx.y;
		dz = nnList[i].dx.z;
		dvx = (-p->vPred().x + q->vPred().x)*vFac;
		dvy = (-p->vPred().y + q->vPred().y)*vFac;
		dvz = (-p->vPred().z + q->vPred().z)*vFac;
		dvxdx += dvx*dx*rs1;
		dvxdy += dvx*dy*rs1;
		dvxdz += dvx*dz*rs1;
		dvydx += dvy*dx*rs1;
		dvydy += dvy*dy*rs1;
		dvydz += dvy*dz*rs1;
		dvzdx += dvz*dx*rs1;
		dvzdy += dvz*dy*rs1;
		dvzdz += dvz*dz*rs1;
                divvnorm += (dx*dx+dy*dy+dz*dz)*rs1;
                /* Grad P estimate */
		/* This used to be:
                   grx += (-p->uPred + q->uPred)*dx*rs1; But that is
                   rho grad u*/
                grx += (q->uPred())*dx*rs1;
                gry += (q->uPred())*dy*rs1;
                grz += (q->uPred())*dz*rs1;

#ifdef CULLENALPHA
                // Special weighting function to reduce noise in R
                // calculation.  See discussion after eq. 29 in
                // Wadsley et al 2017.
                double R_wt = (1-r2*r2*0.0625)* q->mass;
                R_CD += q->dvds_old() * R_wt; 
                R_CDN += R_wt;
                // Convention here dvdx = vxq-vxp, dx = xp-xq so
                // dvdotdr = -dvx*dx ...
                double dvdotdr = -(dvx*dx + dvy*dy + dvz*dz)
                    + fDist2*H; // vFac already in there
                double cavg = (p->c() + q->c())*0.5;
                double vSig = cavg - (dvdotdr < 0 ? dvdotdr/sqrt(fDist2) : 0);
                if (vSig > maxVSignal) maxVSignal = vSig; 

#endif
		}

        if (qiActive)
          TYPESet(p,TYPE_NbrOfACTIVE);

        p->fDensity = fNorm*fDensity;
#ifdef SUPERBUBBLE
    fDensityU *= fNorm;
    if (nSmooth > 1)
    {
        rs = KERNEL(0.0);
        p->fDensityU() = (fDensityU-rs*p->mass*p->uPred()*fNorm)/p->uPred()*p->fDensity/(p->fDensity-rs*p->mass*fNorm);
        CkAssert(p->fDensityU() > 0);
    }
    else
    {
        p->fDensityU() = fDensityU;
    }
    double ih = sqrt(ih2);
    double rhogradu=sqrt(grx*grx+gry*gry+grz*grz)*fNorm*ih2;
    p->fThermalLength() = (rhogradu != 0 ? fDensityU/rhogradu : FLT_MAX);
    if (p->fThermalLength()*ih < 1) p->fThermalLength() = 1/ih;
#endif
        trace = dvxdx+dvydy+dvzdz;
        // keep Norm positive consistent w/ std 1/rho norm
        fNorm1 = (divvnorm != 0 ? 3.0/fabs(divvnorm) : 0.0); 

#if defined(DIFFUSION) || defined(CULLENALPHA)
      double onethirdtrace = (1./3.)*trace;
      /* Build Traceless Strain Tensor (not yet normalized) */
      double sxx = dvxdx - onethirdtrace; /* pure compression/expansion doesn't diffuse */
      double syy = dvydy - onethirdtrace;
      double szz = dvzdz - onethirdtrace;
      double sxy = 0.5*(dvxdy + dvydx); /* pure rotation doesn't diffuse */
      double sxz = 0.5*(dvxdz + dvzdx);
      double syz = 0.5*(dvydz + dvzdy);
#endif
		}
		
#ifdef DIFFUSION
      /* diff coeff., nu ~ C L^2 S (add C via dMetalDiffusionConstant, assume L ~ h) */
      if (bConstantDiffusion) p->diff() = 1;
      else p->diff() = fNorm1*0.25*p->fBall*p->fBall*sqrt(2*(sxx*sxx + syy*syy + szz*szz + 2*(sxy*sxy + sxz*sxz + syz*syz)));
#endif
        p->divv() =  fNorm1*trace + 3.0*H; /* physical */
        p->curlv().x = fNorm1*(dvzdy - dvydz);
        p->curlv().y = fNorm1*(dvxdz - dvzdx);
        p->curlv().z = fNorm1*(dvydx - dvxdy);

#ifdef CULLENALPHA 
        double alphaLoc, tau; 
        double l = 0.1;
        double Hcorr = (fNorm1 != 0 ? H/fNorm1 : 0);
        double gnorm = (grx*grx+gry*gry+grz*grz);
        if (gnorm > 0) gnorm=1/sqrt(gnorm);
        grx *= gnorm;
        gry *= gnorm;
        grz *= gnorm;
        double dvdr = (((dvxdx+Hcorr)*grx+dvxdy*gry+dvxdz*grz)*grx
        +  (dvydx*grx+(dvydy+Hcorr)*gry+dvydz*grz)*gry
        +  (dvzdx*grx+dvzdy*gry+(dvzdz+Hcorr)*grz)*grz)*fNorm1;

        double pdvds_old = p->dvds();
        double dvds = (p->divv() < 0 ? 1.5*(dvdr -(1./3.)*p->divv()) : dvdr );
        double sxxf = dvxdx+Hcorr, syyf = dvydy+Hcorr, szzf = dvzdz+Hcorr;
        double SFull = sqrt(fNorm1*fNorm1*(sxxf*sxxf+syyf*syyf+szzf*szzf 
                                           + 2*(sxy*sxy + sxz*sxz + syz*syz)));
       
        p->dvdsOnSFull() = SFull > 0 ? dvds/SFull : 0; 
#ifdef CD_SFULL
        p->dvds() = p->dvdsOnSFull();
#else
        p->dvds() = dvds;
#endif

        // time interval = current time - last time divv was calculated
        double deltaT = dTime - p->TimeDivV();
        double divVDot = (p->dvds() - pdvds_old)/deltaT;
        p->TimeDivV() = dTime;
        // If we are initializing the simulation, the current time step is zero and we can't compute the time
        // derivative of the velocity divergence in the Cullen & Dehnen formulation
        if (bStarting && !bHaveAlpha){
          // If the divergence of the velocity of the particle is negative and the speed of sound is nonzero
          // we set p->CullenAlpha() using the M&M prescription. Otherwise p->CullenAlpha() is zero
          if ((p->divv() < 0) && (p->c() > 0)){
            tau = p->fBall / (2.0*l*maxVSignal);
            alphaLoc = -dAlphaMax*p->divv()*tau / (1.0 - p->divv()*tau);
          }
          else alphaLoc = 0.0; 
          p->CullenAlpha() = alphaLoc;
        }
        // If the current time step > 0
        else if(!bHaveAlpha) {
          if (p->dvds() < 0 && divVDot < 0){ // Flow is converging and
                                             // convergence is increasing
            double OneMinusR_CD = (R_CDN > 0 ? 1-(R_CD/R_CDN) : 0);

            double xi = (OneMinusR_CD < -1 ? 0 :
                (OneMinusR_CD > 2 ? 1 : 0.0625*OneMinusR_CD*OneMinusR_CD*OneMinusR_CD*OneMinusR_CD));
            // Multiplier for ATerm in CD viscosity.
            const double dAFac = 2.0;
            double Aterm = xi * p->fBall * p->fBall * fabs(divVDot)*dAFac;

            // The local alpha value
            alphaLoc = dAlphaMax* Aterm / (maxVSignal*maxVSignal + Aterm);
            
          }
          else alphaLoc = 0; 
          // Decay parameter
          tau = 1. / (l*maxVSignal*ih);
          // If alphaLoc is larger then the current p->CullenAlpha(), we set p->CullenAlhpa() to be equal to alphaLoc.
          // Otherwise, we decay p->CullenAlpha() to the alphaLoc value

          if (alphaLoc > p->CullenAlpha()) p->CullenAlpha() = alphaLoc;
          else{
            double oldCullenAlpha = p->CullenAlpha();
            p->CullenAlpha() = alphaLoc - (alphaLoc - oldCullenAlpha)*exp(-deltaT/tau);
           }
        }
#endif /* CULLENALPHA */
}
	
void DenDvDxNeighborSmParams::postTreeParticle(GravityParticle *p)
{
#ifdef CULLENALPHA
    if(p->isGas())
	p->dvds_old() = p->dvdsOnSFull();
#endif
}

void 
TreePiece::sphViscosityLimiter(int bOn, int activeRung, const CkCallback& cb)
{
    int i;
    GravityParticle *p;    

    // Pressure will be called next, so check this here.
    CkAssert(bBucketsInited);
    
    if (bOn) {
        for(i=1; i<= myNumParticles; ++i) {
	    p = &myParticles[i];
	    /* Only set values for particles with fresh curlv, divv
	       from smooth */
	    if(TYPETest(p, TYPE_GAS) && p->rung >= activeRung) {
		if (p->divv() != 0.0) {         	 
		    p->BalsaraSwitch() = fabs(p->divv())/
			(fabs(p->divv()) + sqrt(p->curlv().lengthSquared()));
		    }
		else { 
		    p->BalsaraSwitch() = 0.0;
		    }
		}
	    }
        }
    else {
        for(i=1; i<= myNumParticles; ++i) {
	    p = &myParticles[i];
	    if(TYPETest(p, TYPE_GAS)) {
		p->BalsaraSwitch() = 1.0;
		}
	    }
        }
    // Use shadow array to avoid reduction conflict
    smoothProxy[thisIndex].ckLocal()->contribute(cb);
    }

/* Note: Uses uPred */
void TreePiece::getAdiabaticGasPressure(double gamma, double gammam1, double dThermalCondCoeff,
        double dThermalCond2Coeff, double dThermalCondSatCoeff, double dThermalCond2SatCoeff,
        double dEvapMinTemp, double dtFacCourant, const CkCallback &cb)
{
    GravityParticle *p;
    double PoverRho;
    int i;

    for(i=1; i<= myNumParticles; ++i) {
	p = &myParticles[i];
	if (TYPETest(p, TYPE_GAS)) {
	    PoverRho = gammam1*p->uPred();
	    p->PoverRho2() = PoverRho/p->fDensity;
	    p->c() = sqrt(gamma*PoverRho);
#ifdef SUPERBUBBLE
        // Include the hot phase of two phase particles
        double frac = p->massHot()/p->mass;
        PoverRho = gammam1*(p->uHotPred()*frac+p->uPred()*(1-frac));
        p->c() = sqrt(gamma*PoverRho);
        p->PoverRho2() = PoverRho/p->fDensity;
#ifndef COOLING_NONE
        double Tp = CoolCodeEnergyToTemperature(dm->Cool, &p->CoolParticle(), p->uPred(), p->fMetals());
        double fThermalCond = dThermalCondCoeff*pow(p->uPred(),2.5);
        double fThermalCond2 = dThermalCond2Coeff*pow(p->uPred(),0.5);
        if (Tp < dEvapMinTemp)
        {
            fThermalCond = 0;
            fThermalCond2 = 0;
        }
        double fSat = p->fDensity*p->c()*p->fThermalLength();
        double fThermalCondSat = fSat*dThermalCondSatCoeff;
        double fThermalCond2Sat = fSat*dThermalCond2SatCoeff;
        p->fThermalCond() = (fThermalCond < fThermalCondSat ? fThermalCond : fThermalCondSat) +
            (fThermalCond2 < fThermalCond2Sat ? fThermalCond2 : fThermalCond2Sat);
        assert(isfinite(p->fThermalCond()));
#endif
#endif
#ifdef DTADJUST
            {
                double uDot = p->PdV();
                double dt;
                if(uDot > 0.0)
                    dt = dtFacCourant*0.5*p->fBall
                        /sqrt(4.0*(p->c()*p->c() + GAMMA_NONCOOL*uDot*p->dt));
                else
                    dt = dtFacCourant*0.5*p->fBall /(2.0*p->c());
                // Update to scare the neighbors.
                if(dt < p->dtNew()) p->dtNew() = dt;
                }
#endif
	    }
	}
    // Use shadow array to avoid reduction conflict
    smoothProxy[thisIndex].ckLocal()->contribute(cb);
    }

/* Note: Uses uPred */
void TreePiece::getCoolingGasPressure(double gamma, double gammam1, double dThermalCondCoeff,
        double dThermalCond2Coeff, double dThermalCondSatCoeff, double dThermalCond2SatCoeff,
        double dEvapMinTemp, double dtFacCourant, double dResolveJeans, const CkCallback &cb)
{
#ifndef COOLING_NONE
    GravityParticle *p;
    double PoverRho;
    int i;
    COOL *cl = dm->Cool;

    for(i=1; i<= myNumParticles; ++i) {
	p = &myParticles[i];
	if (TYPETest(p, TYPE_GAS)) {
            CkAssert(p->uPred() < LIGHTSPEED*LIGHTSPEED/cl->dErgPerGmUnit);
            double cGas;
	    CoolCodePressureOnDensitySoundSpeed(cl, &p->CoolParticle(),
						p->uPred(), p->fDensity(),
						gamma, gammam1, &PoverRho,
						&cGas);
            double dPoverRhoJeans = PoverRhoFloorJeans(dResolveJeans, p);
            if(PoverRho < dPoverRhoJeans) PoverRho = dPoverRhoJeans;
	    p->PoverRho2() = PoverRho/p->fDensity;
        p->c() = sqrt(cGas*cGas + GAMMA_JEANS*dPoverRhoJeans);
#ifdef SUPERBUBBLE
        double frac = p->massHot()/p->mass;
        PoverRho = gammam1*(p->uHotPred()*frac+p->uPred()*(1-frac));
        p->c() = sqrt(gamma*PoverRho);
        p->PoverRho2() = PoverRho/p->fDensity;
        double fThermalCond = dThermalCondCoeff*pow(p->uPred(),2.5);
        double fThermalCond2 = dThermalCond2Coeff*pow(p->uPred(),0.5);
        double Tp = CoolCodeEnergyToTemperature(cl, &p->CoolParticle(), p->uPred(), p->fMetals());
        if (Tp < dEvapMinTemp) // Only allow conduction & evaporation for particles that are hot
        {
            fThermalCond = 0;
            fThermalCond2 = 0;
        }
        double fSat = p->fDensity*p->c()*p->fThermalLength();
        double fThermalCondSat = fSat*dThermalCondSatCoeff;
        double fThermalCond2Sat = fSat*dThermalCond2SatCoeff;
        p->fThermalCond() = (fThermalCond < fThermalCondSat ? fThermalCond : fThermalCondSat) +
            (fThermalCond2 < fThermalCond2Sat ? fThermalCond2 : fThermalCond2Sat);
        assert(isfinite(p->fThermalCond()));
#endif
#ifdef DTADJUST
            {
                double uDot = p->uDot();
                double dt;
#ifdef SUPERBUBBLE
                uDot = p->uHotDot()*frac + p->uDot()*(1.0-frac);
#endif
                if(uDot > 0.0)
                    dt = dtFacCourant*0.5*p->fBall
                        /sqrt(4.0*(p->c()*p->c() + GAMMA_NONCOOL*uDot*p->dt));
                else
                    dt = dtFacCourant*0.5*p->fBall /(2.0*p->c());
                // Update to scare the neighbors.
                if(dt < p->dtNew()) p->dtNew() = dt;
                }
#endif
	    }
	}
#endif
    // Use shadow array to avoid reduction conflict
    smoothProxy[thisIndex].ckLocal()->contribute(cb);
    }

int PressureSmoothParams::isSmoothActive(GravityParticle *p) 
{
    return (TYPETest(p, TYPE_NbrOfACTIVE));
    }

/* Original Particle */
void PressureSmoothParams::initSmoothParticle(GravityParticle *p)
{
	if (p->rung >= activeRung) {
	    p->mumax() = 0.0;
#ifdef DTADJUST
            p->dtNew() = FLT_MAX;
#endif
	    p->PdV() = 0.0;
#ifdef DIFFUSION
	    p->fMetalsDot() = 0.0;
	    p->fMFracOxygenDot() = 0.0;
	    p->fMFracIronDot() = 0.0;
#endif /* DIFFUSION */
	    }
	}

/* Cached copies of particle */
void PressureSmoothParams::initSmoothCache(GravityParticle *p)
{
	if (p->rung >= activeRung) {
	    p->mumax() = 0.0;
#ifdef DTADJUST
            p->dtNew() = FLT_MAX;
#endif
	    p->PdV() = 0.0;
	    p->treeAcceleration = 0.0;
#ifdef DIFFUSION
	    p->fMetalsDot() = 0.0;
	    p->fMFracOxygenDot() = 0.0;
	    p->fMFracIronDot() = 0.0;
#endif /* DIFFUSION */
	    }
	}

void PressureSmoothParams::combSmoothCache(GravityParticle *p1,
					  ExternalSmoothParticle *p2)
{
	if (p1->rung >= activeRung) {
	    p1->PdV() += p2->PdV;
	    if (p2->mumax > p1->mumax())
		p1->mumax() = p2->mumax;
	    p1->treeAcceleration += p2->treeAcceleration;
#ifdef DIFFUSION
	    p1->fMetalsDot() += p2->fMetalsDot;
	    p1->fMFracOxygenDot() += p2->fMFracOxygenDot;
	    p1->fMFracIronDot() += p2->fMFracIronDot;
#endif /* DIFFUSION */
	    }
#ifdef DTADJUST
        // All neighbors get their rungs adjusted.
        if (p2->dtNew < p1->dtNew())
            p1->dtNew() = p2->dtNew;
#endif
	}

void PressureSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth,
                    pqSmoothNode *nnList)
{
    GravityParticle *q;
    PressSmoothUpdate params;
    PressSmoothParticle pParams;
    PressSmoothParticle qParams;
    double ih2,r2,rs1;
    Vector3D<double> dv;
    double ph,absmu;
    double fNorm1,vFac;
    double fDivv_Corrector;
    double dt;
    int i;

    if(nSmooth < 2) {
        CkError("WARNING: lonely SPH particle\n");
        return;
    }
#ifndef RTFORCE
    pParams.PoverRho2 = p->PoverRho2();
    pParams.PoverRho2f = pParams.PoverRho2;
#endif
    ph = 0.5 * p->fBall;
    ih2 = invH2(p);
    fNorm1 = 0.5*M_1_PI*ih2*ih2/ph;	/* converts to physical u */
    params.aFac = a;        /* comoving acceleration factor */
    vFac = 1./(a*a); /* converts v to xdot */

#ifdef RTFORCE
    double divvi = 0;
    double divvj = 0;
    for (i=0;i<nSmooth;++i) {
        double fDist2 = nnList[i].fKey;
        r2 = fDist2*ih2;
        q = nnList[i].p;
        rs1 = DKERNEL(r2);
        rs1 *= fDist2*q->mass;
        divvi += rs1;
        divvj += rs1/q->fDensity;
    }
    divvi /= p->fDensity;
    fDivv_Corrector = (divvj != 0.0 ? divvi/divvj : 1.0);
#else
    fDivv_Corrector = 1.0;
#endif

    for (i=0;i<nSmooth;++i) {
        q = nnList[i].p;
        if ((p->rung < activeRung) && (q->rung < activeRung)) continue;
        double fDist2 = nnList[i].fKey;
        r2 = fDist2*ih2;
        rs1 = DKERNEL(r2);
        rs1 *= fNorm1;
        rs1 *= fDivv_Corrector;
        pParams.rNorm = rs1 * p->mass;
        qParams.rNorm = rs1 * q->mass;
        params.dx = nnList[i].dx;
        dv = p->vPred() - q->vPred();
        params.dvdotdr = vFac*dot(dv, params.dx) + fDist2*H;
#ifdef RTFORCE
        pParams.PoverRho2 = p->PoverRho2()*p->fDensity/q->fDensity;
        pParams.PoverRho2f = pParams.PoverRho2;
        qParams.PoverRho2 = q->PoverRho2()*q->fDensity/p->fDensity;
        qParams.PoverRho2f = qParams.PoverRho2;
#else
        qParams.PoverRho2 = q->PoverRho2();
        qParams.PoverRho2f = qParams.PoverRho2;
#endif
        /***********************************
         * SPH Pressure Terms Calculation
         ***********************************/
        /* Calculate Artificial viscosity term prefactor terms 
         * 
         * Updates:
         *  dt
         *  params.visc
         */
        { // Begin SPH pressure terms calculation and scope the variables below
        if (params.dvdotdr>=0.0) {
            dt = dtFacCourant*ph/(2*(p->c() > q->c() ? p->c() : q->c()));
            params.visc = 0.0;
        } else {
            #ifdef VSIGVISC /* compile-time flag */
            /* mu multiply by a to be consistent with physical c */
            absmu = -params.dvdotdr*a/sqrt(fDist2);
            /* mu terms for gas time step */
            if (absmu>p->mumax()) p->mumax()=absmu;
            if (absmu>q->mumax()) q->mumax()=absmu;
            /* viscosity terms */
            params.visc = (varAlpha(alpha, p, q)*(p->c() + q->c())
                + varBeta(beta, p, q)*1.5*absmu);
            dt = dtFacCourant*ph/(0.625*(p->c() + q->c())+0.375*params.visc);
            params.visc = switchCombine(p,q)*params.visc*absmu/(p->fDensity + q->fDensity);
            #else
                /* h mean */
                double hav=0.5*(ph+0.5*q->fBall);
                /* mu multiply by a to be consistent with physical c */
                absmu = -hav*params.dvdotdr*a/(fDist2+0.01*hav*hav);
                /* mu terms for gas time step */
                if (absmu>p->mumax()) p->mumax()=absmu;
                if (absmu>q->mumax()) q->mumax()=absmu;
                /* viscosity terms */
                params.visc = (varAlpha(alpha, p, q)*(p->c() + q->c()) \
                        + varBeta(beta, p, q)*2*absmu);
                dt = dtFacCourant*hav/(0.625*(p->c() + q->c())+0.375*params.visc);
                params.visc = switchCombine(p,q)*params.visc*absmu/(p->fDensity + q->fDensity);
            #endif //VSIGVISC
        }
        /* Calculate diffusion terms */
        #ifdef DIFFUSION /* compile-time flag */
            // Diffusion Base term
            #ifdef DIFFUSIONHARMONIC /* compile-time flag */
                double diffSum = (p->diff()+q->diff());
                double diffBase = (diffusionLimitTest(diffSum, dTime, p, q) ? 0
                                    : 4*p->diff()*q->diff()/diffSum);
            #else
                double diffSum = (p->diff()+q->diff());
                double diffBase = (diffusionLimitTest(diffSum, dTime, p, q)
                                   ? 0 : diffSum);
            #endif
            // Metals Base term
            /* massdiff not implemented */
//            #ifdef MASSDIFF /* compile-time flag */
//                double diffMetalsBase = 4*smf->dMetalDiffusionCoeff*diffBase \
//                     /((p->fDensity+q->fDensity)*(p->fMass+q->fMass));
//            #else
                double diffMetalsBase = 2*dMetalDiffusionCoeff*diffBase \
                     /(p->fDensity+q->fDensity); 
//            #endif //MASSDIFF
        
            // Thermal diffusion
            /* 
             * Updates:
             *  dt
             *  params.diffu
             *  diffTh
             *  params.diffuNc
             */
            double diffTh;
//            /* DIFFUSIONPRICE not implemented */
//            #ifdef DIFFUSIONPRICE /* compile-time flag */
//                {
//                double irhobar = 2/(p->fDensity+q->fDensity);
//                double vsig = sqrt(fabs(qParams.PoverRho2*q->fDensity*q->fDensity \
//                                        - pParams.PoverRho2*p->fDensity*p->fDensity)\
//                                        *irhobar);
//                diffTh = smf->dThermalDiffusionCoeff*0.5 \
//                        * (ph+sqrt(0.25*BALL2(q)))*irhobar*vsig;
//                params.diffu = diffTh*(p->uPred-q->uPred);
//                }
//            #else
                #ifndef NODIFFUSIONTHERMAL /* compile-time flag */
                    {
                    diffTh = (2*dThermalDiffusionCoeff*diffBase/(p->fDensity+q->fDensity));
                    double dt_diff;
                    double dThermalCond;
//                    /* THERMALCOND not implemented */
//                    #ifdef THERMALCOND /* compile-time flag */
//                        #if (0)
//                            /* Harmonic average coeff */
//                            double dThermalCondSum = p->fThermalCond + q->fThermalCond;
//                            dThermalCond = ( dThermalCondSum <= 0 ? 0 \
//                                : 4*p->fThermalCond*q->fThermalCond \
//                                /(dThermalCondSum*p->fDensity*q->fDensity) );
//                        #else
//                            /* Arithmetic average coeff */
//                            dThermalCond = (p->fThermalCond + q->fThermalCond) \
//                                    /(p->fDensity*q->fDensity);
//                            if (dThermalCond > 0 && (dt_diff = dtFacDiffusion*ph \
//                                    *p->fThermalLength/(dThermalCond*p->fDensity)) < dt){
//                                dt = dt_diff;
//                            }
//                        #endif
//                    #else
                        dThermalCond = 0.0;
//                    #endif //THERMALCOND
                    if (diffTh > 0 && (dt_diff= dtFacDiffusion*ph*ph/(diffTh*p->fDensity)) < dt) dt = dt_diff;
                    params.diffu = (diffTh+dThermalCond)*(p->uPred()-q->uPred());
                    }
                #endif
//            #endif //DIFFUSIONPRICE
//            /* not implemented */
//            #ifdef UNONCOOL /* compile-time flag */
//                params.diffuNc = diffTh*(p->uNoncoolPred-q->uNoncoolPred);
//            #endif
            // Calculate diffusion pre-factor terms (required for updating particles)
            params.diffMetals = diffMetalsBase*(p->fMetals() - q->fMetals());
            params.diffMetalsOxygen = diffMetalsBase*(p->fMFracOxygen() - q->fMFracOxygen());
            params.diffMetalsIron = diffMetalsBase*(p->fMFracIron() - q->fMFracIron());
//            /* not implemented */
//            #ifdef MASSDIFF /* compile-time flag */
//                params.diffMass = diffMetalsBase*(p->fMass - q->fMass);
//                // To properly implement this in ChaNGa the correct velocity 
//                // should be chosen
//                params.diffVelocity = diffMetalsBase * (p->velocity - q->velocity);
//            #endif
        #endif
        if (p->rung >= activeRung) {
            updateParticle(p, q, &params, &pParams, &qParams, 1);
        }
        if (q->rung >= activeRung) {
            updateParticle(q, p, &params, &qParams, &pParams, -1);
        }
        // Adust dt
        #ifdef DTADJUST /* compile-time flag */
            if (dt < p->dtNew()) p->dtNew() = dt;
            if (dt < q->dtNew()) q->dtNew() = dt;
            if (4*q->dt < p->dtNew()) p->dtNew() = 4*q->dt;
            if (4*p->dt < q->dtNew()) q->dtNew() = 4*p->dt;
        #endif
        } // End SPH Pressure Terms calculations
    }
}

/**
 * @brief updateParticle is used to update particle attributes during the 
 * SPH pressure terms calculations.  
 * 
 * The updating of particle p and the neighbor q during this loop is symmetric
 * (up to a possible sign change).  For example, to update p and its neighbor
 * q is two lines:
 *      updateParticle(p, q, params, pParams, qParams, 1);
 *      updateParticle(q, p, params, qParams, pParams, -1);
 * @param a particle to update
 * @param b interacting neighbor particle
 * @param params prefactor params
 * @param aParams params specific to a
 * @param bParams params specific to b
 * @param sign 1 for a = p (the self particle) and -1 for a = q (the neighbor)
 */
void updateParticle(GravityParticle *a, GravityParticle *b, 
                    PressSmoothUpdate *params, PressSmoothParticle *aParams, 
                    PressSmoothParticle *bParams, int sign) {
    double acc;
    // Update diffusion terms
    #ifdef DIFFUSION /* compile-time flag */
        // Thermal diffusion
//        /* not implemented */
//        #ifdef DIFFUSIONPRICE /* compile-time flag */
//            a->uDotDiff += sign * params->diffu * bParams->rNorm;
//        #else
            #ifndef NODIFFUSIONTHERMAL /* compile-time flag */
                a->PdV() += sign * params->diffu * bParams->rNorm \
                        * massDiffFac(b);
            #endif
//        #endif //DIFFUSIONPRICE
//        /* not implemented */
//        #ifdef UNONCOOL /* compile-time flag */
//            a->uNoncoolDotDiff += sign * params->diffuNc * bParams->rNorm;
//        #endif
        // Metals diffusion
        a->fMetalsDot() += sign * params->diffMetals * bParams->rNorm
                * massDiffFac(b);
        a->fMFracOxygenDot() += sign * params->diffMetalsOxygen
                * bParams->rNorm * massDiffFac(b);
        a->fMFracIronDot() += sign * params->diffMetalsIron * bParams->rNorm
                * massDiffFac(b);
//        /* not implemented */
//        #ifdef MASSDIFF /* compile-time flag */
//        // Note: to implement this in ChaNGa, ACCEL should be properly vectorized
//            a->fMassDot += sign * params->diffMass * a->fMass * bParams->rNorm;
//            ACCEL(a) += sign * params->diffVelocity * bParams->rNorm \
//                    * massDiffFac(b);
//        #endif
    #endif
    // Update pressure/viscosity terms
//    /* not implemented */
//    #ifdef DRHODT /* compile-time flag */
//        a->fDivv_PdV -= bParams->rNorm / params->fDivv_Corrector / \
//                rhoDivv(a->fDensity,b->fDensity) * params->dvdotdr;
//        a->fDivv_PdVcorr -= bParams->rNorm / rhoDivv(a->fDensity,b->fDensity) \
//                * params->dvdotdr;
//    #endif
    a->PdV() += bParams->rNorm*presPdv(aParams->PoverRho2, bParams->PoverRho2)
            * params->dvdotdr;
    a->PdV() += bParams->rNorm * 0.5 * params->visc * params->dvdotdr;
    acc = presAcc(aParams->PoverRho2f, bParams->PoverRho2f) \
            + params->visc;
    acc *= bParams->rNorm * params->aFac;
    a->treeAcceleration -= sign * acc * params->dx;
}

/*
 * Methods to distribute Deleted gas
 */
int DistDeletedGasSmoothParams::isSmoothActive(GravityParticle *p) 
{
    return (TYPETest(p, TYPE_DELETED) && TYPETest(p, iType));
    }

void DistDeletedGasSmoothParams::initSmoothCache(GravityParticle *p) 
{
    if(!TYPETest(p, TYPE_DELETED)) {
	/*
	 * Zero out accumulated quantities.
	 */
	p->mass = 0;
	p->velocity[0] = 0;
	p->velocity[1] = 0;
	p->velocity[2] = 0;
#ifndef COOLING_NONE
	p->u() = 0;
	p->uDot() = 0.0;
#endif
	p->fMetals() = 0.0;
	p->fMFracIron() = 0.0;
	p->fMFracOxygen() = 0.0;
#ifdef SUPERBUBBLE
	p->massHot() = 0;
	p->uHot() = 0;
	p->uHotPred() = 0;
#endif
	}
    }

void DistDeletedGasSmoothParams::combSmoothCache(GravityParticle *p1,
					  ExternalSmoothParticle *p2)
{
    /*
     * Distribute u, v, and fMetals for particles returning from cache
     * so that everything is conserved nicely.  
     */
    if(!TYPETest((p1), TYPE_DELETED)) {
	double delta_m = p2->mass;
	double m_new,f1,f2;
	double fTCool; /* time to cool to zero */
	m_new = p1->mass + delta_m;
	if (delta_m > 0) {
	    f1 = p1->mass /m_new;
	    f2 = delta_m  /m_new;
	    p1->velocity = f1*p1->velocity + f2*p2->velocity;            
	    p1->fMetals() = f1*p1->fMetals() + f2*p2->fMetals;
	    p1->fMFracIron() = f1*p1->fMFracIron() + f2*p2->fMFracIron;
	    p1->fMFracOxygen() = f1*p1->fMFracOxygen() + f2*p2->fMFracOxygen;
#ifdef SUPERBUBBLE
        double mHot_new = p1->massHot() + p2->massHot;
        if (mHot_new > 0) {
            double f1_hot = p1->massHot()/mHot_new;
            double f2_hot = p2->massHot/mHot_new;
            double mCold_new = m_new-mHot_new;
            assert(mCold_new>0);
            double f1_cold = (p1->mass-p1->massHot())/mCold_new;
            double f2_cold = (delta_m-p2->massHot)/mCold_new;
            p1->uHot() = f1_hot*p1->uHot()+f2_hot*p2->uHot;
            p1->uHotPred() = f1_hot*p1->uHotPred()+f2_hot*p2->uHotPred;
            p1->u() = f1_cold*p1->u()+f2_cold*p2->u;
            p1->uPred() = f1_cold*p1->uPred()+f2_cold*p2->uPred;
            p1->massHot() = mHot_new;
        }
        else
#endif
#ifndef COOLING_NONE
        {
            if(p1->uDot() < 0.0) /* margin of 1% to avoid roundoff
                      * problems */
            fTCool = 1.01*p1->uPred()/p1->uDot(); 
            p1->u() = f1*p1->u() + f2*p2->u;
            p1->uPred() = f1*p1->uPred() + f2*p2->uPred;
            if(p1->uDot() < 0.0)
            p1->uDot() = p1->uPred()/fTCool;
        }
#endif
	    p1->mass = m_new;
	    }
	}
    }

void DistDeletedGasSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth,
					   pqSmoothNode *nnList)
{	
    GravityParticle *q;
    double fNorm,ih2,r2,rs,rstot,delta_m,m_new,f1,f2;
    double fTCool; /* time to cool to zero */
    int i;
    CkAssert(TYPETest(p, TYPE_GAS));
    ih2 = invH2(p);
    rstot = 0;        
    for (i=0;i<nSmooth;++i) {
	double fDist2 = nnList[i].fKey;
	q = nnList[i].p;
	if(TYPETest(q, TYPE_DELETED)) continue;
	CkAssert(TYPETest(q, TYPE_GAS));
	r2 = fDist2*ih2;            
        rs = KERNEL(r2, nSmooth);
	rstot += rs;
        }
    if(rstot <= 0.0) {
	if(p->mass == 0.0) /* the particle to be deleted has NOTHING */
	    return;
	/* we have a particle to delete and nowhere to put its mass
	 * => we will keep it around */
	unDeleteParticle(p);
	return;
	}
    CkAssert(rstot > 0.0);
    fNorm = 1./rstot;
    CkAssert(p->mass >= 0.0);
#ifdef SUPERBUBBLE
    if (p->massHot() > 0) {
        pqSmoothNode *massList = (pqSmoothNode *) malloc(sizeof(pqSmoothNode)*nSmooth);
        for (i=0;i<nSmooth;++i) {
            massList[i].p = nnList[i].p;
            massList[i].fKey = -1*nnList[i].p->massHot();
        }
        std::sort_heap(massList, massList+nSmooth);
        for(i=0;i<nSmooth;++i) {
            q = massList[i].p;
            m_new = q->mass + p->mass;
            /* Cached copies can have zero mass: skip them */
            if (m_new == 0) continue;
            double mHot_new = q->massHot()+p->massHot();
            f1 = q->mass/m_new;
            f2 = p->mass/m_new;
            double mCold_new = m_new-mHot_new;
            assert(mCold_new > 0);
            double f1_cold = (q->mass-q->massHot())/mCold_new;
            double f2_cold = (p->mass-p->massHot())/mCold_new;
            double f1_hot = q->massHot()/mHot_new;
            double f2_hot = p->massHot()/mHot_new;
            q->velocity = f1*q->velocity + f2*p->velocity;            
            q->fMetals() = f1*q->fMetals() + f2*p->fMetals();
            q->fMFracIron() = f1*q->fMFracIron() + f2*p->fMFracIron();
            q->fMFracOxygen() = f1*q->fMFracOxygen() + f2*p->fMFracOxygen();
            q->u() = f1_cold*q->u()+f2_cold*p->u();
            q->uPred() = f1_cold*q->uPred()+f2_cold*p->uPred();
            q->uHot() = f1_hot*q->uHot()+f2_hot*p->uHot();
            q->uHotPred() = f1_hot*q->uHotPred()+f2_hot*p->uHotPred();
            q->mass = m_new;
            q->massHot() = mHot_new;
            return;
        }
        free(massList);
    }
#endif
    for (i=0;i<nSmooth;++i) {
	q = nnList[i].p;
	if(TYPETest(q, TYPE_DELETED)) continue;

	double fDist2 = nnList[i].fKey;
	r2 = fDist2*ih2;            
	rs = KERNEL(r2, nSmooth);
	/*
	 * All these quantities are per unit mass.
	 * Exact if only one gas particle being distributed or in serial
	 * Approximate in parallel (small error).
	 */
	delta_m = rs*fNorm*p->mass;
	m_new = q->mass + delta_m;
	/* Cached copies can have zero mass: skip them */
	if (m_new == 0) continue;
	f1 = q->mass /m_new;
	f2 = delta_m  /m_new;
	q->mass = m_new;
	q->velocity = f1*q->velocity + f2*p->velocity;            
	q->fMetals() = f1*q->fMetals() + f2*p->fMetals();
	q->fMFracIron() = f1*q->fMFracIron() + f2*p->fMFracIron();
	q->fMFracOxygen() = f1*q->fMFracOxygen() + f2*p->fMFracOxygen();
#ifndef COOLING_NONE
	if(q->uDot() < 0.0) /* margin of 1% to avoid roundoff error */
	    fTCool = 1.01*q->uPred()/q->uDot(); 
	q->u() = f1*q->u()+f2*p->u();
	q->uPred() = f1*q->uPred()+f2*p->uPred();
	if(q->uDot() < 0.0) /* make sure we don't shorten cooling time */
	    q->uDot() = q->uPred()/fTCool;
#endif
        }
    }

#ifdef SUPERBUBBLE
int PromoteToHotGasSmoothParams::isSmoothActive(GravityParticle *p)
{
    return TYPETest(p, TYPE_FEEDBACK) && TYPETest(p, iType) && !TYPETest(p, TYPE_PROMOTED);
}
void PromoteToHotGasSmoothParams::initSmoothParticle(GravityParticle *p)
{
    /* Initialize the promotion sums */
    TYPEReset(p,TYPE_PROMOTED);
    p->fPromoteSum() = 0;
    p->fPromoteSumuPred() = 0;
    p->fPromoteuPredInit() = p->uPred();

}
void PromoteToHotGasSmoothParams::initTreeParticle(GravityParticle *p)
{
    TYPEReset(p,TYPE_PROMOTED);
    p->fPromoteSum() = 0;
    p->fPromoteSumuPred() = 0;
    p->fPromoteuPredInit() = p->uPred();
}
void PromoteToHotGasSmoothParams::initSmoothCache(GravityParticle *p)
{
    TYPEReset(p,TYPE_PROMOTED);
    p->fPromoteSum() = 0;
    p->fPromoteSumuPred() = 0;
    p->fPromoteuPredInit() = p->uPred();
}
void PromoteToHotGasSmoothParams::combSmoothCache(GravityParticle *p1, ExternalSmoothParticle *p2)
{
    if(TYPETest(p2, TYPE_PROMOTED)) {
        TYPESet(p1,TYPE_PROMOTED);
        if (p2->fTimeCoolIsOffUntil > p1->fTimeCoolIsOffUntil()) p1->fTimeCoolIsOffUntil() = p2->fTimeCoolIsOffUntil;
        }
    p1->fPromoteSum() += p2->fPromoteSum;
    p1->fPromoteSumuPred() += p2->fPromoteSumuPred;
}
void PromoteToHotGasSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth,
        pqSmoothNode *nnList)
{	
#ifdef NOCOOLING
    return;
#endif
    GravityParticle *q;
    double fFactor,ph,ih2,r2,rs,rstot;
    double Tp,Tq,up52,Prob,mPromoted;
    double xc,yc,zc,dotcut2,dot;
	int i,nCold,nHot;

	CkAssert(TYPETest(p, TYPE_GAS));
	CkAssert(TYPETest(p, TYPE_FEEDBACK));
	CkAssert(!TYPETest(p, TYPE_PROMOTED));
    ph = 0.5*p->fBall;
    ih2 = invH2(p);
    /* Exclude cool particles */
    Tp = CoolEnergyToTemperature(tp->Cool(), &p->CoolParticle(), dErgPerGmUnit*p->uPred(), p->fMetals() );
    CkAssert(Tp < 2e11);
    if (Tp <= dEvapMinTemp) return;

    up52 = pow(p->uPred(),2.5);
    rstot = 0;
    xc = 0; yc = 0; zc = 0; 
    nCold = 0;
	for (i=0;i<nSmooth;++i) {
        q = nnList[i].p;
        if (p->iOrder == q->iOrder) continue;
	    if (TYPETest(q, TYPE_DELETED) || (TYPETest(q, TYPE_FEEDBACK) && !TYPETest(q, TYPE_PROMOTED))) continue;
        Tq = CoolEnergyToTemperature(tp->Cool(), &q->CoolParticle(), dErgPerGmUnit*q->uPred(), q->fMetals() );
        CkAssert(Tq < 2e11);
        if (q->uHot() > 0 || Tq >= dEvapMinTemp) continue;  /* Exclude hot particles */
	    CkAssert(TYPETest(q, TYPE_GAS));
        CkAssert(!TYPETest(p, TYPE_STAR));
		r2 = nnList[i].fKey*ih2;            
		rs = KERNEL(r2);
        rstot += rs;
  		xc += rs*nnList[i].dx.x; 
		yc += rs*nnList[i].dx.y;
		zc += rs*nnList[i].dx.z;
        nCold++;
        }

    if (rstot == 0) return;

    /* Check for non-edge hot particle  theta = 45 deg, cos^2 = 0.5 */
    dotcut2 = (xc*xc+yc*yc+zc*zc)*0.5;
    
	for (i=0;i<nSmooth;++i) {
		q = nnList[i].p;
		if (p->iOrder == q->iOrder) continue;
		if (TYPETest(q, TYPE_DELETED)) continue;
        Tq = CoolEnergyToTemperature(tp->Cool(), &q->CoolParticle(), dErgPerGmUnit*q->uPred(), q->fMetals() );
        CkAssert(Tq < 2e11);
        if (q->uHot() == 0 && Tq <= dEvapMinTemp) continue;  
		dot = xc*nnList[i].dx.x + yc*nnList[i].dx.y + zc*nnList[i].dx.y;
		if (dot > 0 && dot*dot > dotcut2*nnList[i].fKey) {
            return;
            }
        }

    /* Area = h^2 4 pi nCold/nSmooth */
	nHot=nSmooth-nCold;
	CkAssert(nHot > 0);
    fFactor = dDeltaStarForm*dEvapCoeff*ph*12.5664*1.5/(nHot)/rstot;

    mPromoted = 0;
	for (i=0;i<nSmooth;++i) {
        q = nnList[i].p;
        if (p->iOrder == q->iOrder) continue;
	    if(TYPETest(q, TYPE_DELETED) || (TYPETest(q, TYPE_FEEDBACK) && !TYPETest(q, TYPE_PROMOTED))) continue;
        Tq = CoolEnergyToTemperature(tp->Cool(), &q->CoolParticle(), dErgPerGmUnit*q->uPred(), q->fMetals() );
        CkAssert(Tq < 2e11);
        if (Tq >= dEvapMinTemp ) continue;  /* Exclude hot particles */
	    CkAssert(TYPETest(q, TYPE_GAS));
		r2 = nnList[i].fKey*ih2;            
		rs = KERNEL(r2);
        q->fPromoteSum() += p->mass;
        q->fPromoteSumuPred() += p->mass*p->uPred();
		
        /* cf. Weaver etal'77 mdot = 4.13d-14 * (dx^2/4 !pi) (Thot^2.5-Tcold^2.5)/dx - 2 udot mHot/(k T/mu) 
           Kernel sets total probability to 1 */
        Prob = fFactor*(up52-pow(q->uPred(),2.5))*rs/q->mass;
        if ( (rand()/((double) RAND_MAX)) < Prob) {
            mPromoted += q->mass; 
            }
        }

    if (mPromoted > 0) {
        double dTimeCool = dTime + 0.9999*dDeltaStarForm;

        std::sort_heap(nnList, nnList+nSmooth);
        for (i=0;i<nSmooth;++i) {
            q = nnList[i].p;
            if (p->iOrder == q->iOrder) continue;
            if (TYPETest(q, TYPE_DELETED) || TYPETest(q, TYPE_FEEDBACK) || TYPETest(q, TYPE_PROMOTED)) continue;
            Tq = CoolEnergyToTemperature(tp->Cool(), &q->CoolParticle(), dErgPerGmUnit*q->uPred(), q->fMetals() );
            CkAssert(Tq < 2e11);
            if (Tq >= dEvapMinTemp ) continue;  /* Exclude hot particles */
            CkAssert(TYPETest(q, TYPE_GAS));

            if (dTimeCool > q->fTimeCoolIsOffUntil()) q->fTimeCoolIsOffUntil() = dTimeCool;
            TYPESet(q, TYPE_PROMOTED|TYPE_FEEDBACK);
            mPromoted -= q->mass;
            if (mPromoted < q->mass*0.1) break;
            }
        }
}

int ShareWithHotGasSmoothParams::isSmoothActive(GravityParticle *p)
{
    return TYPETest(p, TYPE_FEEDBACK) && TYPETest(p, iType) && !TYPETest(p, TYPE_PROMOTED);
}
        
void ShareWithHotGasSmoothParams::initSmoothCache(GravityParticle *p)
{
    if(!TYPETest(p, TYPE_DELETED)) {
        p->u() = 0;  
        p->uPred() = 0;
        }
}

void ShareWithHotGasSmoothParams::combSmoothCache(GravityParticle *p1, ExternalSmoothParticle *p2)
{
	if(!TYPETest((p1), TYPE_DELETED)) {
        p1->u() +=  p2->u;
        p1->uPred() += p2->uPred;
		}
}

void ShareWithHotGasSmoothParams::fcnSmooth(GravityParticle *p,int nSmooth,
        pqSmoothNode *nnList)
{
	GravityParticle *q;
	double uavg,umin;
	double dE,Eadd,factor,Tp;
	int i,nPromoted;

	CkAssert(TYPETest(p, TYPE_GAS));
	CkAssert(TYPETest(p, TYPE_FEEDBACK));
	CkAssert(!TYPETest(p, TYPE_PROMOTED));
    Tp = CoolEnergyToTemperature(tp->Cool(), &p->CoolParticle(), dErgPerGmUnit*p->uPred(), p->fMetals() );
    if (Tp <= dEvapMinTemp) return;

    nPromoted = 0;

    dE = 0;
    umin = FLT_MAX;
	for (i=0;i<nSmooth;++i) {
        q = nnList[i].p;
	    if (TYPETest(q, TYPE_PROMOTED)) {
            nPromoted++;
            uavg = (q->mass*q->fPromoteuPredInit() + q->fPromoteSumuPred())/
                (q->mass + q->fPromoteSum());
            if (uavg < umin) umin=uavg;
            Eadd = (uavg-q->fPromoteuPredInit())*q->mass;
            if (Eadd < 0) Eadd=0;
            dE += p->mass/q->fPromoteSum()*Eadd;
            }
        }

    if (!nPromoted || dE == 0 || p->uPred() <= umin) return;
    factor = ((p->uPred()-umin)*p->mass)/dE;
    if (factor > 1) factor=1;

	for (i=0;i<nSmooth;++i) {
        q = nnList[i].p;
	    if (TYPETest(q, TYPE_PROMOTED)) {
            nPromoted++;
            uavg = (q->mass*q->fPromoteuPredInit() + q->fPromoteSumuPred())/
                (q->mass + q->fPromoteSum());
            if (uavg < umin) umin=uavg;
            Eadd = (uavg-q->fPromoteuPredInit())*q->mass;
            if (Eadd < 0) continue; //Stop evaporating once we have Eadd worth of mass
            dE = factor*p->mass/q->fPromoteSum()*Eadd;
            q->uPred() += dE/q->mass;
            q->u() += dE/q->mass;
            p->uPred() -=  dE/p->mass;
            p->u() -=  dE/p->mass;
            if (!(q->uPred() > 0))
            {
                CkPrintf("SHARE ERROR: %e %e %e %e %e %e\n", q->uPred(), dE, factor, q->fPromoteSum(), q->fPromoteSumuPred(), q->fPromoteuPredInit());
            }
            CkAssert(q->uPred() > 0);
            CkAssert(q->u() > 0);
            CkAssert(p->uPred() > 0);
            CkAssert(p->u() > 0);
            CkAssert(p->u() < LIGHTSPEED*LIGHTSPEED/tp->Cool()->dErgPerGmUnit);
            CkAssert(p->uPred() < LIGHTSPEED*LIGHTSPEED/tp->Cool()->dErgPerGmUnit);
            CkAssert(q->u() < LIGHTSPEED*LIGHTSPEED/tp->Cool()->dErgPerGmUnit);
            CkAssert(q->uPred() < LIGHTSPEED*LIGHTSPEED/tp->Cool()->dErgPerGmUnit);
            }
        }
}
#endif
