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
	// actives, and those who have actives as neighbors.
	DenDvDxSmoothParams pDen(TYPE_GAS, 0, param.csm, dTime, 0,
				 param.bConstantDiffusion);
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
			     CkCallbackResumeThread());
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
	    double T,E;
#ifndef COOLING_NONE
#ifndef COOLING_GRACKLE
	    T = p->u() / dTuFac;
            PERBARYON Y;
            CoolPARTICLEtoPERBARYON(cl, &Y, &p->CoolParticle());
            
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
				 param.bConstantDiffusion);
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
				      param.bConstantDiffusion);
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
				 param.bConstantDiffusion);
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
    if(param.bGasCooling)
	treeProxy.getCoolingGasPressure(param.dConstGamma,
					param.dConstGamma-1,
                                        param.dResolveJeans/csmTime2Exp(param.csm, dTime),
					CkCallbackResumeThread());
    else
	treeProxy.getAdiabaticGasPressure(param.dConstGamma,
					  param.dConstGamma-1,
					  CkCallbackResumeThread());

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
	    double T,E;
#ifndef COOLING_NONE
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
			   const CkCallback& cb)
{
    double dt; // time in seconds
    
#ifndef COOLING_NONE
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
	if (TYPETest(p, TYPE_GAS)
	    && (p->rung == activeRung || (bAll && p->rung >= activeRung))) {
	    dt = CoolCodeTimeToSeconds(dm->Cool, duDelta[p->rung] );
	    double ExternalHeating = p->PdV();
	    ExternalHeating += p->fESNrate();
	    if ( bCool ) {
		COOLPARTICLE cp = p->CoolParticle();
		double E = p->u();
		double r[3];  // For conversion to C
		p->position.array_form(r);
		double dtUse = dt;
		
		if(dStartTime[p->rung] + 0.5*duDelta[p->rung]
		   < p->fTimeCoolIsOffUntil()) {
		    /* This flags cooling shutoff (e.g., from SNe) to
		       the cooling functions. */
		    dtUse = -dt;
		    p->uDot() = ExternalHeating;
		    }

		CoolIntegrateEnergyCode(dm->Cool, CoolData, &cp, &E,
					ExternalHeating, p->fDensity,
					p->fMetals(), r, dtUse);
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

void DenDvDxSmoothParams::initSmoothParticle(GravityParticle *p)
{
    TYPEReset(p, TYPE_NbrOfACTIVE);
    }

void DenDvDxSmoothParams::initTreeParticle(GravityParticle *p)
{
    TYPEReset(p, TYPE_NbrOfACTIVE);
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
	double ih2,r2,rs,rs1,fDensity,fNorm,fNorm1,vFac;
	double dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz;
	double dvx,dvy,dvz,dx,dy,dz,trace;
        double divvnorm = 0.0;
	GravityParticle *q;
	int i;
	unsigned int qiActive;

	ih2 = invH2(p);
	vFac = 1./(a*a); /* converts v to xdot */
	fNorm = M_1_PI*ih2*sqrt(ih2);
	fDensity = 0.0;
	dvxdx = 0; dvxdy = 0; dvxdz= 0;
	dvydx = 0; dvydy = 0; dvydz= 0;
	dvzdx = 0; dvzdy = 0; dvzdz= 0;

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
		rs = KERNEL(r2);
		fDensity += rs*q->mass;
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
		}
	if (qiActive)
	    TYPESet(p,TYPE_NbrOfACTIVE);
		
	p->fDensity = fNorm*fDensity; 
	trace = dvxdx+dvydy+dvzdz;
        /* keep Norm positive consistent w/ std 1/rho norm */
        fNorm1 = (divvnorm != 0 ? 3.0/fabs(divvnorm) : 0.0);
	p->divv() =  fNorm1*trace + 3.0*H; /* physical */
	p->curlv().x = fNorm1*(dvzdy - dvydz); 
	p->curlv().y = fNorm1*(dvxdz - dvzdx);
	p->curlv().z = fNorm1*(dvydx - dvxdy);
#ifdef DIFFUSION
        {
	double onethirdtrace = (1./3.)*trace;
	/* Build Traceless Strain Tensor (not yet normalized) */
	double sxx = dvxdx - onethirdtrace; /* pure compression/expansion doesn't diffuse */
	double syy = dvydy - onethirdtrace;
	double szz = dvzdz - onethirdtrace;
	double sxy = 0.5*(dvxdy + dvydx); /* pure rotation doesn't diffuse */
	double sxz = 0.5*(dvxdz + dvzdx);
	double syz = 0.5*(dvydz + dvzdy);
	/* diff coeff., nu ~ C L^2 S (add C via dMetalDiffusionConstant, assume L ~ h) */
	if (bConstantDiffusion) p->diff() = 1;
	else p->diff() = fNorm1*0.25*p->fBall*p->fBall*sqrt(2*(sxx*sxx + syy*syy + szz*szz + 2*(sxy*sxy + sxz*sxz + syz*syz)));
	}
#endif
	}

/* As above, but no marking */
void DenDvDxNeighborSmParams::fcnSmooth(GravityParticle *p, int nSmooth,
				    pqSmoothNode *nnList)
{
	double ih2,r2,rs,rs1,fDensity,fNorm,fNorm1,vFac;
	double dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz;
	double dvx,dvy,dvz,dx,dy,dz,trace;
	GravityParticle *q;
	int i;

	ih2 = invH2(p);
	vFac = 1./(a*a); /* converts v to xdot */
	fNorm = M_1_PI*ih2*sqrt(ih2);
	fNorm1 = fNorm*ih2;	
	fDensity = 0.0;
	dvxdx = 0; dvxdy = 0; dvxdz= 0;
	dvydx = 0; dvydy = 0; dvydz= 0;
	dvzdx = 0; dvzdy = 0; dvzdz= 0;

	for (i=0;i<nSmooth;++i) {
		double fDist2 = nnList[i].fKey;
		r2 = fDist2*ih2;
		q = nnList[i].p;
		rs = KERNEL(r2);
		fDensity += rs*q->mass;
		rs1 = DKERNEL(r2);
		rs1 *= q->mass;
		dx = nnList[i].dx.x;
		dy = nnList[i].dx.y;
		dz = nnList[i].dx.z;
		dvx = (-p->vPred().x + q->vPred().x)*vFac - dx*H; /* NB: dx = px - qx */
		dvy = (-p->vPred().y + q->vPred().y)*vFac - dy*H;
		dvz = (-p->vPred().z + q->vPred().z)*vFac - dz*H;
		dvxdx += dvx*dx*rs1;
		dvxdy += dvx*dy*rs1;
		dvxdz += dvx*dz*rs1;
		dvydx += dvy*dx*rs1;
		dvydy += dvy*dy*rs1;
		dvydz += dvy*dz*rs1;
		dvzdx += dvz*dx*rs1;
		dvzdy += dvz*dy*rs1;
		dvzdz += dvz*dz*rs1;
		}
		
	p->fDensity = fNorm*fDensity; 
	fNorm1 /= p->fDensity;
	trace = dvxdx+dvydy+dvzdz;
	p->divv() =  fNorm1*trace; /* physical */
	p->curlv().x = fNorm1*(dvzdy - dvydz); 
	p->curlv().y = fNorm1*(dvxdz - dvzdx);
	p->curlv().z = fNorm1*(dvydx - dvxdy);
#ifdef DIFFUSION
        {
	double onethirdtrace = (1./3.)*trace;
	/* Build Traceless Strain Tensor (not yet normalized) */
	double sxx = dvxdx - onethirdtrace; /* pure compression/expansion doesn't diffuse */
	double syy = dvydy - onethirdtrace;
	double szz = dvzdz - onethirdtrace;
	double sxy = 0.5*(dvxdy + dvydx); /* pure rotation doesn't diffuse */
	double sxz = 0.5*(dvxdz + dvzdx);
	double syz = 0.5*(dvydz + dvzdy);
	/* diff coeff., nu ~ C L^2 S (add C via dMetalDiffusionConstant, assume L ~ h) */
	if (bConstantDiffusion) p->diff() = 1;
	else p->diff() = fNorm1*0.25*p->fBall*p->fBall*sqrt(2*(sxx*sxx + syy*syy + szz*szz + 2*(sxy*sxy + sxz*sxz + syz*syz)));
	}
#endif
	}

void 
TreePiece::sphViscosityLimiter(int bOn, int activeRung, const CkCallback& cb)
{
    int i;
    GravityParticle *p;    

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
void TreePiece::getAdiabaticGasPressure(double gamma, double gammam1,
					const CkCallback &cb)
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
	    }
	}
    // Use shadow array to avoid reduction conflict
    smoothProxy[thisIndex].ckLocal()->contribute(cb);
    }

/* Note: Uses uPred */
void TreePiece::getCoolingGasPressure(double gamma, double gammam1,
                                      double dResolveJeans,
                                      const CkCallback &cb)
{
    GravityParticle *p;
    double PoverRho;
    int i;
#ifndef COOLING_NONE
    COOL *cl = dm->Cool;

    for(i=1; i<= myNumParticles; ++i) {
	p = &myParticles[i];
	if (TYPETest(p, TYPE_GAS)) {
            double cGas;
	    CoolCodePressureOnDensitySoundSpeed(cl, &p->CoolParticle(),
						p->uPred(), p->fDensity(),
						gamma, gammam1, &PoverRho,
						&cGas);
            double dPoverRhoJeans = PoverRhoFloorJeans(dResolveJeans, p);
            if(PoverRho < dPoverRhoJeans) PoverRho = dPoverRhoJeans;
	    p->PoverRho2() = PoverRho/p->fDensity;
            p->c() = sqrt(cGas*cGas + GAMMA_JEANS*dPoverRhoJeans);
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
#ifdef DTADJUST
            if (p2->dtNew < p1->dtNew())
                p1->dtNew() = p2->dtNew;
#endif
	    p1->treeAcceleration += p2->treeAcceleration;
#ifdef DIFFUSION
	    p1->fMetalsDot() += p2->fMetalsDot;
	    p1->fMFracOxygenDot() += p2->fMFracOxygenDot;
	    p1->fMFracIronDot() += p2->fMFracIronDot;
#endif /* DIFFUSION */
	    }
	}

void PressureSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth,
				    pqSmoothNode *nnList)
{
	GravityParticle *q;
	double ih2,r2,rs1,rq,rp;
	double dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	double pPoverRho2,pPoverRho2f,pMass;
	double qPoverRho2,qPoverRho2f;
	double ph,pc,pDensity,visc,hav,absmu,Accp,Accq;
	double fNorm,fNorm1,aFac,vFac, divvi, divvj;
	double fDivv_Corrector;
	double dt;
	int i;

	if(nSmooth < 2) {
	    CkError("WARNING: lonely SPH particle\n");
	    return;
	    }
	pc = p->c();
	pDensity = p->fDensity;
	pMass = p->mass;
#ifndef RTFORCE
	pPoverRho2 = p->PoverRho2();
	pPoverRho2f = pPoverRho2;
#endif
	ph = sqrt(0.25*p->fBall*p->fBall);
	ih2 = invH2(p);
	fNorm = 0.5*M_1_PI*ih2/ph;
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = a;        /* comoving acceleration factor */
	vFac = 1./(a*a); /* converts v to xdot */

	divvi = 0;
	divvj = 0;
	for (i=0;i<nSmooth;++i) {
	    double fDist2 = nnList[i].fKey;
	    r2 = fDist2*ih2;
	    q = nnList[i].p;
	    rs1 = DKERNEL(r2);
	    rs1 *= fDist2*q->mass;
	    divvi += rs1/p->fDensity;
	    divvj += rs1/q->fDensity;
	    }
#ifdef RTFORCE
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
	    rp = rs1 * pMass;
	    rq = rs1 * q->mass;

	    dx = nnList[i].dx.x;
	    dy = nnList[i].dx.y;
	    dz = nnList[i].dx.z;
	    dvx = p->vPred()[0] - q->vPred()[0];
	    dvy = p->vPred()[1] - q->vPred()[1];
	    dvz = p->vPred()[2] - q->vPred()[2];
	    dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) + fDist2*H;
#ifdef RTFORCE
	    pPoverRho2 = p->PoverRho2()*pDensity/q->fDensity;
	    pPoverRho2f = pPoverRho2;
	    qPoverRho2 = q->PoverRho2()*q->fDensity/pDensity;
	    qPoverRho2f = qPoverRho2;
#else
	    qPoverRho2 = q->PoverRho2();
	    qPoverRho2f = qPoverRho2;
#endif

	    if (p->rung >= activeRung) {
		if (q->rung >= activeRung) {
#define PACTIVE(xxx) xxx
#define QACTIVE(xxx) xxx
#include "SphPressureTerms.h"
		    }
		else {
#undef QACTIVE
#define QACTIVE(xxx) 
#include "SphPressureTerms.h"
		    }
		}
	    else if (q->rung >= activeRung) {
#undef PACTIVE
#define PACTIVE(xxx) 
#undef QACTIVE
#define QACTIVE(xxx) xxx
#include "SphPressureTerms.h"
		}
	    }
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
	if (m_new > 0) {
	    f1 = p1->mass /m_new;
	    f2 = delta_m  /m_new;
	    p1->mass = m_new;
	    p1->velocity = f1*p1->velocity + f2*p2->velocity;            
	    p1->fMetals() = f1*p1->fMetals() + f2*p2->fMetals;
	    p1->fMFracIron() = f1*p1->fMFracIron() + f2*p2->fMFracIron;
	    p1->fMFracOxygen() = f1*p1->fMFracOxygen() + f2*p2->fMFracOxygen;
#ifndef COOLING_NONE
	    if(p1->uDot() < 0.0) /* margin of 1% to avoid roundoff
				  * problems */
		fTCool = 1.01*p1->uPred()/p1->uDot(); 
	    p1->u() = f1*p1->u() + f2*p2->u;
	    p1->uPred() = f1*p1->uPred() + f2*p2->uPred;
	    if(p1->uDot() < 0.0)
		p1->uDot() = p1->uPred()/fTCool;
#endif
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
        rs = KERNEL(r2);
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
    for (i=0;i<nSmooth;++i) {
	q = nnList[i].p;
	if(TYPETest(q, TYPE_DELETED)) continue;

	double fDist2 = nnList[i].fKey;
	r2 = fDist2*ih2;            
	rs = KERNEL(r2);
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
