/*
 * Image dumping routines for movies from PKDGRAV
 * Original author: James Wadsley, 2002
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <errno.h>
#include "config.h"
#ifdef HAVE_VALUES_H
#include <values.h>
#else
#include <float.h>
#endif

#include <charm.h>

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157E+308
#endif

#include "dumpframe.h"

void dfInitialize( struct DumpFrameContext **pdf, double dYearUnit, double dTime,  
				  double dDumpFrameTime, double dStep, double dDumpFrameStep,
				  double dDelta, int iMaxRung, int bVDetails, char* filename ) {
	double tock=0.0;/*initialized to suppress warning: DCR 12/19/02*/
	struct DumpFrameContext *df;

	if (*pdf!=NULL) (*pdf)->fs = NULL;
	if (dDumpFrameTime <= 0 && dDumpFrameStep <=0) return;

	if (*pdf == NULL) {
		*pdf = (struct DumpFrameContext *) malloc(sizeof(struct DumpFrameContext ));
		(*pdf)->bAllocated = 1;
		}
	else {
		(*pdf)->bAllocated = 0;
		}

	df = *pdf;
	df->iMaxRung = 0;
	df->dYearUnit = dYearUnit;

	df->dTime = 0;
	df->dDumpFrameTime = 0;
	if (dDumpFrameTime > 0) { 
		df->dDumpFrameTime = dDumpFrameTime;
		df->dTime = dTime-dDumpFrameTime*1e-10;
		/* Round to a nearby power of two */
		tock = -log((dDumpFrameTime-floor(dDumpFrameTime))/dDelta*0.50001+1e-60)
			/log(2.0);
		if (tock <= iMaxRung) {
			df->iMaxRung = tock;
			}
		}

	df->dStep = 0;
	df->dDumpFrameStep = 0;
	if (dDumpFrameStep > 0) { 
		df->dDumpFrameStep = dDumpFrameStep;
		df->dStep = dStep-dDumpFrameStep*1e-10;
		/* Round to a nearby power of two */
		tock = -log((dDumpFrameStep-floor(dDumpFrameStep))*0.50001+1e-60)
			/log(2.0);
		if (tock <= iMaxRung && tock >df->iMaxRung) {
			df->iMaxRung = floor(tock);
			}
		}
/* General options */
	df->bVDetails = bVDetails;
	dfParseOptions( *pdf, filename );

/* Parse time dependent camera options (mostly 2D) */
	dfParseCameraDirections( *pdf, filename );

	if (df->bVDetails) CkPrintf("DF Initialized Frame dumping: Time Interval %g [%g] (Step %i) Step Interval %g -> MaxRung %i\n",dDumpFrameTime,floor(dDumpFrameTime)+dDelta*pow(2.0,-floor(tock)),(int) floor(dDumpFrameTime/dDelta),dDumpFrameStep,df->iMaxRung);
	}


void dfFinalize( struct DumpFrameContext *df ) {
	if (df != NULL) {
		if (df->fs != NULL) free( df->fs );
		if (df->bAllocated) free( df );
		}
	}

void *dfAllocateImage( int nxPix, int nyPix ) {
	DFIMAGE *Image;

	assert ( nxPix*nyPix*sizeof(DFIMAGE) <= DF_NBYTEDUMPFRAME );

	Image = (DFIMAGE *) malloc( nxPix*nyPix*sizeof(DFIMAGE) );
	assert (Image != NULL );
	

	return Image;
	}

void dfFreeImage( void *Image ) {
    free( Image );
	}


#define SET( a, b ) { \
							   a[0]=b[0]; \
							   a[1]=b[1]; \
						       a[2]=b[2]; \
} 
#define ADD( a, b ) { \
							   a[0]+=b[0]; \
							   a[1]+=b[1]; \
						       a[2]+=b[2]; \
} 
#define DIFF( a, b, c ) { \
							   c[0] = b[0]-a[0]; \
							   c[1] = b[1]-a[1]; \
							   c[2] = b[2]-a[2]; \
						  }
#define DOT( a, b ) ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] )
#define DIST( a, b ) sqrt( (a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]) )
#define LEN( a ) sqrt( (a[0])*(a[0])+(a[1])*(a[1])+(a[2])*(a[2]) )
#define NORM( a ) { \
						double ilen; \
						ilen = 1./LEN ( a ); \
						a[0] *= ilen; a[1] *= ilen; a[2] *= ilen; \
					}
#define SIZEVEC( a, f ) { \
						double ilen; \
						ilen = f/LEN ( a ); \
						a[0] *= ilen; a[1] *= ilen; a[2] *= ilen; \
					}
#define CURL( a, b, c ) { \
							  c[0] = a[1]*b[2] - a[2]*b[1]; \
							  c[1] = a[2]*b[0] - a[0]*b[2]; \
							  c[2] = a[0]*b[1] - a[1]*b[0]; \
						  }

void dfProjection( struct inDumpFrame *in, struct dfFrameSetup *fs, int nxPix, int nyPix) {
	double width,height;
	double vec[3];
	/*double norm[3]; -- not used: DCR 12/19/02*/

	in->r[0] = fs->target[0];
	in->r[1] = fs->target[1];
	in->r[2] = fs->target[2];

	if ((nxPix == 0) && (nyPix == 0)) {
	  in->nxPix = fs->nxPix;
	  in->nyPix = fs->nyPix;
	}
	else {
	  in->nxPix = nxPix;
	  in->nyPix = nyPix;
	}

	in->bExpansion = fs->bExpansion;
	in->bPeriodic = fs->bPeriodic;
	in->iProject = fs->iProject;
	in->pScale1 = fs->pScale1;
	in->pScale2 = fs->pScale2;
	in->ColStar = fs->ColStar;
	in->ColGas = fs->ColGas;
	in->ColDark = fs->ColDark;
	in->bColMassWeight = fs->bColMassWeight;
	in->bGasSph = fs->bGasSph;
	in->iColStarAge = fs->iColStarAge;
	in->iLogScale = fs->iLogScale;
	in->iTarget = fs->iTarget;
	in->dGasSoftMul = fs->dGasSoftMul;
	in->dDarkSoftMul = fs->dDarkSoftMul;
	in->dStarSoftMul = fs->dStarSoftMul;
	in->iRender = fs->iRender;
	
	if (fs->bEyeAbsolute) {
		DIFF( fs->eye, in->r, in->z );
		}
	else {
		double zero[3]={ 0, 0, 0 };
		DIFF( fs->eye, zero, in->z );
		}

/*	fprintf(stderr,"eye2: %i  %lf %lf %lf,  %lf %lf %lf",fs->bEye2,fs->eye2[0],fs->eye2[1],fs->eye2[2], in->z[0],in->z[1],in->z[2] ); */

	if (fs->bzEye1) SIZEVEC( in->z, fs->zEye1 );

/* 	fprintf(stderr,"eye2: %i  %lf %lf %lf,  %lf %lf %lf",fs->bEye2,fs->eye2[0],fs->eye2[1],fs->eye2[2], in->z[0],in->z[1],in->z[2] ); */


	if (fs->bEye2) {
		SET( vec, fs->eye2 );
		if (fs->bzEye2) SIZEVEC( vec, fs->zEye2 );
		ADD ( in->z, vec );
		}

/*	fprintf(stderr,"eye2: %i  %lf %lf %lf,  %lf %lf %lf",fs->bEye2,fs->eye2[0],fs->eye2[1],fs->eye2[2], in->z[0],in->z[1],in->z[2] ); */

	if (fs->bzEye) {
		in->zEye = fs->zEye;
		}
	else {
		in->zEye = LEN( in->z );
		}

	if (fs->bExpansion) in->zEye/= in->dExp; /* Use physical units for sizing? */
	/* zEye is used to scale viewport */

	assert( in->zEye > 0.0 );
	if (fs->bzClipFrac) {
		in->zClipNear = in->zEye*fs->zClipNear;
		in->zClipFar = in->zEye*fs->zClipFar;
		}
	else {
		in->zClipNear = fs->zClipNear;
		in->zClipFar = fs->zClipFar;
		}

    NORM( in->z );

	width = 2*tan( fs->FOV*M_PI/180.*0.5 )*in->zEye;
	height = width*in->nyPix/in->nxPix;

    CURL( fs->up, in->z, in->x );
    if (fs->iProject == DF_PROJECT_PERSPECTIVE) {
		SIZEVEC( in->x, (in->nxPix*0.5*in->zEye/width) );
		}
	else {
		SIZEVEC( in->x, (in->nxPix*0.5/width) );
		}

	CURL( in->z, in->x, in->y );
    if (fs->iProject == DF_PROJECT_PERSPECTIVE) {
		SIZEVEC( in->y, (in->nyPix*0.5*in->zEye/height) );
		}
	else {
		SIZEVEC( in->y, (in->nyPix*0.5/height) );
		}

	if (fs->bPeriodic) {
		in->nxRepNeg = 0;  /* Replicas for Periodic: Will setup sensibly ultimately */
		in->nxRepPos = 0;
		in->nyRepNeg = 0;  
		in->nyRepPos = 0;
		in->nzRepNeg = 0;  
		in->nzRepPos = 0;
		}
	else {
		in->nxRepNeg = 0;  /* Not Periodic */
		in->nxRepPos = 0;
		in->nyRepNeg = 0;  
		in->nyRepPos = 0;
		in->nzRepNeg = 0;  
		in->nzRepPos = 0;
		}
	
	/* Render helper variables */
	in->xlim = (in->nxPix-1)*.5;
	in->ylim = (in->nyPix-1)*.5;
	in->hmul = 4*sqrt(in->x[0]*in->x[0] + in->x[1]*in->x[1] + in->x[2]*in->x[2]);

	/* in->bNonLocal not set (an internal use only variable) */

	if (in->bVDetails) {
		CkPrintf("DF Projection: %i x %i FOV %f  width %f height %f\n",fs->nxPix,fs->nyPix,fs->FOV,width,height);
		CkPrintf("DF eye %f %f %f, target %f %f %f, (separation %f)\n",fs->eye[0],fs->eye[1],fs->eye[2],fs->target[0],fs->target[1],fs->target[2],in->zEye );
		CkPrintf("DF up %f %f %f  Z-Clipping: Near %f Far %f\n",fs->up[0],fs->up[1],fs->up[2],in->zClipNear,in->zClipFar);
		CkPrintf("DF Vectors: x %f %f %f, y %f %f %f z %f %f %f\n",in->x[0],in->x[1],in->x[2], in->y[0],in->y[1],in->y[2], in->z[0],in->z[1],in->z[2] );
		CkPrintf("DF Colours: star %f %f %f gas %f %f %f dark %f %f %f\n",in->ColStar.r,in->ColStar.g,in->ColStar.b,in->ColGas.r,in->ColGas.g,in->ColGas.b,in->ColDark.r,in->ColDark.g,in->ColDark.b);
		}
	}

void dfParseOptions( struct DumpFrameContext *df, char * filename ) {
	FILE *fp;
	int nitem;
	char line[81],command[40],word[40];
	char FileBaseName[81]="Frame";
	char NumberingFormat[41]="";

	df->iDimension = DF_2D;
	df->iEncode = DF_ENCODE_PPM;
	df->iNumbering = DF_NUMBERING_FRAME;
	df->dMassGasMin  = 0;
	df->dMassGasMax  = DBL_MAX;
	df->dMassDarkMin = 0;
	df->dMassDarkMax = DBL_MAX;
	df->dMassStarMin = 0;
	df->dMassStarMax = DBL_MAX;
	df->nFrame = 0;
	df->bGetCentreOfMass = 0;
	df->bGetOldestStar = 0;
	df->bGetPhotogenic = 0;
	df->bVDetails = 0;
	sprintf(df->FileName, "%s.%%09i.ppm", FileBaseName);

	fp = fopen( filename, "r" );
	if (fp==NULL) return;

	for ( ;; ) {
		if (fgets( line, 81, fp ) == NULL) break;
		sscanf( line, "%s", command );
		if (!strcmp( command, "t" ) || !strcmp( command, "time" )) df->nFrameSetup++; /* count times */
		}

	rewind( fp );

	for ( ;; ) {
		if (fgets( line, 81, fp ) == NULL) break;
		nitem = sscanf( line, "%s", command );
		if (nitem != 1 || command[0]=='#') continue;
		else if (!strcmp( command, "dim") || !strcmp( command, "dimension") ) {
			nitem = sscanf( line, "%s %s", command, word );
			assert( nitem == 2 );
		    if (!strcmp( word, "2") ) {
				df->iDimension = DF_2D;
			}
		    else if (!strcmp( word, "3") ) {
				df->iDimension = DF_3D;
			}
			else {
				fprintf(stderr,"DF Unknown dimensions: %s\n",word);
				assert( 0 );
				} 
			}
		else if (!strcmp( command, "file") ) {
			nitem = sscanf( line, "%s %80s", command, FileBaseName );
			assert( strlen(FileBaseName) > 0 );
			}
		else if (!strcmp( command, "gas") ) {
			nitem = sscanf( line, "%s %s", command, word );
			if (nitem == 2 && (!strcmp( word, "no") ||!strcmp( word, "off" ))) {
				df->dMassGasMin = DBL_MAX;
				}
			else {
				nitem = sscanf( line, "%s %lf %lf", command,
							   &df->dMassGasMin, &df->dMassGasMax );
				assert( nitem == 3 );
				}
			}
		else if (!strcmp( command, "dark") ) {
			nitem = sscanf( line, "%s %s", command, word );
			if (nitem == 2 && (!strcmp( word, "no") ||!strcmp( word, "off" ))) {
				df->dMassDarkMin = DBL_MAX;
				}
			else {
				nitem = sscanf( line, "%s %lf %lf", command,
							   &df->dMassDarkMin, &df->dMassDarkMax );
				assert( nitem == 3 );
				}
			}
		else if (!strcmp( command, "star") ) {
			nitem = sscanf( line, "%s %s", command, word );
			if (nitem == 2 && (!strcmp( word, "no") ||!strcmp( word, "off" ))) {
				df->dMassStarMin = DBL_MAX;
				}
			else {
				nitem = sscanf( line, "%s %lf %lf", command,
							   &df->dMassStarMin, &df->dMassStarMax );
				assert( nitem == 3 );
				}
			}
		else if (!strcmp( command, "encode") ) {
			nitem = sscanf( line, "%s %s", command, word );
			assert( nitem == 2 );
		    if (!strcmp( word, "ppm") ) {
				df->iEncode = DF_ENCODE_PPM;
			}
		    else if (!strcmp( word, "png") ) {
				df->iEncode = DF_ENCODE_PNG;
#ifndef USE_PNG
				fprintf(stderr,"DF PNG encoding support not compiled in\n");
				assert(0);
#endif				
			}
		    else if (!strcmp( word, "rle") ) {
				df->iEncode = DF_ENCODE_RLE;
				fprintf(stderr,"DF RLE encoding not supported yet\n");
				assert(0);
			}
		    else if (!strcmp( word, "treezip") ) {
				df->iEncode = DF_ENCODE_TREEZIP;
				fprintf(stderr,"DF Voxel encoding not supported yet\n");
				assert(0);
			}
			else {
				fprintf(stderr,"DF Unknown encoding: %s\n",word);
				assert(0);
				} 
			}
		else if (!strcmp( command, "numbering") ) {
			nitem = sscanf( line, "%s %s %40s", command, word, NumberingFormat );
			if (nitem == 2 ) {
				nitem = sscanf( line, "%s %s", command, word );
				assert( nitem == 2 );
				}

		    if (!strcmp( word, "frame") ) {
				df->iNumbering = DF_NUMBERING_FRAME;
			}
		    else if (!strcmp( word, "step") ) {
				df->iNumbering = DF_NUMBERING_STEP;
			}
		    else if (!strcmp( word, "time") ) {
				df->iNumbering = DF_NUMBERING_TIME;
			}
			else {
				fprintf(stderr,"DF Unknown numbering type: %s\n",word);
				assert( 0 );
				} 
			}
		else if (!strcmp( command, "frame") ) {
			nitem = sscanf( line, "%s %i", command, &df->nFrame );
			assert( nitem == 2 );
			}
		}
	
	if (df->iDimension == DF_3D) {
		fprintf(stderr,"DF Voxel projection/encoding not supported yet\n");
		assert(0);

		df->iEncode = DF_ENCODE_TREEZIP; /* Encode same thing */
		}
	else assert (df->iEncode != DF_ENCODE_TREEZIP);
	
	if (strlen(NumberingFormat)==0) {
		switch( df->iNumbering ) {
		case DF_NUMBERING_FRAME:
			strcpy(NumberingFormat,"%09i");
			break;
		case DF_NUMBERING_STEP:
		case DF_NUMBERING_TIME:
			strcpy(NumberingFormat,"%012.6f");
			break;
			}
		}
	
	switch ( df->iEncode ) {
	case DF_ENCODE_PPM:
		sprintf(df->FileName,"%s.%s.ppm", FileBaseName, NumberingFormat );
		break;
	case DF_ENCODE_PNG:
		sprintf(df->FileName,"%s.%s.png", FileBaseName, NumberingFormat );
		break;
	default:
		sprintf(df->FileName,"%s.%s.null", FileBaseName, NumberingFormat );
		break;
		}

	fclose(fp);
	}

void dfParseCameraDirections( struct DumpFrameContext *df, char * filename ) {
	FILE *fp;
	/*	struct inDumpFrame in; */
	struct dfFrameSetup fs;
	int n,nitem;
#define CMWOFF 1
#define CMWON  2
	int bCMW = 0;
	char line[81],command[40],word[40];

	df->bLoop = 0;
	df->dTimeMod = 0;
	df->dTimeLoop = 1e20;
	df->dPeriodLoop = 0;
	
/* Defaults -- most of this is for 2D only */
	fs.dTime = 0;
	fs.nxPix = 800;
	fs.nyPix = 600;
	fs.target[0] = 0; 
	fs.target[1] = 0;
	fs.target[2] = 0;
	fs.eye[0] = 0;
	fs.eye[1] = 0;
	fs.eye[2] = -.5;
    fs.up[0] = 0;
	fs.up[1] = 1;
	fs.up[2] = 0;
	fs.zEye = 0.0;
	fs.zEye1 = 0.0;
	fs.zEye2 = 0.0;
	fs.eye2[0] = 0;
	fs.eye2[1] = 0;
	fs.eye2[2] = 0;
	fs.bEye2 = 0;
	fs.bzEye = 0;
	fs.bzEye1 = 0;
	fs.bzEye2 = 0;
	fs.bEyeAbsolute = 0;
	fs.FOV = 90.;
	fs.zClipNear = 0.01;
	fs.zClipFar = 2.0;
	fs.bzClipFrac = 1;
    fs.bExpansion = 0;  /* to physical? */
    fs.bPeriodic = 0;  /* Periodic? */
	fs.iProject = DF_PROJECT_PERSPECTIVE;
	/* Render */
	fs.ColDark.r = 1.0;
	fs.ColDark.g = 0.7;
	fs.ColDark.b = 0.5;
	fs.ColGas.r = 0.5;
	fs.ColGas.g = 1.0;
	fs.ColGas.b = 0.7;
	fs.ColStar.r = 0.5;
	fs.ColStar.g = 0.7;
	fs.ColStar.b = 1.0;
	fs.bColMassWeight = 0;
	fs.bGasSph = 1;
	fs.iColStarAge = DF_STAR_AGE_BASIC;
	fs.iLogScale = DF_LOG_NULL;
	fs.iTarget = DF_TARGET_USER;
	fs.dGasSoftMul = 1;
	fs.dDarkSoftMul = 1;
	fs.dStarSoftMul = 1;
	

	fs.iRender = DF_RENDER_POINT;

	fp = fopen( filename, "r" );
	if (fp==NULL) {
		fprintf(stderr,"DF Could not open camera director file: %s\n",filename );

		df->iFrameSetup = 0;
		df->nFrameSetup = 1;
		df->fs = (struct dfFrameSetup *) malloc(sizeof(struct dfFrameSetup)*df->nFrameSetup);
		assert( df->fs != NULL );

		if (df->bVDetails) CkPrintf("DF Default Frame Setup\n" );
		/* dfProjection( &in, &fs ); */ /* Redundant now -- done again later */
		df->fs[0] = fs;
		return;
		}

	if (df->bVDetails) CkPrintf("DF Reading camera directions from: %s\n",filename );
	
	df->iFrameSetup = -1;
	df->nFrameSetup=1; /* Defaults */
	for ( ;; ) {
		if (fgets( line, 81, fp ) == NULL) break;
		sscanf( line, "%s", command );
		if (!strcmp( command, "t" ) || !strcmp( command, "time" )) df->nFrameSetup++; /* count times */
		}

	rewind( fp );
	df->fs = (struct dfFrameSetup *) malloc(sizeof(struct dfFrameSetup)*df->nFrameSetup);
	assert( df->fs != NULL );

	n=0;
	if (df->bVDetails) CkPrintf("DF Default Frame Setup\n" );
	for ( ;; ) {
		if (fgets( line, 81, fp ) == NULL) break;
		nitem = sscanf( line, "%s", command );
		if (nitem != 1 || command[0]=='#') continue;
		else if (!strcmp( command, "t" ) || !strcmp( command, "time" )) {
			/* dfProjection( &in, &fs );  */ /* Redundant now -- done again later */
			df->fs[n] = fs;
			n++; /* next setup */

			nitem = sscanf( line, "%s %lf", command, &fs.dTime );
			assert( nitem == 2 );
			
			if (df->bVDetails) CkPrintf("DF Frame Setup from time: %f\n",fs.dTime );
			}
		else if (!strcmp( command, "target") ) {
			nitem = sscanf( line, "%s %lf %lf %lf", command, &fs.target[0], &fs.target[1], &fs.target[2] );
			if (nitem != 4) {
			  nitem = sscanf( line, "%s %s", command, word );
			  assert( nitem == 2);
			  if (!strcmp( word, "gas") ) {
				fs.iTarget = DF_TARGET_COM_GAS;
 			    df->bGetCentreOfMass = 1;
  			    }
			  else if (!strcmp( word, "dark") ) {
				fs.iTarget = DF_TARGET_COM_DARK;
			    df->bGetCentreOfMass = 1;
			    }
			  else if (!strcmp( word, "star") ) {
				fs.iTarget = DF_TARGET_COM_STAR;
			    df->bGetCentreOfMass = 1;
			    }
              else if (!strcmp( word, "all") ) {
				fs.iTarget = DF_TARGET_COM_ALL;
			    df->bGetCentreOfMass = 1;
			    }
			  else if (!strcmp( word, "oldeststar" ) ) {
				fs.iTarget = DF_TARGET_OLDEST_STAR;
			    df->bGetOldestStar = 1;
			    }
			  else if (!strcmp( word, "photogenic" ) ) {
				fs.iTarget = DF_TARGET_PHOTOGENIC;
			    df->bGetPhotogenic = 1;
			    }
			  else {
				assert(0);
			    }
		      }
		    }
		else if (!strcmp( command, "eye") || !strcmp( command, "eye1") ) {
			nitem = sscanf( line, "%s %lf %lf %lf", command, &fs.eye[0], &fs.eye[1], &fs.eye[2] );
			if ( nitem != 4 ) {
				nitem = sscanf( line, "%s %s", command, word );
				assert( nitem == 2 );
				if (!strcmp( word, "rel") ) {
					fs.bEyeAbsolute = 0;
					}
				else if (!strcmp( word, "abs") ) {
					fs.bEyeAbsolute = 1;
					}
				else {
					fprintf(stderr,"DF Unknown eye offset type: %s\n",word);
					assert( 0 );
					}
				}
			}
		else if (!strcmp( command, "eye2") ) {
			nitem = sscanf( line, "%s %lf %lf %lf", command, &fs.eye2[0], &fs.eye2[1], &fs.eye2[2] );
			fs.bEye2 = 1;
			assert( nitem == 4 );
			}
		else if (!strcmp( command, "up") ) {
			nitem = sscanf( line, "%s %lf %lf %lf", command, &fs.up[0], &fs.up[1], &fs.up[2] );
			assert( nitem == 4 );
			}
		else if (!strcmp( command, "size") ) {
			nitem = sscanf( line, "%s %i %i", command, &fs.nxPix, &fs.nyPix );
			assert( nitem == 3 );
			assert( fs.nxPix*fs.nyPix <= DF_NXPIXMAX*DF_NYPIXMAX );
			}
		else if (!strcmp( command, "zeye") ) {
			fs.bzEye = 1;
			nitem = sscanf( line, "%s %lf", command, &fs.zEye );
			assert( nitem == 2 );
			}
		else if (!strcmp( command, "zeye1") ) {
			fs.bzEye1 = 1;
			nitem = sscanf( line, "%s %lf", command, &fs.zEye1 );
			assert( nitem == 2 );
			}
		else if (!strcmp( command, "zeye2") ) {
			fs.bzEye2 = 1;
			nitem = sscanf( line, "%s %lf", command, &fs.zEye2 );
			assert( nitem == 2 );
			}
		else if (!strcmp( command, "exp") || !strcmp( command, "physical") ) {
			fs.bExpansion = 1;
			}
		else if (!strcmp( command, "fov") || !strcmp( command, "FOV")  ) {
			nitem = sscanf( line, "%s %lf", command, &fs.FOV );
			assert( nitem == 2 );
			}
		else if (!strcmp( command, "loop") ) {
			nitem = sscanf( line, "%s %lf %lf", command, &df->dTimeLoop, &df->dPeriodLoop );
			assert( nitem == 3 );
			}
		else if (!strcmp( command, "clip") ) {
			df->bLoop = 1;
			nitem = sscanf( line, "%s %lf %lf", command, &fs.zClipNear, &fs.zClipFar );
			assert( nitem == 3 );
			}
		else if (!strcmp( command, "clipabs") ) {
			fs.bzClipFrac = 0;
			nitem = sscanf( line, "%s %lf %lf", command, &fs.zClipNear, &fs.zClipFar );
			assert( nitem == 3 );
			}
		else if (!strcmp( command, "clipfrac") ) {
			fs.bzClipFrac = 1;
			nitem = sscanf( line, "%s %lf %lf", command, &fs.zClipNear, &fs.zClipFar );
			assert( nitem == 3 );
			}
		else if (!strcmp( command, "project") ) {
			nitem = sscanf( line, "%s %s", command, word );
			assert( nitem == 2 );
		    if (!strcmp( word, "ortho") ) {
				fs.iProject = DF_PROJECT_ORTHO;
			}
		    else if (!strcmp( word, "perspective") ) {
				fs.iProject = DF_PROJECT_PERSPECTIVE;
			}
			else {
				fprintf(stderr,"DF Unknown projection: %s\n",word);
				assert( 0 );
				} 
			}
		else if (!strcmp( command, "colgas" )) {
			float scaler;
			nitem = sscanf( line, "%s %f %f %f %f", command, &fs.ColGas.r, &fs.ColGas.g, &fs.ColGas.b, &scaler );
			if (nitem == 5) {
				assert( !(bCMW & CMWOFF) );
				bCMW |= CMWON;
				fs.bColMassWeight = 1;
				fs.ColGas.r/=scaler;
				fs.ColGas.g/=scaler;
				fs.ColGas.b/=scaler;
				}
			else {
				nitem = sscanf( line, "%s %f %f %f", command, &fs.ColGas.r, &fs.ColGas.g, &fs.ColGas.b );
				assert( nitem == 4 );
				assert( !(bCMW & CMWON) );
				bCMW |= CMWOFF;
				}
			}
		else if (!strcmp( command, "coldark" )) {
			float scaler;
			nitem = sscanf( line, "%s %f %f %f %f", command, &fs.ColDark.r, &fs.ColDark.g, &fs.ColDark.b, &scaler );
			if (nitem == 5) {
				assert( !(bCMW & CMWOFF) );
				bCMW |= CMWON;
				fs.bColMassWeight = 1;
				fs.ColDark.r/=scaler;
				fs.ColDark.g/=scaler;
				fs.ColDark.b/=scaler;
				}
			else {
				nitem = sscanf( line, "%s %f %f %f", command, &fs.ColDark.r, &fs.ColDark.g, &fs.ColDark.b );
				assert( nitem == 4 );
				assert( !(bCMW & CMWON) );
				bCMW |= CMWOFF;
				}
			}
		else if (!strcmp( command, "colstar" ) 
				 || !strcmp( command, "colstarbri" ) 
				 || !strcmp( command, "colstarcol" )
				 || !strcmp( command, "colstarbricol" )) {

			float scaler;
			if (!strcmp( command, "colstar" ))
			  fs.iColStarAge = DF_STAR_AGE_BASIC;
			else if(!strcmp( command, "colstarbri" ))
			  fs.iColStarAge = DF_STAR_AGE_BRIGHT;
			else if(!strcmp( command, "colstarcol" ))
			  fs.iColStarAge = DF_STAR_AGE_COLOUR;
			else if(!strcmp( command, "colstarbricol" ))
			  fs.iColStarAge = DF_STAR_AGE_BRIGHT_COLOUR;

			nitem = sscanf( line, "%s %f %f %f %f", command, &fs.ColStar.r, &fs.ColStar.g, &fs.ColStar.b, &scaler );
			if (nitem == 5) {
				assert( !(bCMW & CMWOFF) );
				bCMW |= CMWON;
				fs.bColMassWeight = 1;
				fs.ColStar.r/=scaler;
				fs.ColStar.g/=scaler;
				fs.ColStar.b/=scaler;
				}
			else {
				nitem = sscanf( line, "%s %f %f %f", command, &fs.ColStar.r, &fs.ColStar.g, &fs.ColStar.b );
				assert( nitem == 4 );
				assert( !(bCMW & CMWON) );
				bCMW |= CMWOFF;
				}
			}
		else if (!strcmp( command, "colmass" )) {
			assert( !(bCMW & CMWOFF) );
			bCMW |= CMWON;
			fs.bColMassWeight = 1;
			}
		else if (!strcmp( command, "logscale" )) {
			nitem = sscanf( line, "%s %lf %lf", command, &fs.pScale1, &fs.pScale2 );
			assert( nitem == 3 );
			fs.iLogScale = DF_LOG_SATURATE;
			}
		else if (!strcmp( command, "logscalecoloursafe" )) {
			nitem = sscanf( line, "%s %lf %lf", command, &fs.pScale1, &fs.pScale2 );
			assert( nitem == 3 );
			fs.iLogScale = DF_LOG_COLOURSAFE;
			}
		else if (!strcmp( command, "softgassph" )) {
			fs.bGasSph = 1;
			}
		else if (!strcmp( command, "softgas" )) {
			fs.bGasSph = 0;
			nitem = sscanf( line, "%s %lf", command, &fs.dGasSoftMul );
			assert( nitem == 2 );
			}
		else if (!strcmp( command, "softdark" )) {
			nitem = sscanf( line, "%s %lf", command, &fs.dDarkSoftMul );
			assert( nitem == 2 );
			}
		else if (!strcmp( command, "softstar" )) {
			nitem = sscanf( line, "%s %lf", command, &fs.dStarSoftMul );
			assert( nitem == 2 );
			}
		else if (!strcmp( command, "render") ) {
			nitem = sscanf( line, "%s %s", command, word );
			assert( nitem == 2 );
		    if (!strcmp( word, "point") ) {
				fs.iRender = DF_RENDER_POINT;
			}
		    else if (!strcmp( word, "tsc") ) {
				fs.iRender = DF_RENDER_TSC;
			}
		    else if (!strcmp( word, "solid") ) {
				fs.iRender = DF_RENDER_SOLID;
			}
		    else if (!strcmp( word, "shine") ) {
				fs.iRender = DF_RENDER_SHINE;
			}
			else {
				fprintf(stderr,"DF Unknown rendering: %s\n",word);
				assert( 0 );
				} 
			}
		}
	
	/* Redundant now -- done again later */
	/*	
	in.bVDetails = df->bVDetails;
	dfProjection( &in, &fs );
	*/ 

	df->fs[n] = fs;
	n++; /* next setup */
	if (n==1) df->iFrameSetup = 0;
	assert ( n == df->nFrameSetup );

	fclose(fp);
	}

void dfSetupFrame( struct DumpFrameContext *df, double dTime, double dStep, double dExp, double *com, struct inDumpFrame *vin, int nxPix, int nyPix ) {
	struct dfFrameSetup fs;

	int ifs = df->iFrameSetup;

	vin->dTime = dTime;
	vin->dYearUnit = df->dYearUnit;
	vin->dStep = dStep;
	vin->dExp = dExp;
	vin->bVDetails = df->bVDetails;
	vin->dMassStarMin = df->dMassStarMin;
	vin->dMassStarMax = df->dMassStarMax;
	vin->dMassGasMin  = df->dMassGasMin;
	vin->dMassGasMax  = df->dMassGasMax;
	vin->dMassDarkMin = df->dMassDarkMin;
	vin->dMassDarkMax = df->dMassDarkMax;

	if (df->bLoop) {
		dTime -= df->dTimeMod;
		if (dTime > df->dPeriodLoop+df->dTimeLoop) {
			dTime -= df->dPeriodLoop;
			df->dTimeMod += df->dPeriodLoop;
			ifs = -1;
			}
		}
	
	if (ifs == -1) {
		assert( df->nFrameSetup > 1 );
		if (dTime < df->fs[1].dTime) { /* Outside range */
			fprintf(stderr,"DF WARNING Initial Time outside camera direction table: %g < %g\n",dTime,df->fs[1].dTime);
			df->iFrameSetup = ifs = 0;
			}
		else {
			ifs = 1;
			while (ifs < df->nFrameSetup && dTime >= df->fs[ifs+1].dTime ) ifs++;
			if (ifs >= df->nFrameSetup-1) { /* Outside Range */
				if (dTime == df->fs[df->nFrameSetup-1].dTime && ifs > 1) ifs--;
				else {
					fprintf(stderr,"DF WARNING Initial Time outside camera direction table: %g > %g\n",dTime,df->fs[df->nFrameSetup-1].dTime);
					df->iFrameSetup = ifs = 0;
					}
				}
			df->iFrameSetup = ifs;
			dfGetCoeff( df, ifs );
			}
		}
	else if (df->nFrameSetup > 1) {
		while (dTime > df->fs[ifs+1].dTime || (!ifs && dTime >= df->fs[ifs+1].dTime)) { 
			ifs++;
			if (ifs >= df->nFrameSetup-1) {
				fprintf(stderr,"DF WARNING Time outside camera direction table: %g > %g\n",dTime,df->fs[df->nFrameSetup-1].dTime);
				df->iFrameSetup = ifs = 0;
				break;
				}
			df->iFrameSetup = ifs;
			dfGetCoeff( df, ifs );
			}
		}

	
	if (df->bVDetails) CkPrintf("DF Interpolating at t=%g Setups: %i (t=%g) %i (t=%g)\n",dTime,ifs,df->fs[ifs].dTime,ifs+1,df->fs[ifs+1].dTime);

	fs = df->fs[ifs];

	/* Nothing to interpolate? */
	if (ifs) {
	/* 
	   Interpolate Eye position, FOV etc... 
	   from df->fs[ifs].dTime <= dTime <= df->fs[ifs+1].dTime
	   */
	  dfInterp( df, &fs, (dTime-fs.dTime)*df->rdt );
	  }

	switch (fs.iTarget) {
	case DF_TARGET_COM_GAS:
	  if (com[3]>0) {
		fs.target[0] = com[0]/com[3];
		fs.target[1] = com[1]/com[3];
		fs.target[2] = com[2]/com[3];
	  }
	  break;
	case DF_TARGET_COM_DARK:
	  if (com[7]>0) {
		fs.target[0] = com[4]/com[7];
		fs.target[1] = com[5]/com[7];
		fs.target[2] = com[6]/com[7];
	  }
	  break;
	case DF_TARGET_COM_STAR:
	  if (com[11]>0) {
		fs.target[0] = com[8]/com[11];
		fs.target[1] = com[9]/com[11];
		fs.target[2] = com[10]/com[11];
	  }
	  break;
	case DF_TARGET_COM_ALL:
	  if (com[3]+com[7]+com[11]>0) {
		fs.target[0] = (com[0]+com[4]+com[8])/(com[3]+com[7]+com[11]);
		fs.target[1] = (com[1]+com[5]+com[9])/(com[3]+com[7]+com[11]);
		fs.target[2] = (com[2]+com[6]+com[10])/(com[3]+com[7]+com[11]);
	  }
	  break;
	case DF_TARGET_OLDEST_STAR:
	  if (com[3] < FLT_MAX) {
		fs.target[0] = com[0];
		fs.target[1] = com[1];
		fs.target[2] = com[2];
	  }
	  break;
	case DF_TARGET_PHOTOGENIC:
	  if (com[3] < FLT_MAX) {
		fs.target[0] = com[0]/com[3];
		fs.target[1] = com[1]/com[3];
		fs.target[2] = com[2]/com[3];
	  }
	  break;
	}
	
	dfProjection( vin, &fs, nxPix, nyPix ); 
    
    }


void dfClearImage( struct inDumpFrame *in, void *vImage, int *nImage ) {
	DFIMAGE *Image = vImage;
    int i;
	DFIMAGE blank;

	blank.r = 0;
	blank.g = 0;
    blank.b = 0;

	*nImage = in->nxPix * in->nyPix * sizeof(DFIMAGE);

	for (i=in->nxPix*in->nyPix-1;i>=0;i--) Image[i] = blank;
	}

#ifdef DEBUGTRACK
int trackp = 0;
int trackcnt = 0;
#endif

void dfRenderParticlePoint( struct inDumpFrame *in, void *vImage, 
						 double *r, double fMass, double fSoft, double fBall2, int iActive, double fAge ) {

	DFIMAGE *Image = vImage;
	DFIMAGE col; /* Colour */
	double x,y,z,dr[3];
	int j;
	int xp,yp;

	if ((iActive & in->iTypeGas)) {
		 if (fMass < in->dMassGasMin || fMass > in->dMassGasMax) return;
		 col = in->ColGas;
		 }
	else if ((iActive & in->iTypeDark)) {
		if (fMass < in->dMassDarkMin || fMass > in->dMassDarkMax) return;
		col = in->ColDark;
		}
	else if ((iActive & in->iTypeStar)) {
		if (fMass < in->dMassStarMin || fMass > in->dMassStarMax) return;
		col = in->ColStar;
		}

	for (j=0;j<3;j++) {
		dr[j] = r[j]-in->r[j];
		}
	z = dr[0]*in->z[0] + dr[1]*in->z[1] + dr[2]*in->z[2] + in->zEye;
	if (z >= in->zClipNear && z <= in->zClipFar) {
		x = dr[0]*in->x[0] + dr[1]*in->x[1] + dr[2]*in->x[2];
		if (in->iProject == DF_PROJECT_PERSPECTIVE) x/=z;
		if (fabs(x)<in->xlim) {
			y = dr[0]*in->y[0] + dr[1]*in->y[1] + dr[2]*in->y[2];
			if (in->iProject == DF_PROJECT_PERSPECTIVE) y/=z;
			if (fabs(y)<in->ylim) {
				xp = x+in->xlim;
				yp = in->ylim-y; /* standard screen convention */
				Image[ xp + yp*in->nxPix ].r += col.r;
				Image[ xp + yp*in->nxPix ].g += col.g;
				Image[ xp + yp*in->nxPix ].b += col.b;
				}
			}
		}
	}

/*
 * Change for ChaNGa: "fBall2" is not the square of fBall.
 */
void dfRenderParticleTSC( struct inDumpFrame *in, void *vImage, 
						 double *r, double fMass, double fSoft, double fBall2, int iActive, double fAge ) {
    DFIMAGE *Image = vImage;
	DFIMAGE col; /* Colour */
	double h;
	double x,y,z,dr[3],br0;
	int hint;
	int j;
	int xp,yp;

	br0=1;
	if ((iActive & in->iTypeGas)) {
		if (fMass < in->dMassGasMin || fMass > in->dMassGasMax) return;
		if (in->bGasSph) h = fBall2*0.5*in->dGasSoftMul;
		else h = fSoft*in->dGasSoftMul;
		col = in->ColGas;
		}
	else if ((iActive & in->iTypeDark)) {
		if (fMass < in->dMassDarkMin || fMass > in->dMassDarkMax) return;
		h = fSoft*in->dDarkSoftMul;
		col = in->ColDark;
		}
	else if ((iActive & in->iTypeStar)) {
		if (fMass < in->dMassStarMin || fMass > in->dMassStarMax) return;
		h = fSoft*in->dStarSoftMul;
		switch (in->iColStarAge) {
		case DF_STAR_AGE_BRIGHT_COLOUR:
		case DF_STAR_AGE_COLOUR:
		  {
		  float al;
		  al = (log(fabs(fAge+1e6))-13.815)*(1/9.6);
		  col.b = in->ColStar.b*fabs(1-0.7*al);
		  col.g = in->ColStar.g*0.4;
		  col.r = in->ColStar.r*(0.4+0.32*al);
		  /*		  printf("star: %g %g %g %g %g  %g %g %g\n",fAge,al,col.r,col.g,col.b,in->ColStar.r,in->ColStar.g,in->ColStar.b);*/
		  }
		  break;
		default:
		  col = in->ColStar;
		}
		switch (in->iColStarAge) {
		case DF_STAR_AGE_BASIC:
		  break;
		case DF_STAR_AGE_BRIGHT:
		case DF_STAR_AGE_BRIGHT_COLOUR:
		  br0 = 1./(fabs(fAge)/1e6+1.);
		  break;
		case DF_STAR_AGE_COLOUR:
		  break;
		}
		}

	if (in->bColMassWeight) br0*=fMass;
	
	for (j=0;j<3;j++) {
		dr[j] = r[j]-in->r[j];
		}

#ifdef DEBUGTRACK
	trackp++;
	if ((trackp<DEBUGTRACK)) {
		printf("p %f %f %f  %f %f\n",r[0],r[1],r[2],fMass,fSoft);
		}
#endif

	z = dr[0]*in->z[0] + dr[1]*in->z[1] + dr[2]*in->z[2] + in->zEye;
	if (z >= in->zClipNear && z <= in->zClipFar) {
		if (in->iProject == DF_PROJECT_PERSPECTIVE) h = h*in->hmul/z;
		else h = h*in->hmul;
		hint = h;
		x = dr[0]*in->x[0] + dr[1]*in->x[1] + dr[2]*in->x[2];
		if (in->iProject == DF_PROJECT_PERSPECTIVE) x/=z;
		if (fabs(x)<in->xlim+hint) {
			y = dr[0]*in->y[0] + dr[1]*in->y[1] + dr[2]*in->y[2];
			if (in->iProject == DF_PROJECT_PERSPECTIVE) y/=z;
			if (fabs(y)<in->ylim+hint) {
				x = x+in->xlim;
				xp = x;
				y = in->ylim-y;  /* standard screen convention */
				yp = y;
				if (hint < 1) {
					Image[ xp + yp*in->nxPix ].r += br0*col.r;
					Image[ xp + yp*in->nxPix ].g += br0*col.g;
					Image[ xp + yp*in->nxPix ].b += br0*col.b;
#ifdef DEBUGTRACK
					trackcnt++;
					if ((trackcnt<DEBUGTRACK)) {
						printf("%i %i %f %f %f*\n",xp,yp,col.r,col.g,col.b );
						}
#endif
					}
				else {
					int xpmin,xpmax,ypmin,ypmax,ix,iy;
					DFIMAGE *Imagey;
					double br,br1,r2,ih2;
					ih2 = 1./(h*h);
					br1 = br0*(6/(2.0*3.1412))*ih2;
					xpmin = xp - hint; if (xpmin<0) xpmin=0;
					xpmax = xp + hint; if (xpmax>=in->nxPix) xpmax=in->nxPix-1;
					ypmin = yp - hint; if (ypmin<0) ypmin=0;
					ypmax = yp + hint; if (ypmax>=in->nyPix) ypmax=in->nyPix-1;
					for (iy=ypmin,Imagey = Image + iy*in->nxPix;iy<=ypmax;iy++,Imagey += in->nxPix) {
						for (ix=xpmin;ix<=xpmax;ix++) {
							r2 = ((ix-x)*(ix-x)+(iy-y)*(iy-y))*ih2;
							if (r2 > 1) continue;
							br = br1*(1.0-sqrt(r2));
							Imagey[ ix ].r += br*col.r;
							Imagey[ ix ].g += br*col.g;
							Imagey[ ix ].b += br*col.b;
#ifdef DEBUGTRACK
							trackcnt++;
							if ((trackcnt<DEBUGTRACK)) {
								printf("%i %i %f %f %f\n",ix,iy,br*col.r,br*col.g,br*col.b );
								}
#endif
							}
						}
					}
				}
			}
		}
	}

void dfRenderParticleSolid( struct inDumpFrame *in, void *vImage, 
						 double *r, double fMass, double fSoft, double fBall2, int iActive, double fAge ) {

	DFIMAGE *Image = vImage;
	DFIMAGE col; /* Colour */
	double h;
	double x,y,z,dr[3],br0;
	int hint;
	int j;
	int xp,yp;
		
	br0=1;

	if ((iActive & in->iTypeGas )) {
		if (fMass < in->dMassGasMin || fMass > in->dMassGasMax) return;
		if (in->bGasSph) h = fBall2*0.5*in->dGasSoftMul;
		else h = fSoft*in->dGasSoftMul;
		col = in->ColGas;
		}
	else if ((iActive  & in->iTypeDark )) {
		if (fMass < in->dMassDarkMin || fMass > in->dMassDarkMax) return;
		h = fSoft*in->dDarkSoftMul;
		col = in->ColDark;
		}
	else if ((iActive & in->iTypeStar )) {
		if (fMass < in->dMassStarMin || fMass > in->dMassStarMax) return;
		h = fSoft*in->dStarSoftMul;
		switch (in->iColStarAge) {
		case DF_STAR_AGE_BRIGHT_COLOUR:
		case DF_STAR_AGE_COLOUR:
		  {
		  float al;
		  al = (log(fabs(fAge+1e6))-13.815)*(1/9.6);
		  col.b = in->ColStar.b*fabs(1-0.7*al);
		  col.g = in->ColStar.g*0.4;
		  col.r = in->ColStar.r*(0.4+0.32*al);
		  /*		  printf("star: %g %g %g %g %g  %g %g %g\n",fAge,al,col.r,col.g,col.b,in->ColStar.r,in->ColStar.g,in->ColStar.b); */
		  }
		  break;
		default:
		  col = in->ColStar;
		}
		switch (in->iColStarAge) {
		case DF_STAR_AGE_BASIC:
		  break;
		case DF_STAR_AGE_BRIGHT:
		case DF_STAR_AGE_BRIGHT_COLOUR:
		  br0 = 1./(fabs(fAge)/1e6+1.);
		  break;
		}
		}
	
	if (in->bColMassWeight) br0*=fMass;
	
	for (j=0;j<3;j++) {
		dr[j] = r[j]-in->r[j];
		}
	z = dr[0]*in->z[0] + dr[1]*in->z[1] + dr[2]*in->z[2] + in->zEye;
	if (z >= in->zClipNear && z <= in->zClipFar) {
		if (in->iProject == DF_PROJECT_PERSPECTIVE) h = h*in->hmul/z;
		else h = h*in->hmul;
		hint = h;
		x = dr[0]*in->x[0] + dr[1]*in->x[1] + dr[2]*in->x[2];
		if (in->iProject == DF_PROJECT_PERSPECTIVE) x/=z;
		if (fabs(x)<in->xlim+hint) {
			y = dr[0]*in->y[0] + dr[1]*in->y[1] + dr[2]*in->y[2];
			if (in->iProject == DF_PROJECT_PERSPECTIVE) y/=z;
			if (fabs(y)<in->ylim+hint) {
				xp = x+in->xlim;
				yp = in->ylim-y; /* standard screen convention */
				if (hint < 1) {
					double br = 0.523599*br0*h*h; /* integral of TSC to h */
					Image[ xp + yp*in->nxPix ].r += br*col.r;
					Image[ xp + yp*in->nxPix ].g += br*col.g;
					Image[ xp + yp*in->nxPix ].b += br*col.b;
					}
				else {
					int xpmin,xpmax,ypmin,ypmax,ix,iy;
					DFIMAGE *Imagey;
					double br,r2,ih2;
					ih2 = 1./(h*h);
					xpmin = xp - hint; if (xpmin<0) xpmin=0;
					xpmax = xp + hint; if (xpmax>=in->nxPix) xpmax=in->nxPix-1;
					ypmin = yp - hint; if (ypmin<0) ypmin=0;
					ypmax = yp + hint; if (ypmax>=in->nyPix) ypmax=in->nyPix-1;
					for (iy=ypmin,Imagey = Image + iy*in->nxPix;iy<=ypmax;iy++,Imagey += in->nxPix) {
						for (ix=xpmin;ix<=xpmax;ix++) {
							if (ix==xp && iy==yp) {
								br = br0*(1.57080-1.04720/h); /* Integral of tsc disc to r=1 */
								}
							else {
								r2 = ((ix-xp)*(ix-xp)+(iy-yp)*(iy-yp))*ih2;
								if (r2 > 1) continue;
								br = br0*(1.0-sqrt(r2));
								}
							Imagey[ ix ].r += br*col.r;
							Imagey[ ix ].g += br*col.g;
							Imagey[ ix ].b += br*col.b;
							}
						}
					}
				}
			}
		}
	}

void dfRenderParticleInit( struct inDumpFrame *in, int iTypeGas, int iTypeDark, int iTypeStar )
{
	in->iTypeGas = iTypeGas;
	in->iTypeDark = iTypeDark;
	in->iTypeStar = iTypeStar;
	}

void dfRenderParticle( struct inDumpFrame *in, void *vImage, 
					  double r[3], double fMass, double fSoft, double fBall2, int iActive, double fTimeForm ) 
{

	switch (in->iRender) {
	case DF_RENDER_POINT:
		dfRenderParticlePoint( in, vImage, r, fMass, fSoft, fBall2, iActive, (in->dTime-fTimeForm)*in->dYearUnit );
		break;
	case DF_RENDER_TSC:
		dfRenderParticleTSC( in, vImage, r, fMass, fSoft, fBall2, iActive, (in->dTime-fTimeForm)*in->dYearUnit );
		break;
	case DF_RENDER_SOLID:
		dfRenderParticleSolid( in, vImage, r, fMass, fSoft, fBall2, iActive, (in->dTime-fTimeForm)*in->dYearUnit );
		break;
	case DF_RENDER_SHINE: /* Not implemented -- just does point */
		dfRenderParticlePoint( in, vImage, r, fMass, fSoft, fBall2, iActive, (in->dTime-fTimeForm)*in->dYearUnit );
		break;
		}
	}


/* 
   This approach allows type checking
   Notes: 
     Structures and arrays are both supported --> a fixed stride is all that matters
     Floats must be double precision
     position must be three contiguous memory locations
*/
void dfRenderParticlesInit( struct inDumpFrame *in, int iTypeGas, int iTypeDark, int iTypeStar,
			    cosmoType *pr, cosmoType *pfMass, cosmoType *pfSoft, double *pfBall2, unsigned int *piActive, double *pfTimeForm,
						   void *p, int sizeofp )
{
	in->offsetp_r = ((char *) pr)-((char *) p);
	in->offsetp_fMass = ((char *) pfMass)-((char *) p);
	in->offsetp_fSoft = ((char *) pfSoft)-((char *) p);
	in->offsetp_fBall2 = ((char *) pfBall2)-((char *) p);
	in->offsetp_iActive = ((char *) piActive)-((char *) p);
	in->offsetp_fTimeForm = ((char *) pfTimeForm)-((char *) p);
    in->sizeofp = sizeofp;
	in->iTypeGas = iTypeGas;
	in->iTypeDark = iTypeDark;
	in->iTypeStar = iTypeStar;
	}


void dfRenderParticles( struct inDumpFrame *in, void *vImage, void *pStore, int n ) {
	int i;
	char *p = pStore;
	int offsetp_r = in->offsetp_r;
	int offsetp_fMass = in->offsetp_fMass;
	int offsetp_fSoft = in->offsetp_fSoft;
	int offsetp_fBall2 = in->offsetp_fBall2;
	int offsetp_iActive = in->offsetp_iActive;
	int offsetp_fTimeForm = in->offsetp_fTimeForm;
	int sizeofp = in->sizeofp;

	/*	printf("DF: ColStarAge Type %d, LogScaleType %d\n",in->iColStarAge,in->iLogScale);*/

	switch (in->iRender) {
	case DF_RENDER_POINT:
		for (i=0;i<n;i++) {
			dfRenderParticlePoint( in, vImage, (double *) (p+offsetp_r), *((double *) (p+offsetp_fMass)), 0, 0, *((int *) (p+offsetp_iActive)), (in->dTime - *((double *) (p+offsetp_fTimeForm)))*in->dYearUnit );
			p += sizeofp;
			}
		break;
	case DF_RENDER_TSC:
		for (i=0;i<n;i++) {
			dfRenderParticleTSC( in, vImage, (double *) (p+offsetp_r), *((double *) (p+offsetp_fMass)), *((double *) (p+offsetp_fSoft)), *((double *) (p+offsetp_fBall2)), *((int *) (p+offsetp_iActive)), (in->dTime - *((double *) (p+offsetp_fTimeForm)))*in->dYearUnit );
			/*			if (!(i%1000)) printf("%d %d : %g %g %g %g\n",i,*((int *) (p+offsetp_iActive)), in->dTime, *((double *) (p+offsetp_fTimeForm)), in->dYearUnit,  (in->dTime - *((double *) (p+offsetp_fTimeForm)))*in->dYearUnit);*/
			p += sizeofp;
			}
		break;
	case DF_RENDER_SOLID:
		for (i=0;i<n;i++) {	
			dfRenderParticleSolid( in, vImage, (double *) (p+offsetp_r), *((double *) (p+offsetp_fMass)), *((double *) (p+offsetp_fSoft)), *((double *) (p+offsetp_fBall2)), *((int *) (p+offsetp_iActive)), (in->dTime - *((double *) (p+offsetp_fTimeForm)))*in->dYearUnit );
			p += sizeofp;
			}
		break;
	case DF_RENDER_SHINE: /* Not implemented -- just does point */
		for (i=0;i<n;i++) {
			dfRenderParticlePoint( in, vImage, (double *) (p+offsetp_r), *((double *) (p+offsetp_fMass)), 0, 0, *((int *) (p+offsetp_iActive)), (in->dTime - *((double *) (p+offsetp_fTimeForm)))*in->dYearUnit );
			p += sizeofp;
			}
		break;
		}
	}

/* Relies on PKD definition -- not portable */
#ifdef PKD_HINCLUDED
void dfRenderImage( PKD pkd, struct inDumpFrame *in, void *vImage ) {
	PARTICLE *p;
	int i;

	p = pkd->pStore;
	if (in->iRender == DF_RENDER_POINT) {
		for (i=0;i<pkd->nLocal;i++) {
			double r[3];
			r[0] = p[i].r[0];
			r[1] = p[i].r[1];
			r[2] = p[i].r[2];
			dfRenderParticlePoint( in, vImage, r, p[i].fMass, 0, 0, p[i].iActive, 0 );
			}
		}
	else if (in->iRender == DF_RENDER_TSC) {
		
		for (i=0;i<pkd->nLocal;i++) {
			double r[3];
			r[0] = p[i].r[0];
			r[1] = p[i].r[1];
			r[2] = p[i].r[2];
			dfRenderParticleTSC( in, vImage, r, p[i].fMass, p[i].fSoft, p[i].fBall2, p[i].iActive, 0 );
			}
		}
	else if (in->iRender == DF_RENDER_SOLID) {
		for (i=0;i<pkd->nLocal;i++) {
			double r[3];
			r[0] = p[i].r[0];
			r[1] = p[i].r[1];
			r[2] = p[i].r[2];
			dfRenderParticleSolid( in, vImage, r, p[i].fMass, p[i].fSoft, p[i].fBall2, p[i].iActive, 0 );
			}
		}
	else if (in->iRender == DF_RENDER_SHINE) { /* Not implemented -- just does point */
		for (i=0;i<pkd->nLocal;i++) {
			double r[3];
			r[0] = p[i].r[0];
			r[1] = p[i].r[1];
			r[2] = p[i].r[2];
			dfRenderParticlePoint( in, vImage, r, p[i].fMass, 0, 0, p[i].iActive, 0 );
			}
		}
	}

void dfRenderImageOld( PKD pkd, struct inDumpFrame *in, void *vImage ) {
	PARTICLE *p;
	DFIMAGE *Image = vImage;
	DFIMAGE col; /* Colour */
	int i,j;
	double x,y,z,dr[3],br0;
	double xlim = (in->nxPix-1)*.5;
	double ylim = (in->nyPix-1)*.5;
	int xp,yp;

	p = pkd->pStore;
	if (in->iRender == DF_RENDER_POINT) {
		for (i=0;i<pkd->nLocal;i++) {
			if (TYPETest( &p[i], TYPE_GAS )) {
				if (p[i].fMass < in->dMassGasMin || p[i].fMass > in->dMassGasMax) continue;
				col = in->ColGas;
				}
			if (TYPETest( &p[i], TYPE_DARK )) {
				if (p[i].fMass < in->dMassDarkMin || p[i].fMass > in->dMassDarkMax) continue;
				col = in->ColDark;
				}
			if (TYPETest( &p[i], TYPE_STAR )) {
				if (p[i].fMass < in->dMassStarMin || p[i].fMass > in->dMassStarMax) continue;
				col = in->ColStar;
				}
			for (j=0;j<3;j++) {
				dr[j] = p[i].r[j]-in->r[j];
				}
			z = dr[0]*in->z[0] + dr[1]*in->z[1] + dr[2]*in->z[2] + in->zEye;
			if (z >= in->zClipNear && z <= in->zClipFar) {
				x = dr[0]*in->x[0] + dr[1]*in->x[1] + dr[2]*in->x[2];
				if (in->iProject == DF_PROJECT_PERSPECTIVE) x/=z;
				if (fabs(x)<xlim) {
					y = dr[0]*in->y[0] + dr[1]*in->y[1] + dr[2]*in->y[2];
					if (in->iProject == DF_PROJECT_PERSPECTIVE) y/=z;
					if (fabs(y)<ylim) {
						xp = x+xlim;
						yp = ylim-y; /* standard screen convention */
						Image[ xp + yp*in->nxPix ].r += col.r;
						Image[ xp + yp*in->nxPix ].g += col.g;
						Image[ xp + yp*in->nxPix ].b += col.b;
						}
					}
				}
			}
		}
	else if (in->iRender == DF_RENDER_TSC) {
		double hmul = 4*sqrt(in->x[0]*in->x[0] + in->x[1]*in->x[1] + in->x[2]*in->x[2]),h=0.0/*initialized to suppress warning: DCR 12/19/02*/;
		int hint;
		
		for (i=0;i<pkd->nLocal;i++) {
			if (TYPETest( &p[i], TYPE_GAS )) {
				if (p[i].fMass < in->dMassGasMin || p[i].fMass > in->dMassGasMax) continue;
				if (in->bGasSph) h = sqrt(p[i].fBall2)*0.5*in->dGasSoftMul;
				else h = p[i].fSoft*in->dGasSoftMul;
				col = in->ColGas;
				}
			if (TYPETest( &p[i], TYPE_DARK )) {
				if (p[i].fMass < in->dMassDarkMin || p[i].fMass > in->dMassDarkMax) continue;
				h = p[i].fSoft*in->dDarkSoftMul;
				col = in->ColDark;
				}
			if (TYPETest( &p[i], TYPE_STAR )) {
				if (p[i].fMass < in->dMassStarMin || p[i].fMass > in->dMassStarMax) continue;
				h = p[i].fSoft*in->dStarSoftMul;
				col = in->ColStar;
				}

			if (in->bColMassWeight) br0=p[i].fMass;

			for (j=0;j<3;j++) {
				dr[j] = p[i].r[j]-in->r[j];
				}

#ifdef DEBUGTRACK
	trackp++;
	if ((trackp<DEBUGTRACK)) {
		printf("p %f %f %f  %f %f\n",p[i].r[0],p[i].r[1],p[i].r[2],p[i].fMass,p[i].fSoft);
		}
#endif

			z = dr[0]*in->z[0] + dr[1]*in->z[1] + dr[2]*in->z[2] + in->zEye;
			if (z >= in->zClipNear && z <= in->zClipFar) {
				if (in->iProject == DF_PROJECT_PERSPECTIVE) h = h*hmul/z;
				else h = h*hmul;
				hint = h;
				x = dr[0]*in->x[0] + dr[1]*in->x[1] + dr[2]*in->x[2];
				if (in->iProject == DF_PROJECT_PERSPECTIVE) x/=z;
				if (fabs(x)<xlim+hint) {
					y = dr[0]*in->y[0] + dr[1]*in->y[1] + dr[2]*in->y[2];
					if (in->iProject == DF_PROJECT_PERSPECTIVE) y/=z;
					if (fabs(y)<ylim+hint) {
						x = x+xlim;
						xp = x;
						y = ylim-y;  /* standard screen convention */
						yp = y;
						if (hint < 1) {
							Image[ xp + yp*in->nxPix ].r += col.r;
							Image[ xp + yp*in->nxPix ].g += col.g;
							Image[ xp + yp*in->nxPix ].b += col.b;
#ifdef DEBUGTRACK
					trackcnt++;
					if ((trackcnt<DEBUGTRACK)) {
						printf("%i %f %f %f*\n",xp,col.r,col.g,col.b );
						}
#endif
							}
						else {
							int xpmin,xpmax,ypmin,ypmax,ix,iy;
							DFIMAGE *Imagey;
							double br,br1,r2,ih2;
							ih2 = 1./(h*h);
							br1 = (6/(2.0*3.1412))*ih2;
							xpmin = xp - hint; if (xpmin<0) xpmin=0;
							xpmax = xp + hint; if (xpmax>=in->nxPix) xpmax=in->nxPix-1;
							ypmin = yp - hint; if (ypmin<0) ypmin=0;
							ypmax = yp + hint; if (ypmax>=in->nyPix) ypmax=in->nyPix-1;
							for (iy=ypmin,Imagey = Image + iy*in->nxPix;iy<=ypmax;iy++,Imagey += in->nxPix) {
								for (ix=xpmin;ix<=xpmax;ix++) {
									r2 = ((ix-x)*(ix-x)+(iy-y)*(iy-y))*ih2;
									if (r2 > 1) continue;
									br = br1*(1.0-sqrt(r2));
									Imagey[ ix ].r += br*col.r;
									Imagey[ ix ].g += br*col.g;
									Imagey[ ix ].b += br*col.b;
#ifdef DEBUGTRACK
									trackcnt++;
									if ((trackcnt<DEBUGTRACK)) {
										printf("%i %f %f %f\n",ix,br*col.r,br*col.g,br*col.b );
										}
#endif
									}
								}
							}
						}
					}
				}
			}
		}
	else if (in->iRender == DF_RENDER_SOLID) {
		double hmul = 4*sqrt(in->x[0]*in->x[0] + in->x[1]*in->x[1] + in->x[2]*in->x[2]),h=0.0/*initialized to suppress warning: DCR 12/19/02*/;
		int hint;
		
		br0=1;

		for (i=0;i<pkd->nLocal;i++) {
		        if (in->bColMassWeight) br0=p[i].fMass;
			if (TYPETest( &p[i], TYPE_GAS )) {
				if (p[i].fMass < in->dMassGasMin || p[i].fMass > in->dMassGasMax) continue;
				if (in->bGasSph) h = sqrt(p[i].fBall2)*0.5*in->dGasSoftMul;
				else h = p[i].fSoft*in->dGasSoftMul;
				col = in->ColGas;
				}
			if (TYPETest( &p[i], TYPE_DARK )) {
				if (p[i].fMass < in->dMassDarkMin || p[i].fMass > in->dMassDarkMax) continue;
				h = p[i].fSoft*in->dDarkSoftMul;
				col = in->ColDark;
				}
			if (TYPETest( &p[i], TYPE_STAR )) {
				if (p[i].fMass < in->dMassStarMin || p[i].fMass > in->dMassStarMax) continue;
				h = p[i].fSoft*in->dStarSoftMul;
				col = in->ColStar;
				}


			for (j=0;j<3;j++) {
				dr[j] = p[i].r[j]-in->r[j];
				}
			z = dr[0]*in->z[0] + dr[1]*in->z[1] + dr[2]*in->z[2] + in->zEye;
			if (z >= in->zClipNear && z <= in->zClipFar) {
				if (in->iProject == DF_PROJECT_PERSPECTIVE) h = h*hmul/z;
				else h = h*hmul;
				hint = h;
				x = dr[0]*in->x[0] + dr[1]*in->x[1] + dr[2]*in->x[2];
				if (in->iProject == DF_PROJECT_PERSPECTIVE) x/=z;
				if (fabs(x)<xlim+hint) {
					y = dr[0]*in->y[0] + dr[1]*in->y[1] + dr[2]*in->y[2];
					if (in->iProject == DF_PROJECT_PERSPECTIVE) y/=z;
					if (fabs(y)<ylim+hint) {
						xp = x+xlim;
						yp = ylim-y; /* standard screen convention */
						if (hint < 1) {
							double br = 0.523599*br0*h*h; /* integral of TSC to h */
							Image[ xp + yp*in->nxPix ].r += br*col.r;
							Image[ xp + yp*in->nxPix ].g += br*col.g;
							Image[ xp + yp*in->nxPix ].b += br*col.b;
							}
						else {
							int xpmin,xpmax,ypmin,ypmax,ix,iy;
							DFIMAGE *Imagey;
							double br,r2,ih2;
							ih2 = 1./(h*h);
							xpmin = xp - hint; if (xpmin<0) xpmin=0;
							xpmax = xp + hint; if (xpmax>=in->nxPix) xpmax=in->nxPix-1;
							ypmin = yp - hint; if (ypmin<0) ypmin=0;
							ypmax = yp + hint; if (ypmax>=in->nyPix) ypmax=in->nyPix-1;
							for (iy=ypmin,Imagey = Image + iy*in->nxPix;iy<=ypmax;iy++,Imagey += in->nxPix) {
								for (ix=xpmin;ix<=xpmax;ix++) {
									if (ix==xp && iy==yp) {
										br = br0*(1.57080-1.04720/h); /* Integral of tsc disc to r=1 */
										}
									else {
										r2 = ((ix-xp)*(ix-xp)+(iy-yp)*(iy-yp))*ih2;
										if (r2 > 1) continue;
										br = br0*(1.0-sqrt(r2));
										}
									Imagey[ ix ].r += br*col.r;
									Imagey[ ix ].g += br*col.g;
									Imagey[ ix ].b += br*col.b;
									}
								}
							}
						}
					}
				}
			}
		}
	else if (in->iRender == DF_RENDER_SHINE) {
		for (i=0;i<pkd->nLocal;i++) {
			if (TYPETest( &p[i], TYPE_GAS )) {
				if (p[i].fMass < in->dMassGasMin || p[i].fMass > in->dMassGasMax) continue;
				col = in->ColGas;
				}
			if (TYPETest( &p[i], TYPE_DARK )) {
				if (p[i].fMass < in->dMassDarkMin || p[i].fMass > in->dMassDarkMax) continue;
				col = in->ColDark;
				}
			if (TYPETest( &p[i], TYPE_STAR )) {
				if (p[i].fMass < in->dMassStarMin || p[i].fMass > in->dMassStarMax) continue;
				col = in->ColStar;
				}
			for (j=0;j<3;j++) {
				dr[j] = p[i].r[j]-in->r[j];
				}
			z = dr[0]*in->z[0] + dr[1]*in->z[1] + dr[2]*in->z[2] + in->zEye;
			if (z >= in->zClipNear && z <= in->zClipFar) {
				x = dr[0]*in->x[0] + dr[1]*in->x[1] + dr[2]*in->x[2];
				if (in->iProject == DF_PROJECT_PERSPECTIVE) x/=z;
				if (fabs(x)<xlim) {
					y = dr[0]*in->y[0] + dr[1]*in->y[1] + dr[2]*in->y[2];
					if (in->iProject == DF_PROJECT_PERSPECTIVE) y/=z;
					if (fabs(y)<ylim) {
						xp = x+xlim;
						yp = ylim-y; /* standard screen convention */
						if (col.r > Image[ xp + yp*in->nxPix ].r) Image[ xp + yp*in->nxPix ].r = col.r;
						if (col.g > Image[ xp + yp*in->nxPix ].r) Image[ xp + yp*in->nxPix ].g = col.g;
						if (col.b > Image[ xp + yp*in->nxPix ].r) Image[ xp + yp*in->nxPix ].b = col.b;
						}
					}
				}
			}
		}
	}
#endif

void dfMergeImage( struct inDumpFrame *in, void *vImage1, int *nImage1, void *vImage2, int *nImage2 ) {
	int i;
	DFIMAGE *Image1 = vImage1, *Image2 = vImage2;

	switch (in->iRender) {
	case DF_RENDER_POINT:
	case DF_RENDER_TSC:
	case DF_RENDER_SOLID:
		assert( *nImage1 == in->nxPix*in->nyPix*sizeof(DFIMAGE) );
		assert( *nImage1 == *nImage2 );

		for (i=in->nxPix*in->nyPix-1;i>=0;i--) {
			Image1[i].r += Image2[i].r;
			Image1[i].g += Image2[i].g;
			Image1[i].b += Image2[i].b;
			}
		break;
	case DF_RENDER_SHINE:
		assert( *nImage1 == in->nxPix*in->nyPix*sizeof(DFIMAGE) );
		assert( *nImage1 == *nImage2 );

		for (i=in->nxPix*in->nyPix-1;i>=0;i--) {
			if (Image2[i].r > Image1[i].r) Image1[i].r = Image2[i].r;
			if (Image2[i].g > Image1[i].g) Image1[i].g = Image2[i].g;
			if (Image2[i].b > Image1[i].b) Image1[i].b = Image2[i].b;
			}
		break;
	default:
		break;
		}
	}

void dfFinishFrame( struct DumpFrameContext *df, double dTime, double dStep, struct inDumpFrame *in, void  *vImage, int liveViz, unsigned char **outgray) {
	DFIMAGE *Image = vImage;
	char fileout[160];
	FILE *fp;
	int i;
	int iMax;
	unsigned char *g;
	unsigned char *gray;
	/*char number[40]; -- not used: DCR 12/19/02*/

	iMax = in->nxPix*in->nyPix;
	gray = (unsigned char *) malloc(sizeof(unsigned char)*3*iMax);
	assert( gray != NULL );
	*outgray = gray;

	if (in->iRender == DF_RENDER_POINT) {
		for (i=0,g=gray;i<iMax;i++) {
			if ( Image[i].r > 0.0 ) {
				*g=255; g++; 
				*g=255; g++; 
				*g=255; g++; 
				}
			else {
				*g=0; g++; 
				*g=0; g++; 
				*g=0; g++; 
				}
			}
		}
	else if (in->iRender == DF_RENDER_TSC || in->iRender == DF_RENDER_SOLID || in->iRender == DF_RENDER_SHINE) {
		if (in->iLogScale == DF_LOG_COLOURSAFE) {
			double lmin,factor;
			lmin = log(in->pScale1);
			factor = 255.999/(log(in->pScale2)-lmin);
			for (i=0,g=gray;i<iMax;i++) {
				int bing;
				float tot = Image[i].r;
				if (Image[i].g > tot) tot = Image[i].g;
				if (Image[i].b > tot) tot = Image[i].b;
				if (tot <= 0) {
				  *g = 0; g++;
				  *g = 0; g++;
				  *g = 0; g++;
				}
				else {
				  tot = factor/tot*(log(tot)-lmin);
				  bing = Image[i].r*tot;
				  *g = (bing < 255 ? (bing < 0 ? 0 : bing) : 255 );
				  g++;
				  bing = Image[i].g*tot;
				  *g = (bing < 255 ? (bing < 0 ? 0 : bing) : 255 );
				  g++;
				  bing = Image[i].b*tot;
				  *g = (bing < 255 ? (bing < 0 ? 0 : bing) : 255 );
				  g++;
				  }
			    }
		    }
		else if (in->iLogScale == DF_LOG_SATURATE) {
			double lmin,factor;
			lmin = log(in->pScale1);
			factor = 255.999/(log(in->pScale2)-lmin);
			for (i=0,g=gray;i<iMax;i++) {
				int bing;
				if (Image[i].r <= 0) *g = 0;
				else { 
					bing = factor*(log(Image[i].r)-lmin);
					*g = (bing < 255 ? (bing < 0 ? 0 : bing) : 255 );
					}
				g++;
				if (Image[i].g <= 0) *g = 0;
				else { 
					bing = factor*(log(Image[i].g)-lmin);
					*g = (bing < 255 ? (bing < 0 ? 0 : bing) : 255 );
					}
				g++;
				if (Image[i].b <= 0) *g = 0;
				else { 
					bing = factor*(log(Image[i].b)-lmin);
					*g = (bing < 255 ? (bing < 0 ? 0 : bing) : 255 );
					}
				g++;
				}
			}
		else {
			for (i=0,g=gray;i<iMax;i++) {
				int bing;
				bing = 260*(1.-1./(0.1*Image[i].r+1));
				*g = (bing < 255 ? bing : 255 );
				g++;
				bing = 260*(1.-1./(0.1*Image[i].g+1));
				*g = (bing < 255 ? bing : 255 );
				g++;
				bing = 260*(1.-1./(0.1*Image[i].b+1));
				*g = (bing < 255 ? bing : 255 );
				g++;
				}
			}
	}

	if(!liveViz) {
	  switch( df->iNumbering ) {
	  case DF_NUMBERING_FRAME:
	    sprintf(fileout, df->FileName, df->nFrame );
	    break;
	  case DF_NUMBERING_STEP:
	    sprintf(fileout, df->FileName, dStep );
	    break;
	  case DF_NUMBERING_TIME:
	    sprintf(fileout, df->FileName, dTime );
	    break;
	  }

	df->nFrame++; /* NB: need to sort out something for restarts */
        fp = CmiFopen(fileout,"w");
        if(fp == NULL) {
            fprintf(stderr, "Bad Frame file open: %s\n", fileout);
            if(errno == ENOENT)
                fprintf(stderr, "directory of %s does not exist\n", fileout);
            }
        assert(fp!=NULL);

        if (df->iEncode == DF_ENCODE_PPM) {
		fprintf(fp,"P6\n#T=%20.10f\n%5i %5i\n255\n",dTime,in->nxPix,in->nyPix);
		CmiFwrite( gray, 3*iMax, sizeof(char), fp);
		}

	else if (df->iEncode == DF_ENCODE_PNG) {
#ifdef USE_PNG
		static mainprog_info wpng_info;
		int rowbytes;
		int iErr;
		
		wpng_info.outfile = fp;
	    wpng_info.infile = NULL; 
		wpng_info.image_data = gray; 
		wpng_info.pnmtype = 6; /* RGB */

		wpng_info.sample_depth = 8;  /* <==> maxval 255 */

		wpng_info.width = in->nxPix;
		wpng_info.height = in->nyPix;
		rowbytes = 3*wpng_info.width; /* 8 bit RGB */
		wpng_info.row_pointers = (uch **)malloc(wpng_info.height*sizeof(uch *)); 

		wpng_info.filter = FALSE; 
		wpng_info.interlaced = FALSE; 
		wpng_info.have_bg = FALSE; 
		wpng_info.have_time = FALSE; 
		wpng_info.have_text = 0; 
		wpng_info.gamma = 0.0; 

		iErr =  writepng_init(&wpng_info);
		assert (!iErr);

		for (i = 0;  i < wpng_info.height;  ++i) 
            wpng_info.row_pointers[i] = wpng_info.image_data + i*rowbytes; 
	
        iErr = writepng_encode_image(&wpng_info);
		assert (!iErr);

		writepng_cleanup(&wpng_info);
#endif
	}

	CmiFclose(fp);

	if (df->dDumpFrameTime > 0 && dTime >= df->dTime)
	  df->dTime = df->dTime + df->dDumpFrameTime;
	
	if (df->dDumpFrameStep > 0 && dStep >= df->dStep) 
	  df->dStep = df->dStep + df->dDumpFrameStep;
	
	free(gray);
	}
}
