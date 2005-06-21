/** @file SPH_Requests.cpp
 Register the PUP routines for the SPH request/response pairs.
 @author Graeme Lufkin (gwl@u.washington.edu)
 */
 
#include <numeric>

#include "SPH_Requests.h"

PUPable_def(BallCountResponse);
PUPable_def(BallCountRequest);
PUPable_def(NeighborSearchResponse);
PUPable_def(NeighborSearchRequest);
PUPable_def(DensityResponse);
PUPable_def(DensityRequest);
PUPable_def(VelDivCurlResponse);
PUPable_def(VelDivCurlRequest);
PUPable_def(VelDispResponse);
PUPable_def(VelDispRequest);

void initSPHPUP() {
	PUPable_reg(BallCountResponse);
	PUPable_reg(BallCountRequest);
	PUPable_reg(NeighborSearchResponse);
	PUPable_reg(NeighborSearchRequest);
	PUPable_reg(DensityResponse);
	PUPable_reg(DensityRequest);
	PUPable_reg(VelDivCurlResponse);
	PUPable_reg(VelDivCurlRequest);
	PUPable_reg(VelDispResponse);
	PUPable_reg(VelDispRequest);
}
