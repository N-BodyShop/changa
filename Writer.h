/** @file Writer.h
    The Writer group holds a copy of the particles before they are
    written to disk.  It sorts them and writes the requested files.
*/
#ifndef _WRITER_H_
#define _WRITER_H_

#include <vector>
#include "ParallelGravity.decl.h"

class Writer : public CBase_Writer {
public:
    
private:
    std::vector<GravityParticle> vMyParticles;
    std::vector<extraSPHData> vMySPHParticles;
    std::vector<extraStarData> vMyStarParticles;
    int64_t nTotalParticles;
    int64_t nTotalSPH;
    int64_t nTotalDark;
    int64_t nTotalStar;
    /// A local pointer to my DataManager.
    DataManager* dm;
    /// Setup for writing
    int nSetupWriteStage;
    int64_t nStartWrite;	// Particle number at which this piece starts
				// to write file.
    /// used to determine which coordinate we are outputting, while printing
    /// accelerations in ASCII format
    int cnt;
    /// used to determine if the x, y, z coordinates should be printed
    /// together while outputting accelerations in ASCII format
    int packed;
public:
    Writer() {
        nStartWrite = -1;
        nSetupWriteStage = -1;
        packed=0;			/* output accelerations in x,y,z
					   packed format */
        cnt=0;			/* Counter over x,y,z when not packed */
    }
    Writer(CkMigrateMessage *m) : CBase_Writer(m) {
        nStartWrite = -1;
        nSetupWriteStage = -1;
        packed=0;
        cnt=0;
    }
    void pup(PUP::er &p){
        CBase_Writer::pup(p);
    }
    
    void init(const CkCallback& cb);
    void receive(std::vector<GravityParticle> vPart, std::vector<extraSPHData>,
                 std::vector<extraStarData>);
    void clear(int64_t nTotal, int64_t nSPH, int64_t nDark, int64_t nStar,
               const CkCallback& cb);
    void reOrder(const CkCallback &cb);
	// control parallelism in the write
	void parallelWrite(int iPass, const CkCallback& cb,
			   const std::string& filename, const double dTime,
			   const double dvFac, // scale velocities
			   const double duTFac, // convert temperature
                           const bool bDoublePos,
                           const bool bDoubleVel,
			   const int bCool);
	// serial output
	void serialWrite(u_int64_t iPrevOffset, const std::string& filename,
			 const double dTime,
			 const double dvFac, const double duTFac,
                         const bool bDoublePos,
                         const bool bDoubleVel,
			 const int bCool, const CkCallback& cb);
	// setup for serial output
	void oneNodeWrite(int iIndex,
                          int iOutSPH,
                          int iOutStar,
                          std::vector<GravityParticle> vParticles, // particles to write
                          std::vector<extraSPHData> vGas, // SPH data
                          std::vector<extraStarData> vStar, // Star data
                          int *piSPH, // SPH data offsets
                          int *piStar, // Star data offsets
                          const u_int64_t iPrevOffset,
                          const std::string& filename,  // output file
                          const double dTime,      // time or expansion
                          const double dvFac,  // velocity conversion
                          const double duTFac, // temperature
						  // conversion
                          const bool bDoublePos,
                          const bool bDoubleVel,
			  const int bCool,
			  const CkCallback &cb);
    void setupWrite(int iStage, // stage of scan
                    u_int64_t iPrevOffset,
                    const std::string& filename,
                    const double dTime,
                    const double dvFac,
                    const double duTFac,
                    const bool bDoublePos,
                    const bool bDoubleVel,
                    const int bCool,
                    const CkCallback& cb);
    void outputASCII(OutputParams& params, int bParaWrite,
                     const CkCallback& cb);
    void oneNodeOutVec(OutputParams& params, Vector3D<double>* avOut,
                       int nPart, int iIndex, int bDone,
                       const CkCallback& cb) ;
    void oneNodeOutArr(OutputParams& params, double* adOut,
                       int nPart, int iIndex, int bDone,
                       const CkCallback& cb) ;
    void oneNodeOutIntArr(OutputParams& params, int *aiOut,
                          int nPart, int iIndex, const CkCallback& cb);
    void outputBinary(Ck::IO::Session session, OutputParams& params);
    void outputBinaryStart(OutputParams& params, int64_t nStart,
                           const CkCallback& cb);
    void writeTipsy(Tipsy::TipsyWriter& w,
                    std::vector<GravityParticle>& vParticles,
                    const double dvFac, // scale velocities
                    const double duTFac, // convert temperature
                    const bool bDoublePos,
                    const bool bDoubleVel,
                    const int bCool);
};


#endif /* _WRITER_H_ */
