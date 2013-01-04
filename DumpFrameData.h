#ifndef DUMP_FRAME_DATA_H
#define DUMP_FRAME_DATA_H

#include "ParallelGravity.h"
#include "Reductions.h"

/// @brief Group to hold the frame buffer memory.
class DumpFrameData : public CBase_DumpFrameData 
{
    void *bufImage;
    int nImage;
  public:
    void *Image;

    DumpFrameData() { bufImage = NULL; }
    DumpFrameData(CkMigrateMessage *m) : CBase_DumpFrameData(m) { bufImage = NULL;}
    void pup(PUP::er &p) {}
    ~DumpFrameData() { if(bufImage) free(bufImage); }
    /// @brief
    void clearFrame(InDumpFrame in, const CkCallback& cb) {
	if(bufImage == NULL)
	    bufImage = malloc(sizeof(in) + in.nxPix*in.nyPix*sizeof(DFIMAGE));
	Image = ((char *)bufImage) + sizeof(in);
	*((struct inDumpFrame *)bufImage) = in; //start of reduction
						//message is the parameters
	dfClearImage( &in, Image, &nImage);
	contribute(cb);
	}
    
    void combineFrame(InDumpFrame in, const CkCallback& cb) {
	contribute(sizeof(in) + nImage, bufImage, dfImageReduction, cb);
	free(bufImage);
	bufImage = NULL;
	}
    };

#endif
