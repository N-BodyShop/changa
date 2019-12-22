#include "PEUnshuffle.h"

/// @brief create unshuffle messages based on all TreePieces on this
/// PE.
/// @param treePiece Pointer to calling TreePiece.
///
/// This is a replacement for TreePiece::sendParticlesDuringDD(),
/// except that all the TreePieces on a PE are processed in one go.
/// The aim is to reduce the number of unShuffle messages by a factor
/// of the mean number of TreePieces per PE.

void PEUnshuffle::sendParticlesDuringDD(TreePiece *treePiece) {

    vtpLocal.push_back(treePiece);

    // On first call, find the total number of non-empty pieces on this PE.
    // The charm++ location manager gives us this count in cTreePieces
    if(vtpLocal.length() == 1) {
        CkLocMgr *locMgr = treeProxy.ckLocMgr();  
        locMgr->iterate(cTreePieces);
    }

    // check if we have everyone
    if(vtpLocal.length() < cTreePieces.count)
        return;
            
    /// Needed to get boundary keys
    DataManager *dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    /// Store total load on this PE for load balancing calculations
    double tpLoad = 0.0;
    
    /// The last domain
    vector<SFC::Key>::const_iterator endKeys = dm->boundaryKeys.end();
    /// for each TreePiece, the beginning of the domain
    vector<GravityParticle *>vBinBegin;
    /// for each TreePiece, the end of the domain
    vector<GravityParticle *>vBinEnd;

    for(int i = 0; i < vtpLocal.length(); i++){
        TreePiece *pTreePiece = vtpLocal[i];
        tpLoad += pTreePiece->getObjTime();
        vBinBegin.push_back(&(pTreePiece->getParticles()[1]));
    }

    /// Iterator over all domains
    /// N.B. There could be an efficiency gained by looking for the
    /// minimum over all present TreePieces.
    vector<SFC::Key>::iterator iter = dm->boundaryKeys.begin();
    /// There is an offset between the boundaryKey index and the index
    /// into the array of TreePieces.
    int offset = iter - dm->boundaryKeys.begin() - 1;
    /// Iterator for the TreePieces associated with the domains.
    vector<int>::iterator responsibleIter = dm->responsibleIndex.begin() + offset;

    // Loop over domains
    for( ; iter != endKeys; ++iter, ++responsibleIter) {
        /// accumulators particles.
        int nPartOut = 0;
        int nGasOut = 0;
        int nStarOut = 0;
        int saved_phase_len = 0;
        /// Needed for searches over particles.
        GravityParticle dummy;
        dummy.key = *iter;
        // First loop over TreePieces is to find totals
        for(int iTp = 0; iTp < vtpLocal.length(); iTp++){
            TreePiece *pTreePiece = vtpLocal[iTp];
            GravityParticle *binEnd;
            // find particles in this domain
            binEnd = upper_bound(vBinBegin[iTp],
                                 &(pTreePiece->getParticles()[pTreePiece->getNumParticles()+1]),
                                 dummy);
            vBinEnd.push_back(binEnd);
            nPartOut += binEnd - vBinBegin[iTp];
            if(saved_phase_len < pTreePiece->savedPhaseLoad.size())
                saved_phase_len = pTreePiece->savedPhaseLoad.size();
            
            for(GravityParticle *pPart = vBinBegin[iTp]; pPart < binEnd; pPart++) {
                if(pPart->isGas())
                    nGasOut++;
                if(pPart->isStar())
                    nStarOut++;
            }
        }
        
        if(nPartOut > 0) {
            // We have particles to send for this domain.  Allocate
            // message, populate it and send it out.
            ParticleShuffleMsg *shuffleMsg
                = new (saved_phase_len, saved_phase_len, nPartOut, nGasOut, nStarOut)
                ParticleShuffleMsg(saved_phase_len, nPartOut, nGasOut, nStarOut, 0.0);
            memset(shuffleMsg->parts_per_phase, 0, saved_phase_len*sizeof(unsigned int));

            /// Total particles on this PE are needed for load calculations.
            int64_t nPEParticles = 0;
            for(int iTp = 0; iTp < vtpLocal.length(); iTp++){
                TreePiece *pTreePiece = vtpLocal[iTp];
                nPEParticles += pTreePiece->getNumParticles();
                // Calculate the number of particles leaving the
                // treepiece per phase
                for(GravityParticle *pPart = vBinBegin[iTp];
                    pPart < vBinEnd[iTp]; pPart++) {
                    for(int i = 0; i < saved_phase_len; i++) {
                        if (pPart->rung >= i) {
                            shuffleMsg->parts_per_phase[i] += 1;
                        }
                    }
                    if(pTreePiece->havePhaseData(PHASE_FEEDBACK)
                       && (pPart->isGas() || pPart->isStar()))
                        shuffleMsg->parts_per_phase[PHASE_FEEDBACK] += 1;
                }

            }
            // The total load represented by the particles in this message.
            shuffleMsg->load = tpLoad * nPartOut / nPEParticles;
            memset(shuffleMsg->loads, 0.0, saved_phase_len*sizeof(double));

            // Calculate the partial load per phase
            for (int i = 0; i < saved_phase_len; i++) {
                /// The number of particles in this phase (rung)
                int64_t iPESavedPhasePart = 0;
                double dPELoad = 0.0;
                /// Load from the last time Rung0 was calculated
                double dPELoad0 = 0.0;
                for(int iTp = 0; iTp < vtpLocal.length(); iTp++){
                    TreePiece *pTreePiece = vtpLocal[iTp];
                    if (pTreePiece->havePhaseData(i)) {
                        iPESavedPhasePart += pTreePiece->savedPhaseParticle[i];
                        dPELoad += pTreePiece->savedPhaseLoad[i];
                        dPELoad0 += pTreePiece->savedPhaseLoad[0];
                    }
                }
                if(iPESavedPhasePart > 0) {
                    double dLoadFrac = shuffleMsg->parts_per_phase[i]
                        / (double) iPESavedPhasePart;
                    /*
                     * The following can happen if the number of particles on
                     * a given rung increases significantly because of a
                     * timestep adjustment.
                     */
                    if (dLoadFrac > 1.0) dLoadFrac = 1.0;
                    shuffleMsg->loads[i] = dPELoad * dLoadFrac;
                } else if(dPELoad0 > 0.0) {
                    shuffleMsg->loads[i] = dPELoad0 *
                        (shuffleMsg->parts_per_phase[i] / (double) nPEParticles);
                }
            }

            if (verbosity>=3)
                CkPrintf("me:%d to:%d nPart :%d, nGas:%d, nStar: %d\n",
                         CkMyPe(), *responsibleIter,nPartOut, nGasOut,
                         nStarOut);

            int iGasOut = 0;
            int iStarOut = 0;
            GravityParticle *pPartOut = shuffleMsg->particles;
            // Now fill in the particle data for this message
            for(int iTp = 0; iTp < vtpLocal.length(); iTp++){
                TreePiece *pTreePiece = vtpLocal[iTp];
                for(GravityParticle *pPart = vBinBegin[iTp]; pPart < vBinEnd[iTp];
                    pPart++, pPartOut++) {
                    *pPartOut = *pPart;
                    if(pPart->isGas()) {
                        shuffleMsg->pGas[iGasOut]
                            = *(extraSPHData *)pPart->extraData;
                        iGasOut++;
                    }
                    if(pPart->isStar()) {
                        shuffleMsg->pStar[iStarOut]
                            = *(extraStarData *)pPart->extraData;
                        iStarOut++;
                    }
                }
            }
            treeProxy[*responsibleIter].acceptSortedParticles(shuffleMsg);
        }
        // Update domain boundaries for next iteration
        for(int iTp = 0; iTp < vtpLocal.length(); iTp++) {
            vBinBegin[iTp] = vBinEnd[iTp];
        }
        vBinEnd.clear();
    }
    // Mark that the TreePieces are ready to process incoming
    // particles and make sure the processing happens.
    for(int iTp = 0; iTp < vtpLocal.length(); iTp++) {
        vtpLocal[iTp]->incomingParticlesSelf = true;
        vtpLocal[iTp]->acceptSortedParticles(NULL);
    }
    // reset for next call.
    cTreePieces.reset();
    vtpLocal.length() = 0;
}
