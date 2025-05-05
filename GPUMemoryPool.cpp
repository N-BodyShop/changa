#include <stdio.h>        // Standard C I/O
#include <cstdio>         // C++ wrapper for C stdio
#include <string>         // std::string
#include <iostream>       // std::cerr

#include <cuda_runtime.h> // Core CUDA runtime

#include "charm++.h"         // Core Charm++ definitions MUST COME EARLY
#include "ParallelGravity.h" // Main application header (includes decl.h, param.h etc.) MUST COME EARLY

#include "GPUMemoryPool.h" // Header for this file's functions/macros
#include "DataManager.h"     // Needs Charm++, ParallelGravity types
#include "memlog.h"          // Needs basic types, maybe string
#include "ParallelGravity.decl.h" // Generated header, include AFTER dependencies

cudaError_t gpuMallocTracked(void** devPtr, size_t size, const char* tag, const char* functionTag, const char* file, int line) {
    // Comment out old logging
    // fprintf(stdout, "[GPU Memory][Alloc] %s:%d - %s (%zu bytes)\n", file, line, tag, size);
    // fflush(stdout);

    double timestamp = CkWallTimer(); // Get timestamp before allocation
    cudaError_t result = cudaMalloc(devPtr, size);
    double timestamp_after = CkWallTimer(); // Get timestamp after (more accurate for event time)

    // Access node-local DataManager
    DataManager* dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    if (dm && dm->memLog) { // Check if dm and memLog are valid
        MemLogOpType opType = (result == cudaSuccess) ? MEMLOG_ALLOC : MEMLOG_ALLOC_FAIL;
        uintptr_t address = (result == cudaSuccess && *devPtr != NULL) ? (uintptr_t)(*devPtr) : 0;
        MemLogEvent event(CkMyNode(), opType, size, address, timestamp_after, file, line, tag, functionTag);

        CmiLock(dm->lockMemLog);
        dm->memLog->meTab.push_back(event);
        CmiUnlock(dm->lockMemLog);
    }

    return result;
}

cudaError_t gpuFreeTracked(void* devPtr, const char* tag, const char* functionTag, const char* file, int line) {
    double timestamp = CkWallTimer();
    MemLogOpType opType;
    uintptr_t address = (uintptr_t)devPtr;

    // Access node-local DataManager
    DataManager* dm = (DataManager*)CkLocalNodeBranch(dataManagerID);

    if (devPtr == NULL) {
        // Comment out old logging
        // fprintf(stdout, "[GPU Memory][Free SKIP] %s:%d - %s (NULL pointer)\n", file, line, tag);
        // fflush(stdout);
        opType = MEMLOG_FREE_SKIP;
        if (dm && dm->memLog) {
             MemLogEvent event(CkMyNode(), opType, 0, address, timestamp, file, line, tag, functionTag);
             CmiLock(dm->lockMemLog);
             dm->memLog->meTab.push_back(event);
             CmiUnlock(dm->lockMemLog);
        }
        return cudaSuccess; // Freeing NULL is a no-op in CUDA
    }

    // Comment out old logging
    // fprintf(stdout, "[GPU Memory][Free ] %s:%d - %s (%p)\n", file, line, tag, devPtr);
    // fflush(stdout);

    cudaError_t result = cudaFree(devPtr);
    double timestamp_after = CkWallTimer();

    opType = (result == cudaSuccess) ? MEMLOG_FREE : MEMLOG_FREE_FAIL;

    if (dm && dm->memLog) {
        MemLogEvent event(CkMyNode(), opType, 0, address, timestamp_after, file, line, tag, functionTag);
        CmiLock(dm->lockMemLog);
        dm->memLog->meTab.push_back(event);
        CmiUnlock(dm->lockMemLog);
    }

    return result;
}

// Implementation for the DataManager entry method to set the log filename
// This follows the starlog pattern where the implementation is in the feature's file.
void DataManager::initMemLog(std::string _fileName, const CkCallback &cb) {
    CmiLock(lockMemLog);
    if (memLog != nullptr) { 
        memLog->fileName = _fileName;
    } else {
        CkPrintf("WARNING PE %d: memLog is NULL in initMemLog! Cannot set filename.\n", CkMyPe());
    }
    CmiUnlock(lockMemLog);
    // Signal completion for this PE for the barrier in Main::initMemLog
    contribute(cb);
}

/// @brief Initializes the memory log file on PE 0 and distributes the filename.
/// Mimics the Main::initStarLog pattern.
void Main::initMemLog() {
    std::string memLogFile = "memlog.out"; // Hardcoded filename

    // Create a callback to wait for all DMs to receive the filename
    // Call the init function on the DataManager on all PEs
    dMProxy.initMemLog(memLogFile, CkCallbackResumeThread());
    // Implicit wait for all PEs to finish DataManager::initMemLog happens here

    // PE 0 creates/truncates the file and writes a header AFTER the barrier completes
    if (CkMyPe() == 0) {
        FILE* fpLog = CmiFopen(memLogFile.c_str(), "w");
        if (fpLog == NULL) {
            // Use CkPrintf or CkError for output in Charm++ applications
            CkPrintf("ERROR: PE 0 Cannot create memlog file: %s\n", memLogFile.c_str());
            // Decide on error handling: CkAbort? or just print warning?
            // For now, print and continue, logging might just fail silently later.
        } else {
            fprintf(fpLog, "# ChaNGa Memory Log v1.1\n");
            fprintf(fpLog, "# NodeID OpType Size Address Timestamp File:Line PointerID FunctionTag\n");
            int close_err = CmiFclose(fpLog);
             if (close_err != 0) {
                 CkPrintf("ERROR: PE 0 failed to close memlog file: %s (Error %d)\n", memLogFile.c_str(), close_err);
             }
        }
    }
}

/// @brief Flush memlog table to disk sequentially across nodes.
/// This coordinates the flush; the actual writing happens in MemLog::flush().
void DataManager::flushMemLog(const CkCallback& cb) {
    // Call the actual file writing implementation in MemLog
    // Assumes memLog->flush() handles checking if the buffer is empty,
    // file opening/writing/closing, error checking, and buffer clearing.
    if (memLog) { // Ensure memLog is not null
         memLog->flush();
    } else {
         // Log a warning if memLog is unexpectedly null
         // Use CkPrintf for Charm++ compatible output
         CkPrintf("WARNING Node %d: memLog is NULL in flushMemLog! Skipping flush.\n", thisIndex);
    }

    // Sequential node logic: pass the call along until the last node
    if(thisIndex != CkNumNodes()-1) {
        // Pass the call to the next node, forwarding the final callback
        thisProxy[thisIndex + 1].flushMemLog(cb);
    } else {
        // We are the last node, signal completion of the entire sequence
        cb.send();
    }
}

/// @brief Flush buffered memory log events to the designated file.
/// This function performs the actual file I/O for the memory log.
void MemLog::flush() {
    if (meTab.empty()) {
        return; // Nothing to flush
    }

    FILE* outfile = CmiFopen(fileName.c_str(), "a"); // Open in append mode

    if (outfile == NULL) {
        // Use CkPrintf for Charm++ compatible error output. Avoid aborting for logging failures.
        CkPrintf("ERROR: Could not open memlog file '%s' for appending.\n", fileName.c_str());
        return; // Cannot proceed without the file
    }

    // Iterate through the buffered events and write them to the file
    for (const auto& event : meTab) {
        const char* opTypeStr;
        switch (event.opType) {
            case MEMLOG_ALLOC:      opTypeStr = "ALLOC";      break;
            case MEMLOG_FREE:       opTypeStr = "FREE ";      break; // Padded for alignment
            case MEMLOG_ALLOC_FAIL: opTypeStr = "ALLOC_F";    break;
            case MEMLOG_FREE_FAIL:  opTypeStr = "FREE_F ";    break;
            case MEMLOG_FREE_SKIP:  opTypeStr = "FREE_S ";    break;
            default:                opTypeStr = "UNKNOWN";    break;
        }

        // Format: NodeID OpType Size Address Timestamp File:Line Tag
        // Use %d for NodeID, %zu for size_t, %p for pointer (address), %.6f for timestamp
        // Ensure tag does not contain problematic characters if needed (currently assuming ok)
        fprintf(outfile, "%d %s %zu %p %.6f %s %s %s\n",
                event.nodeId,
                opTypeStr,
                event.size,
                (void*)event.address, // Cast uintptr_t back to void* for %p
                event.timestamp,
                event.location.c_str(),
                event.pointerId.c_str(),
                event.functionTag.c_str());
    }

    int result = CmiFclose(outfile);
    if (result != 0) {
        CkPrintf("WARNING: Failed to close memlog file '%s' properly (Error %d).\n", fileName.c_str(), result);
        // Continue even if close fails, data might still be flushed
    }

    // Clear the buffer now that events are written
    meTab.clear();
} 