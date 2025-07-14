#include "charm++.h"         // Core Charm++ definitions MUST COME EARLY
#include "ParallelGravity.h" // Main application header (includes decl.h, param.h etc.) MUST COME EARLY

#include "GPUMemoryPool.h" // Header for this file's functions/macros
#include "ParallelGravity.decl.h" // Generated header, include AFTER dependencies

#ifdef CUDA

/**
 * @brief Wrapper for cudaMalloc that logs allocation events for memory tracking.
 * @param devPtr Pointer to allocated device memory (output)
 * @param size Size of memory to allocate in bytes
 * @param tag Identifier for the memory allocation (typically variable name)
 * @param functionTag Name of the calling function for context
 * @param file Source file where allocation occurs
 * @param line Line number where allocation occurs
 * @return cudaError_t result from cudaMalloc
 */
cudaError_t gpuMallocTracked(void** devPtr, size_t size, const char* tag, const char* functionTag, const char* file, int line) {

    double timestamp = CkWallTimer(); // Get timestamp before allocation
    cudaError_t result = cudaMalloc(devPtr, size);
    double timestamp_after = CkWallTimer(); // Get timestamp after allocation completes

    // Access node-local DataManager
    DataManager* dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    if (dm && dm->memLog && dm->bGpuMemLogger) { // Check if dm, memLog are valid and logging is enabled
        // Determine operation type based on allocation success/failure
        MemLogOpType opType = (result == cudaSuccess) ? MEMLOG_ALLOC : MEMLOG_ALLOC_FAIL;
        // Get allocated address if successful, otherwise record 0
        uintptr_t address = (result == cudaSuccess && *devPtr != NULL) ? (uintptr_t)(*devPtr) : 0;
        MemLogEvent event(CkMyNode(), opType, size, address, timestamp_after, file, line, tag, functionTag);

        CmiLock(dm->lockMemLog);
        dm->memLog->meTab.push_back(event);
        CmiUnlock(dm->lockMemLog);
    }

    return result;
}

/**
 * @brief Wrapper for cudaFree that logs deallocation events for memory tracking.
 * @param devPtr Pointer to device memory to free
 * @param tag Identifier for the memory deallocation (typically variable name)
 * @param functionTag Name of the calling function for context
 * @param file Source file where deallocation occurs
 * @param line Line number where deallocation occurs
 * @return cudaError_t result from cudaFree (cudaSuccess for NULL pointer)
 */
cudaError_t gpuFreeTracked(void* devPtr, const char* tag, const char* functionTag, const char* file, int line) {
    double timestamp = CkWallTimer();
    MemLogOpType opType;
    uintptr_t address = (uintptr_t)devPtr;

    // Access node-local DataManager
    DataManager* dm = (DataManager*)CkLocalNodeBranch(dataManagerID);

    if (devPtr == NULL) {
        // NULL pointer free is a no-op in CUDA, log as skipped operation
        opType = MEMLOG_FREE_SKIP;
        if (dm && dm->memLog && dm->bGpuMemLogger) {
             MemLogEvent event(CkMyNode(), opType, 0, address, timestamp, file, line, tag, functionTag);
             CmiLock(dm->lockMemLog);
             dm->memLog->meTab.push_back(event);
             CmiUnlock(dm->lockMemLog);
        }
        return cudaSuccess; 
    }

    cudaError_t result = cudaFree(devPtr);
    double timestamp_after = CkWallTimer();

    // Determine operation type based on free success/failure
    opType = (result == cudaSuccess) ? MEMLOG_FREE : MEMLOG_FREE_FAIL;

    if (dm && dm->memLog && dm->bGpuMemLogger) {
        MemLogEvent event(CkMyNode(), opType, 0, address, timestamp_after, file, line, tag, functionTag);
        CmiLock(dm->lockMemLog);
        dm->memLog->meTab.push_back(event);
        CmiUnlock(dm->lockMemLog);
    }

    return result;
}

// Implementation for the DataManager entry method to set the log filename
// This follows the starlog pattern where the implementation is in the feature's file.
void DataManager::initMemLog(std::string _fileName, int bGpuMemLoggerFlag, const CkCallback &cb) {
    CmiLock(lockMemLog);
    if (memLog != nullptr) { 
        memLog->fileName = _fileName;
    } else {
        CkPrintf("WARNING PE %d: memLog is NULL in initMemLog! Cannot set filename.\n", CkMyPe());
    }
    bGpuMemLogger = bGpuMemLoggerFlag; // Set the logging flag
    CmiUnlock(lockMemLog);
    // Signal completion for this PE in the collective operation
    contribute(cb);
}

/// @brief Initializes the memory log file on PE 0 and sends filename to all DataManagers.
/// Mimics the Main::initStarLog pattern.
void Main::initMemLog() {
    std::string memLogFile = "memlog.out"; // Hardcoded filename

    // Send filename and flag to all DataManagers and wait for completion
    // Call the init function on the DataManager on all PEs
    dMProxy.initMemLog(memLogFile, param.bGpuMemLogger, CkCallbackResumeThread());
    // Implicit wait for all PEs to finish DataManager::initMemLog happens here

    // PE 0 creates/truncates the file and writes a header AFTER collective operation completes
    if (CkMyPe() == 0) {
        FILE* fpLog = CmiFopen(memLogFile.c_str(), "w");
        fprintf(fpLog, "# ChaNGa Memory Log v1.1\n");
        fprintf(fpLog, "# NodeID OpType Size Address Timestamp File:Line PointerID FunctionTag\n");
        int close_err = CmiFclose(fpLog);
        if (close_err != 0) {
            CkPrintf("WARNING: PE 0 failed to close memlog file: %s (Error %d)\n", memLogFile.c_str(), close_err);
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

    // Sequential node flushing: ensures ordered writes to avoid file corruption
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
        CkPrintf("WARNING: Could not open memlog file '%s' for appending.\n", fileName.c_str());
        return; 
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

        // Format: NodeID OpType Size Address Timestamp File:Line PointerID FunctionTag
        // Use %d for NodeID, %zu for size_t, %p for pointer (address), %.6f for timestamp
        // Assumes strings do not contain problematic characters (spaces, newlines)
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
#endif