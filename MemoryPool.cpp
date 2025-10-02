#ifdef CUDA

#include <cuda_runtime.h>
#include "MemoryPool.h"
#include "ParallelGravity.h"
#include "DataManager.h"
#include "memlog.h"
#include "ParallelGravity.decl.h"
#include <unordered_map>
#include <vector>

/**
 * @brief Initialize the GPU memory pool
 * 
 * Configures the default device memory pool to retain freed memory,
 * minimizing allocations from the OS and reducing allocation overhead.
 */
void poolInit() {
    int dev = 0;
    cudaGetDevice(&dev);
    
    cudaMemPool_t pool = nullptr;
    cudaDeviceGetDefaultMemPool(&pool, dev);
    
    // Set threshold to maximum to keep memory cached
    unsigned long long threshold = ~0ull;
    cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &threshold);
}

/**
 * @brief Warm up GPU pool with typical ChaNGa allocation patterns
 * Based on profiling data showing common allocation sizes
 * Aggressive warmup to fully utilize available GPU memory
 */
void poolInitWithWarmup() {
    poolInit();
    
    // Allocation sizes based on typical ChaNGa workload patterns
    // Extended range from 64 KB to 1 GB for comprehensive coverage
    size_t sizes[] = {
        64 * 1024,          // 64 KB - small bucket metadata
        128 * 1024,         // 128 KB - larger bucket metadata
        256 * 1024,         // 256 KB - small interaction lists
        512 * 1024,         // 512 KB
        1 * 1024 * 1024,    // 1 MB - medium interaction lists
        2 * 1024 * 1024,    // 2 MB
        4 * 1024 * 1024,    // 4 MB - common interaction lists
        8 * 1024 * 1024,    // 8 MB
        16 * 1024 * 1024,   // 16 MB - large interaction lists
        32 * 1024 * 1024,   // 32 MB - missed nodes/parts
        64 * 1024 * 1024,   // 64 MB - large missed data
        128 * 1024 * 1024,  // 128 MB
        256 * 1024 * 1024,  // 256 MB
        512 * 1024 * 1024,  // 512 MB
        1024 * 1024 * 1024  // 1 GB - very large allocations
    };
    
    // Aggressive pre-allocation: 10x blocks per size
    // Total warmup will pre-allocate and free significant memory
    // to fully populate the pool's free-lists
    int counts[] = {
        80,  // 64 KB - very common (bucket metadata)
        80,  // 128 KB
        60,  // 256 KB
        60,  // 512 KB
        60,  // 1 MB - common
        60,  // 2 MB
        40,  // 4 MB - frequently used
        40,  // 8 MB
        30,  // 16 MB
        20,  // 32 MB
        20,  // 64 MB
        15,  // 128 MB
        10,  // 256 MB
        8,   // 512 MB
        5    // 1 GB
    };
    
    int numSizes = sizeof(sizes) / sizeof(sizes[0]);
    
    cudaError_t result = poolWarmup(sizes, counts, numSizes, 0);
    if (result == cudaSuccess) {
        CkPrintf("GPU memory pool warmed up\n");
    } else {
        CkPrintf("WARNING: GPU memory pool warmup failed: %s\n", cudaGetErrorString(result));
    }
}

/**
 * @brief Pre-warm the GPU memory pool
 * 
 * Allocates and immediately frees blocks to populate the pool's free-list.
 * This avoids allocation overhead during the critical computation path.
 * 
 * @param sizes Array of allocation sizes in bytes
 * @param counts Array of how many allocations per size
 * @param numSizes Number of different sizes
 * @param stream CUDA stream to use (0 for default)
 */
cudaError_t poolWarmup(const size_t* sizes, const int* counts, int numSizes, cudaStream_t stream) {
    if (sizes == nullptr || counts == nullptr || numSizes <= 0) {
        return cudaErrorInvalidValue;
    }
    
    // Use default stream if none provided
    if (stream == nullptr) {
        cudaStreamCreate(&stream);
    }
    
    // Allocate and free blocks to populate the pool
    for (int i = 0; i < numSizes; i++) {
        size_t size = sizes[i];
        int count = counts[i];
        
        for (int j = 0; j < count; j++) {
            void* ptr = nullptr;
            
            // Allocate block
            cudaError_t result = cudaMallocAsync(&ptr, size, stream);
            if (result != cudaSuccess) {
                // Warmup failure is not critical, just return
                return result;
            }
            
            // Immediately free it back to the pool
            cudaFreeAsync(ptr, stream);
        }
    }
    
    // Synchronize to ensure all warmup operations complete
    cudaStreamSynchronize(stream);
    
    return cudaSuccess;
}

/**
 * @brief Allocate GPU memory from the pool (stream-aware)
 * 
 * Uses cudaMallocAsync for stream-ordered allocation when a stream is provided,
 * falling back to synchronous cudaMalloc when stream is NULL.
 * 
 * @param devPtr Output pointer to allocated device memory
 * @param size Size of memory to allocate in bytes
 * @param stream CUDA stream for ordering (NULL for synchronous)
 * @return cudaError_t result from allocation
 */
cudaError_t poolMalloc(void** devPtr, size_t size, cudaStream_t stream) {
    if (devPtr == nullptr) {
        return cudaErrorInvalidValue;
    }
    
    if (size == 0) {
        *devPtr = nullptr;
        return cudaSuccess;
    }
    
    if (stream != nullptr) {
        // Stream-ordered async allocation
        return cudaMallocAsync(devPtr, size, stream);
    }
    
    // Fallback to synchronous allocation
    return cudaMalloc(devPtr, size);
}

/**
 * @brief Free GPU memory back to the pool (stream-aware)
 * 
 * Uses cudaFreeAsync for stream-ordered deallocation when a stream is provided,
 * falling back to synchronous cudaFree when stream is NULL.
 * 
 * @param ptr Device pointer to free
 * @param stream CUDA stream for ordering (NULL for synchronous)
 * @return cudaError_t result from free operation
 */
cudaError_t poolFree(void* ptr, cudaStream_t stream) {
    if (ptr == nullptr) {
        return cudaSuccess;
    }
    
    if (stream != nullptr) {
        // Stream-ordered async free
        return cudaFreeAsync(ptr, stream);
    }
    
    // Fallback to synchronous free
    return cudaFree(ptr);
}

/**
 * @brief Destroy the memory pool
 * 
 * Cleans up host pool allocations.
 */
void poolDestroy() {
    hostMemoryPoolCleanup();
}

// ============================================================================
// Host Memory Pool (Pinned Memory)
// ============================================================================

namespace {
    static CmiNodeLock g_hostPoolLock = nullptr;
    static std::unordered_map<size_t, std::vector<void*>> g_hostFreePool;
    static std::unordered_map<void*, size_t> g_hostPointerSizes;
    
    /**
     * @brief Round up size to appropriate bucket using tiered granularity
     * Small allocations use fine-grained buckets to reduce waste
     */
    static inline size_t hostBucketSize(size_t requested) {
        if (requested == 0) return 0;
        
        // Fine-grained buckets for small allocations
        if (requested <= 1024) return 1024;              // ≤1KB → 1KB
        if (requested <= 2048) return 2048;              // ≤2KB → 2KB
        if (requested <= 4096) return 4096;              // ≤4KB → 4KB
        if (requested <= 8 * 1024) return 8 * 1024;      // ≤8KB → 8KB
        if (requested <= 16 * 1024) return 16 * 1024;    // ≤16KB → 16KB
        if (requested <= 32 * 1024) return 32 * 1024;    // ≤32KB → 32KB
        
        // Medium allocations: 64KB granularity
        const size_t step = 64u * 1024u;
        return ((requested + step - 1) / step) * step;
    }
}

/**
 * @brief Initialize host memory pool (lazy init, no-op)
 */
cudaError_t hostMemoryPoolInit(size_t /*initialSize*/) {
    return cudaSuccess;
}

/**
 * @brief Warm up host pool with typical ChaNGa allocation patterns
 * 
 * Pre-allocates common sizes based on typical buffer allocations:
 * - Interaction list buffers (1-64 MB)
 * - Moment and particle transfer buffers (1-128 MB)
 * - Bucket metadata (small allocations)
 */
void hostPoolInitWithWarmup() {
    // Lazy init the lock first
    if (g_hostPoolLock == nullptr) {
        g_hostPoolLock = CmiCreateLock();
    }
    
    // Comprehensive host allocation sizes based on profiled ChaNGa workloads
    // Covering tiny (1KB) to very large (640MB) allocations
    // Sized to prevent ANY pool exhaustion and host thread blocking
    size_t sizes[] = {
        1024,            // 1 KB - tiny allocations
        2048,            // 2 KB
        4096,            // 4 KB
        8 * 1024,        // 8 KB
        16 * 1024,       // 16 KB
        32 * 1024,       // 32 KB
        64 * 1024,       // 64 KB - very common small buffers
        128 * 1024,      // 128 KB - very common
        256 * 1024,      // 256 KB
        384 * 1024,      // 384 KB - gap filler
        512 * 1024,      // 512 KB
        640 * 1024,      // 640 KB - observed
        768 * 1024,      // 768 KB - gap filler
        1 * 1024 * 1024, // 1 MB
        1536 * 1024,     // 1.5 MB - observed
        2 * 1024 * 1024, // 2 MB
        2560 * 1024,     // 2.5 MB - high frequency in large runs
        3 * 1024 * 1024, // 3 MB - high frequency in large runs
        4 * 1024 * 1024, // 4 MB
        5 * 1024 * 1024, // 5 MB - very high frequency in large runs
        6 * 1024 * 1024, // 6 MB - observed
        8 * 1024 * 1024, // 8 MB
        10 * 1024 * 1024,// 10 MB - very high frequency in large runs
        12 * 1024 * 1024,// 12 MB - gap filler
        16 * 1024 * 1024,// 16 MB - high frequency
        20 * 1024 * 1024,// 20 MB
        24 * 1024 * 1024,// 24 MB
        28 * 1024 * 1024,// 28 MB
        32 * 1024 * 1024,// 32 MB - very high frequency
        40 * 1024 * 1024,// 40 MB
        48 * 1024 * 1024,// 48 MB
        56 * 1024 * 1024,// 56 MB
        64 * 1024 * 1024,// 64 MB - high frequency in large runs
        80 * 1024 * 1024,// 80 MB
        96 * 1024 * 1024,// 96 MB
        128 * 1024 * 1024,// 128 MB
        160 * 1024 * 1024,// 160 MB
        192 * 1024 * 1024,// 192 MB - observed in very large runs
        200 * 1024 * 1024,// 200 MB - observed (~190-200MB range)
        256 * 1024 * 1024,// 256 MB
        512 * 1024 * 1024,// 512 MB - observed (~490-520MB range)
        640 * 1024 * 1024 // 640 MB - observed (~590-630MB range)
    };
    
    // ULTRA-AGGRESSIVE pre-allocation: 3x observed usage for all buckets <64MB
    // Ensuring ZERO host thread blocking during execution
    // Based on large-scale run profiling with pool exhaustion data
    int counts[] = {
        1000,  // 1 KB - small metadata allocations (+50 safety)
        1000,  // 2 KB (+50 safety)
        1000,  // 4 KB (+50 safety, observed 3 exhaustions)
        1000,  // 8 KB (+50 safety, observed 6 exhaustions)
        1000,  // 16 KB - had 500, exhausted +39, now 700 (40% margin)
        1100, // 32 KB - had 800, exhausted +114, now 1100 (20% margin)
        4400, // 64 KB - had 3500, exhausted +135, now 4400 (20% margin)
        5600, // 128 KB - had 4500, exhausted +148, now 5600 (20% margin)
        100,   // 256 KB - (+10 safety)
        100,   // 384 KB - (+5 safety, observed 1 exhaustion)
        100,   // 512 KB - (+5 safety, observed 2 exhaustions)
        100,   // 640 KB - (+5 safety, observed 1 exhaustion)
        100,   // 768 KB - (+5 safety, observed 6 exhaustions)
        100,   // 1 MB - had 50, exhausted +8, now 65 (30% margin)
        100,   // 1.5 MB - had 50, exhausted +10, now 75 (50% margin)
        140,  // 2 MB - had 100, exhausted +19, now 140 (20% margin)
        140,  // 2.5 MB - had 100, exhausted +18, now 140 (20% margin)
        110,  // 3 MB - had 80, exhausted +22, now 110 (10% margin)
        190,  // 4 MB - had 120, exhausted +31, now 190 (25% margin)
        1730, // 5 MB - had 1400, exhausted +30, now 1730 (20% margin)
        75,   // 6 MB - had 60, exhausted +11, now 75 (25% margin)
        240,  // 8 MB - had 160, exhausted +32, now 240 (25% margin)
        2000, // 10 MB - had 1600, exhausted +38, now 2000 (25% margin)
        15,   // 12 MB - (+5 safety, observed 1 exhaustion)
        2200, // 16 MB - had 1750, exhausted +44, now 2200 (25% margin)
        280,  // 20 MB - had 180, exhausted +44, now 280 (25% margin)
        35,   // 24 MB - (+5 safety)
        35,   // 28 MB - (+5 safety)
        2800, // 32 MB - had 2250, exhausted +66, now 2800 (25% margin)
        35,   // 40 MB - (+5 safety)
        20,   // 48 MB - (+5 safety)
        15,   // 56 MB - (+5 safety)
        2150, // 64 MB - had 1700, exhausted +51, now 2150 (25% margin)
        15,   // 80 MB
        10,   // 96 MB
        15,   // 128 MB
        10,   // 160 MB
        10,   // 192 MB
        10,   // 200 MB
        10,   // 256 MB
        10,   // 512 MB
        8     // 640 MB
    };
    
    int numSizes = sizeof(sizes) / sizeof(sizes[0]);
    
    // Warm up the pool
    for (int i = 0; i < numSizes; i++) {
        for (int j = 0; j < counts[i]; j++) {
            void* ptr = nullptr;
            cudaError_t res = hostMallocPool(&ptr, sizes[i]);
            if (res != cudaSuccess) {
                CkPrintf("WARNING: Host pool warmup failed for size %zu: %s\n", 
                         sizes[i], cudaGetErrorString(res));
                return;
            }
            // Immediately free back to pool
            hostFreePool(ptr);
        }
    }
    
    CkPrintf("Host memory pool warmed up (1KB - 640MB, 42 bucket sizes, ~25000 blocks)\n");
}

/**
 * @brief Allocate pinned host memory from the pool
 * 
 * Uses bucketed free-list for memory reuse. Bucket size is 64KB granularity.
 */
cudaError_t hostMallocPool(void** ptr, size_t size) {
    if (ptr == nullptr) return cudaErrorInvalidValue;
    if (size == 0) { *ptr = nullptr; return cudaSuccess; }
    
    // Lazy initialization of lock
    if (g_hostPoolLock == nullptr) {
        g_hostPoolLock = CmiCreateLock();
    }
    
    size_t bucket = hostBucketSize(size);
    
    // Try to reuse from pool
    CmiLock(g_hostPoolLock);
    auto it = g_hostFreePool.find(bucket);
    if (it != g_hostFreePool.end() && !it->second.empty()) {
        void* p = it->second.back();
        it->second.pop_back();
        g_hostPointerSizes[p] = bucket;
        CmiUnlock(g_hostPoolLock);
        *ptr = p;
        // DEBUG: Pool hit
        static int hitCount = 0;
        if (++hitCount % 100 == 0) {
            CkPrintf("[HostPool] HIT %d: reused %zu byte block (bucket %zu)\n", hitCount, size, bucket);
        }
        return cudaSuccess;
    }
    CmiUnlock(g_hostPoolLock);
    
    // Allocate new pinned block at bucket size
    // DEBUG: Pool miss
    static int missCount = 0;
    if (++missCount % 100 == 0) {
        CkPrintf("[HostPool] MISS %d: allocating new %zu byte block (bucket %zu)\n", missCount, size, bucket);
    }
    
    double startTime = CkWallTimer();
    void* newPtr = nullptr;
    cudaError_t res = cudaMallocHost(&newPtr, bucket);
    double allocTime = CkWallTimer() - startTime;
    
    if (allocTime > 0.001) {  // > 1ms
        CkPrintf("[HostPool] WARNING: cudaMallocHost took %.3f ms for %zu bytes\n", allocTime * 1000, bucket);
    }
    
    if (res != cudaSuccess) return res;
    
    CmiLock(g_hostPoolLock);
    g_hostPointerSizes[newPtr] = bucket;
    CmiUnlock(g_hostPoolLock);
    *ptr = newPtr;
    return cudaSuccess;
}

/**
 * @brief Stream-aware host allocation (currently ignores stream)
 */
cudaError_t hostMallocPoolStream(void** ptr, size_t size, cudaStream_t /*stream*/) {
    return hostMallocPool(ptr, size);
}

/**
 * @brief Free pinned host memory back to the pool
 */
cudaError_t hostFreePool(void* ptr) {
    if (ptr == nullptr) return cudaSuccess;
    
    if (g_hostPoolLock == nullptr) {
        // Pool not initialized; free directly
        return cudaFreeHost(ptr);
    }
    
    CmiLock(g_hostPoolLock);
    auto it = g_hostPointerSizes.find(ptr);
    if (it == g_hostPointerSizes.end()) {
        CmiUnlock(g_hostPoolLock);
        // Unknown pointer - free directly to be safe
        return cudaFreeHost(ptr);
    }
    
    size_t bucket = it->second;
    g_hostPointerSizes.erase(it);
    g_hostFreePool[bucket].push_back(ptr);
    CmiUnlock(g_hostPoolLock);
    return cudaSuccess;
}

/**
 * @brief Stream-aware host free (currently ignores stream)
 */
cudaError_t hostFreePoolStream(void* ptr, cudaStream_t /*stream*/) {
    return hostFreePool(ptr);
}

// ============================================================================
// Host Memory Allocation with Logging
// ============================================================================

/**
 * @brief Host memory allocation with optional logging
 * 
 * Allocates from the host pool and logs the event if logging is enabled.
 */
cudaError_t hostMalloc(void** ptr, size_t size,
                       const char* tag, const char* functionTag,
                       const char* file, int line) {
    double timestamp = CkWallTimer();
    cudaError_t result = hostMallocPool(ptr, size);
    double timestamp_after = CkWallTimer();
    
    // Log if enabled
    DataManager* dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    if (dm && dm->cpuMemLog && dm->bCpuMemLogger) {
        MemLogOpType opType = (result == cudaSuccess) ? MEMLOG_ALLOC : MEMLOG_ALLOC_FAIL;
        uintptr_t address = (result == cudaSuccess && *ptr != NULL) ? (uintptr_t)(*ptr) : 0;
        MemLogEvent event(CkMyNode(), opType, size, address, timestamp_after, file, line, tag, functionTag);
        
        CmiLock(dm->lockCpuMemLog);
        dm->cpuMemLog->meTab.push_back(event);
        CmiUnlock(dm->lockCpuMemLog);
    }
    
    return result;
}

/**
 * @brief Host memory deallocation with optional logging
 * 
 * Frees to the pool and logs the event if logging is enabled.
 */
cudaError_t hostFree(void* ptr,
                     const char* tag, const char* functionTag,
                     const char* file, int line) {
    double timestamp = CkWallTimer();
    MemLogOpType opType;
    uintptr_t address = (uintptr_t)ptr;
    
    // Access node-local DataManager
    DataManager* dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    
    if (ptr == NULL) {
        // NULL pointer free is a no-op
        opType = MEMLOG_FREE_SKIP;
        if (dm && dm->cpuMemLog && dm->bCpuMemLogger) {
            MemLogEvent event(CkMyNode(), opType, 0, address, timestamp, file, line, tag, functionTag);
            CmiLock(dm->lockCpuMemLog);
            dm->cpuMemLog->meTab.push_back(event);
            CmiUnlock(dm->lockCpuMemLog);
        }
        return cudaSuccess;
    }
    
    cudaError_t result = hostFreePool(ptr);
    double timestamp_after = CkWallTimer();
    
    // Determine operation type based on free success/failure
    opType = (result == cudaSuccess) ? MEMLOG_FREE : MEMLOG_FREE_FAIL;
    
    if (dm && dm->cpuMemLog && dm->bCpuMemLogger) {
        MemLogEvent event(CkMyNode(), opType, 0, address, timestamp_after, file, line, tag, functionTag);
        CmiLock(dm->lockCpuMemLog);
        dm->cpuMemLog->meTab.push_back(event);
        CmiUnlock(dm->lockCpuMemLog);
    }
    
    return result;
}

/**
 * @brief Cleanup all host pool allocations
 */
cudaError_t hostMemoryPoolCleanup() {
    if (g_hostPoolLock == nullptr) return cudaSuccess;
    
    CmiLock(g_hostPoolLock);
    for (auto &kv : g_hostFreePool) {
        for (void* p : kv.second) {
            cudaFreeHost(p);
        }
    }
    g_hostFreePool.clear();
    g_hostPointerSizes.clear();
    CmiUnlock(g_hostPoolLock);
    return cudaSuccess;
}

// ============================================================================
// GPU Memory Allocation with Logging
// ============================================================================

/**
 * @brief GPU memory allocation with optional logging
 * 
 * Allocates from the pool and logs the event if logging is enabled.
 */
cudaError_t gpuMalloc(void** devPtr, size_t size, cudaStream_t stream,
                      const char* tag, const char* functionTag,
                      const char* file, int line) {
    double timestamp = CkWallTimer();
    cudaError_t result = poolMalloc(devPtr, size, stream);
    double timestamp_after = CkWallTimer();
    
    // Log if enabled
    DataManager* dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    if (dm && dm->memLog && dm->bGpuMemLogger) {
        MemLogOpType opType = (result == cudaSuccess) ? MEMLOG_ALLOC : MEMLOG_ALLOC_FAIL;
        uintptr_t address = (result == cudaSuccess && *devPtr != NULL) ? (uintptr_t)(*devPtr) : 0;
        MemLogEvent event(CkMyNode(), opType, size, address, timestamp_after, file, line, tag, functionTag);
        
        CmiLock(dm->lockMemLog);
        dm->memLog->meTab.push_back(event);
        CmiUnlock(dm->lockMemLog);
    }
    
    return result;
}

/**
 * @brief GPU memory deallocation with optional logging
 * 
 * Frees to the pool and logs the event if logging is enabled.
 */
cudaError_t gpuFree(void* devPtr, cudaStream_t stream,
                    const char* tag, const char* functionTag,
                    const char* file, int line) {
    double timestamp = CkWallTimer();
    MemLogOpType opType;
    uintptr_t address = (uintptr_t)devPtr;
    
    // Access node-local DataManager
    DataManager* dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    
    if (devPtr == NULL) {
        // NULL pointer free is a no-op
        opType = MEMLOG_FREE_SKIP;
        if (dm && dm->memLog && dm->bGpuMemLogger) {
            MemLogEvent event(CkMyNode(), opType, 0, address, timestamp, file, line, tag, functionTag);
            CmiLock(dm->lockMemLog);
            dm->memLog->meTab.push_back(event);
            CmiUnlock(dm->lockMemLog);
        }
        return cudaSuccess;
    }
    
    cudaError_t result = poolFree(devPtr, stream);
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

// ============================================================================
// Memory Log Management
// ============================================================================

/**
 * @brief Initialize memory log on DataManager
 * Called collectively across all DataManagers
 */
void DataManager::initMemLog(std::string _fileName, int bGpuMemLoggerFlag, const CkCallback &cb) {
    CmiLock(lockMemLog);
    if (memLog != nullptr) { 
        memLog->fileName = _fileName;
    } else {
        CkPrintf("WARNING PE %d: memLog is NULL in initMemLog! Cannot set filename.\n", CkMyPe());
    }
    bGpuMemLogger = bGpuMemLoggerFlag;
    CmiUnlock(lockMemLog);
    contribute(cb);
}

/**
 * @brief Initialize CPU memory log on DataManager
 * Called collectively across all DataManagers
 */
void DataManager::initCpuMemLog(std::string _fileName, int bCpuMemLoggerFlag, const CkCallback &cb) {
    CmiLock(lockCpuMemLog);
    if (cpuMemLog != nullptr) { 
        cpuMemLog->fileName = _fileName;
    } else {
        CkPrintf("WARNING PE %d: cpuMemLog is NULL in initCpuMemLog! Cannot set filename.\n", CkMyPe());
    }
    bCpuMemLogger = bCpuMemLoggerFlag;
    CmiUnlock(lockCpuMemLog);
    contribute(cb);
}

/**
 * @brief Initialize memory log file on PE 0
 * Called from Main chare
 */
void Main::initMemLog() {
    std::string memLogFile = "memlog.out";
    
    // Send filename and flag to all DataManagers and wait
    dMProxy.initMemLog(memLogFile, param.bGpuMemLogger, CkCallbackResumeThread());
    
    // PE 0 creates the file and writes header
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

/**
 * @brief Initialize CPU memory log file on PE 0
 * Called from Main chare
 */
void Main::initCpuMemLog() {
    std::string cpuMemLogFile = "cpumemlog.out";
    
    // Send filename and flag to all DataManagers and wait
    dMProxy.initCpuMemLog(cpuMemLogFile, param.bCpuMemLogger, CkCallbackResumeThread());
    
    // PE 0 creates the file and writes header
    if (CkMyPe() == 0) {
        FILE* fpLog = CmiFopen(cpuMemLogFile.c_str(), "w");
        fprintf(fpLog, "# ChaNGa CPU Memory Log v1.0\n");
        fprintf(fpLog, "# NodeID OpType Size Address Timestamp File:Line PointerID FunctionTag\n");
        int close_err = CmiFclose(fpLog);
        if (close_err != 0) {
            CkPrintf("WARNING: PE 0 failed to close CPU memlog file: %s (Error %d)\n", cpuMemLogFile.c_str(), close_err);
        }
    }
}

/**
 * @brief Flush memlog sequentially across nodes
 */
void DataManager::flushMemLog(const CkCallback& cb) {
    if (memLog) {
        memLog->flush();
    } else {
        CkPrintf("WARNING Node %d: memLog is NULL in flushMemLog! Skipping flush.\n", thisIndex);
    }
    
    // Sequential flushing to avoid file corruption
    if(thisIndex != CkNumNodes()-1) {
        thisProxy[thisIndex + 1].flushMemLog(cb);
    } else {
        cb.send();
    }
}

/**
 * @brief Flush CPU memlog sequentially across nodes
 */
void DataManager::flushCpuMemLog(const CkCallback& cb) {
    if (cpuMemLog) {
        cpuMemLog->flush();
    } else {
        CkPrintf("WARNING Node %d: cpuMemLog is NULL in flushCpuMemLog! Skipping flush.\n", thisIndex);
    }
    
    // Sequential flushing to avoid file corruption
    if(thisIndex != CkNumNodes()-1) {
        thisProxy[thisIndex + 1].flushCpuMemLog(cb);
    } else {
        cb.send();
    }
}

/**
 * @brief Write buffered memory log events to file
 */
void MemLog::flush() {
    if (meTab.empty()) {
        return;
    }
    
    FILE* outfile = CmiFopen(fileName.c_str(), "a");
    if (outfile == NULL) {
        CkPrintf("WARNING: Could not open memlog file '%s' for appending.\n", fileName.c_str());
        return;
    }
    
    // Write all buffered events
    for (const auto& event : meTab) {
        const char* opTypeStr;
        switch (event.opType) {
            case MEMLOG_ALLOC:      opTypeStr = "ALLOC";      break;
            case MEMLOG_FREE:       opTypeStr = "FREE ";      break;
            case MEMLOG_ALLOC_FAIL: opTypeStr = "ALLOC_F";    break;
            case MEMLOG_FREE_FAIL:  opTypeStr = "FREE_F ";    break;
            case MEMLOG_FREE_SKIP:  opTypeStr = "FREE_S ";    break;
            default:                opTypeStr = "UNKNOWN";    break;
        }
        
        fprintf(outfile, "%d %s %zu %p %.6f %s %s %s\n",
                event.nodeId,
                opTypeStr,
                event.size,
                (void*)event.address,
                event.timestamp,
                event.location.c_str(),
                event.pointerId.c_str(),
                event.functionTag.c_str());
    }
    
    int result = CmiFclose(outfile);
    if (result != 0) {
        CkPrintf("WARNING: Failed to close memlog file '%s' properly (Error %d).\n", fileName.c_str(), result);
    }
    
    meTab.clear();
}

#endif // CUDA
