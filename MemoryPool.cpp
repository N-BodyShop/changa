#ifdef CUDA

#include <cuda_runtime.h>
#include "MemoryPool.h"
#include "ParallelGravity.h"
#include "DataManager.h"
#include "memlog.h"
#include "ParallelGravity.decl.h"
#include <unordered_map>
#include <vector>

// ============================================================================
// GPU MEMORY POOL
// ============================================================================

// ----------------------------------------------------------------------------
// GPU Pool Initialization & Warmup
// ----------------------------------------------------------------------------

/**
 * @brief Initialize the GPU memory pool
 * 
 * Configures the default device memory pool with a 5GB limit to prevent
 * excessive memory consumption when multiple ranks share a GPU device.
 */
void gpuPoolInit() {
    int dev = 0;
    cudaGetDevice(&dev);
    
    cudaMemPool_t pool = nullptr;
    cudaDeviceGetDefaultMemPool(&pool, dev);
    
    // Set 5GB limit (5 * 1024^3 bytes) to prevent excessive memory usage
    // with multiple ranks per GPU device
    unsigned long long threshold = 5ULL * 1024 * 1024 * 1024;
    cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &threshold);

    // Allocation sizes based on typical ChaNGa workload patterns
    // Extended range from 64 KB to 1 GB for comprehensive coverage
    size_t sizes[] = {
        64 * 1024,          // 64 KB 
        128 * 1024,         // 128 KB 
        256 * 1024,         // 256 KB 
        512 * 1024,         // 512 KB
        1 * 1024 * 1024,    // 1 MB 
        2 * 1024 * 1024,    // 2 MB
        4 * 1024 * 1024,    // 4 MB 
        8 * 1024 * 1024,    // 8 MB
        16 * 1024 * 1024,   // 16 MB 
        32 * 1024 * 1024,   // 32 MB 
        64 * 1024 * 1024,   // 64 MB 
        128 * 1024 * 1024,  // 128 MB
        256 * 1024 * 1024,  // 256 MB
        512 * 1024 * 1024,  // 512 MB
        1024 * 1024 * 1024  // 1 GB 
    };
    
    // Total warmup will pre-allocate and free significant memory
    // to fully populate the pool's free-lists
    int counts[] = {
        80,  // 64 KB 
        80,  // 128 KB
        60,  // 256 KB
        60,  // 512 KB
        60,  // 1 MB 
        60,  // 2 MB
        40,  // 4 MB 
        40,  // 8 MB
        30,  // 16 MB
        20,  // 32 MB
        10,  // 64 MB
        5,   // 128 MB
        5,   // 256 MB
        2,   // 512 MB
        1    // 1 GB
    };
    
    int numSizes = sizeof(sizes) / sizeof(sizes[0]);
    
    cudaError_t result = gpuPoolWarmup(sizes, counts, numSizes, 0);
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
 * @return cudaError_t from allocation operations, or cudaSuccess
 */
cudaError_t gpuPoolWarmup(const size_t* sizes, const int* counts, int numSizes, cudaStream_t stream) {
    if (sizes == nullptr || counts == nullptr || numSizes <= 0) {
        return cudaErrorInvalidValue;
    }
    
    // Use default stream if none provided
    bool createdStream = false;
    if (stream == nullptr) {
        cudaError_t createResult = cudaStreamCreate(&stream);
        if (createResult != cudaSuccess) {
            return createResult;
        }
        createdStream = true;
    }
    
    // Allocate and free blocks to populate the pool
    cudaError_t finalResult = cudaSuccess;
    for (int i = 0; i < numSizes; i++) {
        size_t size = sizes[i];
        int count = counts[i];
        
        for (int j = 0; j < count; j++) {
            void* ptr = nullptr;
            
            // Allocate block
            cudaError_t result = cudaMallocAsync(&ptr, size, stream);
            if (result != cudaSuccess) {
                finalResult = result;
                break;
            }
            
            // Immediately free it back to the pool
            cudaFreeAsync(ptr, stream);
        }
        
        if (finalResult != cudaSuccess) {
            break;
        }
    }
    
    // Synchronize to ensure all warmup operations complete
    if (finalResult == cudaSuccess) {
        finalResult = cudaStreamSynchronize(stream);
    }
    
    // Clean up stream if we created it
    if (createdStream) {
        cudaStreamDestroy(stream);
    }
    
    return finalResult;
}

// ----------------------------------------------------------------------------
// GPU Pool Raw Operations (Internal - use macros for public API)
// ----------------------------------------------------------------------------

/**
 * @brief Allocate GPU memory from the pool (raw, internal, stream-aware)
 * 
 * Uses cudaMallocAsync for stream-ordered allocation when a stream is provided,
 * falling back to synchronous cudaMalloc when stream is NULL.
 * 
 * @param devPtr Output pointer to allocated device memory
 * @param size Size of memory to allocate in bytes
 * @param stream CUDA stream for ordering (NULL for synchronous)
 * @return cudaError_t result from allocation
 */
cudaError_t gpuPoolMallocRaw(void** devPtr, size_t size, cudaStream_t stream) {
    if (devPtr == nullptr) {
        return cudaErrorInvalidValue;
    }
    
    // Initialize to NULL in case of failure
    *devPtr = nullptr;
    
    if (size == 0) {
        return cudaSuccess;
    }
    
    cudaError_t result;
    if (stream != nullptr) {
        // Stream-ordered async allocation
        result = cudaMallocAsync(devPtr, size, stream);
    } else {
        // Fallback to synchronous allocation
        result = cudaMalloc(devPtr, size);
    }
    
    // Ensure pointer is NULL on failure
    if (result != cudaSuccess) {
        *devPtr = nullptr;
    }
    
    return result;
}

/**
 * @brief Free GPU memory back to the pool (raw, internal, stream-aware)
 * 
 * Uses cudaFreeAsync for stream-ordered deallocation when a stream is provided,
 * falling back to synchronous cudaFree when stream is NULL.
 * 
 * @param ptr Device pointer to free
 * @param stream CUDA stream for ordering (NULL for synchronous)
 * @return cudaError_t result from free operation
 */
cudaError_t gpuPoolFreeRaw(void* ptr, cudaStream_t stream) {
    if (ptr == nullptr) {
        return cudaSuccess;
    }
    
    if (stream != nullptr) {
        // Stream-ordered async free
        cudaError_t result = cudaFreeAsync(ptr, stream);
        // Synchronize to ensure operations complete before memory is returned to pool.
        // Without this, cudaMallocAsync can return the same memory while operations
        // using it are still in flight, causing race conditions and non-determinism.
        if (result == cudaSuccess) {
            cudaStreamSynchronize(stream);
        }
        return result;
    }
    
    // Fallback to synchronous free
    return cudaFree(ptr);
}

// ----------------------------------------------------------------------------
// GPU Pool Logging Helpers (Private)
// ----------------------------------------------------------------------------

namespace {
    /**
     * @brief Log GPU allocation event if logging is enabled
     * 
     * @param result Result code from allocation operation
     * @param devPtr Pointer to allocated device memory (or NULL on failure)
     * @param size Size of allocation in bytes
     * @param timestamp Wall-clock timestamp of operation
     * @param file Source file where allocation occurred
     * @param line Line number where allocation occurred
     * @param tag Pointer identifier tag (e.g. "devBuf")
     * @param functionTag Function context tag
     */
    inline void logGpuAlloc(cudaError_t result, void** devPtr, size_t size,
                            double timestamp, const char* file, int line,
                            const char* tag, const char* functionTag) {
        DataManager* dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
        if (dm && dm->memLog && dm->bGpuMemLogger) {
            MemLogOpType opType = (result == cudaSuccess) ? MEMLOG_ALLOC : MEMLOG_ALLOC_FAIL;
            uintptr_t address = (result == cudaSuccess && *devPtr != NULL) ? (uintptr_t)(*devPtr) : 0;
            MemLogEvent event(CkMyNode(), opType, size, address, timestamp, file, line, tag, functionTag);
            
            CmiLock(dm->lockMemLog);
            dm->memLog->meTab.push_back(event);
            CmiUnlock(dm->lockMemLog);
        }
    }
    
    /**
     * @brief Log GPU free event if logging is enabled
     * 
     * @param result Result code from free operation
     * @param devPtr Device pointer being freed
     * @param timestamp Wall-clock timestamp of operation
     * @param file Source file where free occurred
     * @param line Line number where free occurred
     * @param tag Pointer identifier tag
     * @param functionTag Function context tag
     */
    inline void logGpuFree(cudaError_t result, void* devPtr, double timestamp,
                           const char* file, int line, const char* tag,
                           const char* functionTag) {
        DataManager* dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
        if (dm && dm->memLog && dm->bGpuMemLogger) {
            MemLogOpType opType;
            if (devPtr == NULL) {
                opType = MEMLOG_FREE_SKIP;
            } else {
                opType = (result == cudaSuccess) ? MEMLOG_FREE : MEMLOG_FREE_FAIL;
            }
            uintptr_t address = (uintptr_t)devPtr;
            MemLogEvent event(CkMyNode(), opType, 0, address, timestamp, file, line, tag, functionTag);
            
            CmiLock(dm->lockMemLog);
            dm->memLog->meTab.push_back(event);
            CmiUnlock(dm->lockMemLog);
        }
    }
}

// ----------------------------------------------------------------------------
// GPU Pool API with Logging - Pool-based Allocation
// ----------------------------------------------------------------------------

/**
 * @brief GPU memory allocation from pool with logging
 * 
 * Allocates from the pool and logs the event if logging is enabled.
 * Use the gpuPoolMalloc macro instead of calling this directly.
 * 
 * @param devPtr Output pointer to allocated device memory
 * @param size Size to allocate in bytes
 * @param stream CUDA stream for ordering
 * @param tag Pointer identifier tag
 * @param functionTag Function context tag
 * @param file Source file
 * @param line Line number
 * @return cudaError_t from allocation operation
 */
cudaError_t gpuPoolMallocImpl(void** devPtr, size_t size, cudaStream_t stream,
                              const char* tag, const char* functionTag,
                              const char* file, int line) {
    cudaError_t result = gpuPoolMallocRaw(devPtr, size, stream);
    logGpuAlloc(result, devPtr, size, CkWallTimer(), file, line, tag, functionTag);
    return result;
}

/**
 * @brief GPU memory deallocation to pool with logging
 * 
 * Frees to the pool and logs the event if logging is enabled.
 * Use the gpuPoolFree macro instead of calling this directly.
 * 
 * @param devPtr Device pointer to free
 * @param stream CUDA stream for ordering
 * @param tag Pointer identifier tag
 * @param functionTag Function context tag
 * @param file Source file
 * @param line Line number
 * @return cudaError_t from free operation
 */
cudaError_t gpuPoolFreeImpl(void* devPtr, cudaStream_t stream,
                            const char* tag, const char* functionTag,
                            const char* file, int line) {
    if (devPtr == NULL) {
        logGpuFree(cudaSuccess, devPtr, CkWallTimer(), file, line, tag, functionTag);
        return cudaSuccess;
    }
    
    cudaError_t result = gpuPoolFreeRaw(devPtr, stream);
    logGpuFree(result, devPtr, CkWallTimer(), file, line, tag, functionTag);
    return result;
}

// ----------------------------------------------------------------------------
// GPU API with Logging - Direct Allocation (Bypass Pool)
// ----------------------------------------------------------------------------

/**
 * @brief Direct GPU memory allocation with logging (bypasses pool)
 * 
 * Allocates directly using cudaMalloc and logs the event if logging is enabled.
 * Use the gpuMalloc macro instead of calling this directly.
 * 
 * @param devPtr Output pointer to allocated device memory
 * @param size Size to allocate in bytes
 * @param stream CUDA stream (currently unused, reserved for future use)
 * @param tag Pointer identifier tag
 * @param functionTag Function context tag
 * @param file Source file
 * @param line Line number
 * @return cudaError_t from allocation operation
 */
cudaError_t gpuMallocImpl(void** devPtr, size_t size, cudaStream_t stream,
                          const char* tag, const char* functionTag,
                          const char* file, int line) {
    if (devPtr == nullptr) {
        return cudaErrorInvalidValue;
    }
    
    if (size == 0) {
        *devPtr = nullptr;
        return cudaSuccess;
    }
    
    // Use cudaMalloc directly (ignores stream for now)
    cudaError_t result = cudaMalloc(devPtr, size);
    logGpuAlloc(result, devPtr, size, CkWallTimer(), file, line, tag, functionTag);
    return result;
}

/**
 * @brief Direct GPU memory deallocation with logging (bypasses pool)
 * 
 * Frees directly using cudaFree and logs the event if logging is enabled.
 * Use the gpuFree macro instead of calling this directly.
 * 
 * @param devPtr Device pointer to free
 * @param stream CUDA stream (currently unused, reserved for future use)
 * @param tag Pointer identifier tag
 * @param functionTag Function context tag
 * @param file Source file
 * @param line Line number
 * @return cudaError_t from free operation
 */
cudaError_t gpuFreeImpl(void* devPtr, cudaStream_t stream,
                        const char* tag, const char* functionTag,
                        const char* file, int line) {
    if (devPtr == NULL) {
        logGpuFree(cudaSuccess, devPtr, CkWallTimer(), file, line, tag, functionTag);
        return cudaSuccess;
    }
    
    // Direct free bypassing pool
    cudaError_t result = cudaFree(devPtr);
    logGpuFree(result, devPtr, CkWallTimer(), file, line, tag, functionTag);
    return result;
}

// ============================================================================
// Host Memory Pool (Pinned Memory)
// ============================================================================

// ----------------------------------------------------------------------------
// Host Pool Internals (Private State and Helpers)
// ----------------------------------------------------------------------------

namespace {
    // --- Host Pool Global State ---
    
    static CmiNodeLock g_hostPoolLock = nullptr;
    
    // Initialization lock - used to protect lazy initialization of g_hostPoolLock
    // Function-local static ensures thread-safe initialization (C++11+)
    static CmiNodeLock& getInitLock() {
        static CmiNodeLock initLock = CmiCreateLock();
        return initLock;
    }
    static std::unordered_map<size_t, std::vector<void*>> g_hostFreePool;
    static std::unordered_map<void*, size_t> g_hostPointerSizes;
    
    // Adaptive refill tracking
    static std::unordered_map<size_t, int> g_bucketMissCount;
    static int g_timestepCount = 0;
    
    // Per-bucket usage tracking for capacity-based trimming
    struct BucketUsageInfo {
        int lastUsedTimestep = 0;
        int recentHits = 0;  // Hits in last 10 timesteps
    };
    static std::unordered_map<size_t, BucketUsageInfo> g_bucketUsage;
    
    // --- Host Pool Analytics ---
    
    struct PoolAnalytics {
        // Cumulative data (never reset)
        size_t totalAllocated = 0;      // Total memory allocated from OS
        size_t totalFreed = 0;           // Total memory freed back to OS
        size_t peakMemory = 0;          // Peak memory usage
        size_t currentMemory = 0;       // Current memory usage
        int totalAllocations = 0;        // Total allocation calls
        int totalFrees = 0;              // Total free calls
        int totalHits = 0;               // Pool hits
        int totalMisses = 0;            // Pool misses
        int totalFreedToOS = 0;          // Blocks actually freed to OS (not pool)
        double totalAllocTime = 0.0;     // Total time spent in allocations
        double totalFreeTime = 0.0;     // Total time spent in frees
        
        // Per-timestep data (reset each timestep)
        int timestepAllocations = 0;
        int timestepFrees = 0;
        int timestepHits = 0;
        int timestepMisses = 0;
        int timestepFreedToOS = 0;
        
        // Per-bucket cumulative tracking
        std::unordered_map<size_t, int> bucketTotalAllocations;
        std::unordered_map<size_t, int> bucketTotalMisses;  // Cumulative misses per bucket
    };
    static PoolAnalytics g_poolAnalytics;
    
    // Historical tracking for memory leak detection
    struct TimestepSnapshot {
        int timestep;
        size_t totalBuckets;
        size_t totalFreeBlocks;
        size_t totalFreeBytes;
        size_t liveBlocks;
        size_t totalAllocated;
        size_t totalFreed;
        double timestamp;
    };
    static std::vector<TimestepSnapshot> g_timestepHistory;
    
    // --- Host Pool Helper Functions ---
    
    /**
     * @brief Round up size to appropriate bucket using tiered granularity
     * 
     * Tiered approach to minimize bucket proliferation:
     * - Tiny (≤64KB): Fine-grained buckets (256B, 512B, 1KB, 2KB, 4KB, ..., 64KB)
     * - Small (64KB-1MB): 64KB granularity
     * - Medium (1MB-64MB): 1MB granularity  
     * - Large (>64MB): 50MB granularity
     * 
     * This drastically reduces unique bucket count for large allocations
     * while maintaining efficiency for small allocations.
     * 
     * @param requested Requested allocation size in bytes
     * @return Bucket size that will accommodate the requested size
     */
    static inline size_t hostBucketSize(size_t requested) {
        if (requested == 0) return 0;
        
        // Tiny allocations: use 256B and 512B buckets
        if (requested <= 256) return 256;
        if (requested <= 512) return 512;
        if (requested <= 1024) return 1024;
        if (requested <= 2048) return 2048;
        if (requested <= 4 * 1024) return 4 * 1024;
        if (requested <= 8 * 1024) return 8 * 1024;
        if (requested <= 16 * 1024) return 16 * 1024;
        if (requested <= 32 * 1024) return 32 * 1024;
        if (requested <= 64 * 1024) return 64 * 1024;
        
        // Small allocations (64KB - 1MB): 64KB granularity
        if (requested <= 1024 * 1024) {
            const size_t step = 64 * 1024;
            return ((requested + step - 1) / step) * step;
        }
        
        // Medium allocations (1MB - 64MB): 1MB granularity
        if (requested <= 64 * 1024 * 1024) {
            const size_t step = 1024 * 1024;
            return ((requested + step - 1) / step) * step;
        }
        
        // Large allocations (>64MB): 50MB granularity
        // This prevents bucket proliferation for very large buffers
        const size_t step = 50 * 1024 * 1024;
        return ((requested + step - 1) / step) * step;
    }
    
    // --- Host Pool Logging Helpers ---
    
    /**
     * @brief Log host allocation event if logging is enabled
     * 
     * @param result Result code from allocation operation
     * @param ptr Pointer to allocated memory (or NULL on failure)
     * @param size Size of allocation in bytes
     * @param timestamp Wall-clock timestamp of operation
     * @param file Source file where allocation occurred
     * @param line Line number where allocation occurred
     * @param tag Pointer identifier tag
     * @param functionTag Function context tag
     */
    inline void logHostAlloc(cudaError_t result, void** ptr, size_t size,
                             double timestamp, const char* file, int line,
                             const char* tag, const char* functionTag) {
        DataManager* dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
        if (dm && dm->cpuMemLog && dm->bCpuMemLogger) {
            MemLogOpType opType = (result == cudaSuccess) ? MEMLOG_ALLOC : MEMLOG_ALLOC_FAIL;
            uintptr_t address = (result == cudaSuccess && *ptr != NULL) ? (uintptr_t)(*ptr) : 0;
            MemLogEvent event(CkMyNode(), opType, size, address, timestamp, file, line, tag, functionTag);
            
            CmiLock(dm->lockCpuMemLog);
            dm->cpuMemLog->meTab.push_back(event);
            CmiUnlock(dm->lockCpuMemLog);
        }
    }
    
    /**
     * @brief Log host free event if logging is enabled
     * 
     * @param result Result code from free operation
     * @param ptr Pointer being freed
     * @param timestamp Wall-clock timestamp of operation
     * @param file Source file where free occurred
     * @param line Line number where free occurred
     * @param tag Pointer identifier tag
     * @param functionTag Function context tag
     */
    inline void logHostFree(cudaError_t result, void* ptr, double timestamp,
                            const char* file, int line, const char* tag,
                            const char* functionTag) {
        DataManager* dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
        if (dm && dm->cpuMemLog && dm->bCpuMemLogger) {
            MemLogOpType opType;
            if (ptr == NULL) {
                opType = MEMLOG_FREE_SKIP;
            } else {
                opType = (result == cudaSuccess) ? MEMLOG_FREE : MEMLOG_FREE_FAIL;
            }
            uintptr_t address = (uintptr_t)ptr;
            MemLogEvent event(CkMyNode(), opType, 0, address, timestamp, file, line, tag, functionTag);
            
            CmiLock(dm->lockCpuMemLog);
            dm->cpuMemLog->meTab.push_back(event);
            CmiUnlock(dm->lockCpuMemLog);
        }
    }
}

// ----------------------------------------------------------------------------
// Host Pool Initialization & Warmup
// ----------------------------------------------------------------------------

/**
 * @brief Initialize the host memory pool
 * 
 * Creates the pool lock and pre-warms the pool with common allocation sizes
 * ranging from 256B to 64MB across 19 bucket sizes. Total warmup allocates
 * and immediately frees blocks to populate the pool's free-lists.
 */
// Shared warmup configuration - single source of truth
// This configuration is used by:
// 1. hostPoolInit() - for actual warmup allocation
// 2. hostPoolAnalyzeWarmupEffectiveness() - for analysis targets
// 3. hostPoolAdaptiveRefill() - for refill target buckets
// Modify only this array to change warmup behavior everywhere
namespace {
    struct WarmupConfig {
        size_t size;
        int count;
    };
    
    // Warmup configuration - single place to modify bucket sizes and counts
    const WarmupConfig g_warmupConfig[] = {
        {256, 1000},             // 256 B 
        {512, 1500},             // 512 B 
        {1024, 1000},            // 1 KB 
        {2048, 1000},            // 2 KB 
        {4096, 1500},            // 4 KB 
        {8192, 1500},            // 8 KB 
        {16384, 500},            // 16 KB 
        {32768, 500},            // 32 KB 
        {64 * 1024, 500},        // 64 KB 
        {128 * 1024, 800},       // 128 KB 
        {256 * 1024, 300},       // 256 KB 
        {512 * 1024, 300},       // 512 KB 
        {1 * 1024 * 1024, 1},    // 1 MB 
        {2 * 1024 * 1024, 1},    // 2 MB 
        {4 * 1024 * 1024, 1},    // 4 MB 
        {8 * 1024 * 1024, 1},    // 8 MB 
        {16 * 1024 * 1024, 5},   // 16 MB 
        {32 * 1024 * 1024, 1},   // 32 MB 
        {64 * 1024 * 1024, 5}    // 64 MB 
    };
    
    const int g_numWarmupSizes = sizeof(g_warmupConfig) / sizeof(g_warmupConfig[0]);
}

void hostPoolInit() {
    // Initialize lock before concurrent access to avoid race conditions
    // in lazy initialization
    if (g_hostPoolLock == nullptr) {
        g_hostPoolLock = CmiCreateLock();
    }
    
    size_t totalWarmupBytes = 0;
    
    for (int i = 0; i < g_numWarmupSizes; i++) {
        for (int j = 0; j < g_warmupConfig[i].count; j++) {
            void* ptr = nullptr;
            cudaError_t res = hostPoolMallocRaw(&ptr, g_warmupConfig[i].size);
            if (res != cudaSuccess) {
                CkPrintf("WARNING: Host pool warmup failed for size %zu: %s\n", 
                         g_warmupConfig[i].size, cudaGetErrorString(res));
                return;
            }
            hostPoolFreeRaw(ptr);
            totalWarmupBytes += g_warmupConfig[i].size;
        }
    }
    
    CkPrintf("Host memory pool warmed up\n");
}

// ----------------------------------------------------------------------------
// Host Pool Raw Operations (Internal - use macros for public API)
// ----------------------------------------------------------------------------

/**
 * @brief Allocate host memory from pool (raw, internal)
 * 
 * Attempts to reuse a block from the pool. If no block is available,
 * allocates a new pinned memory block from the OS at the bucket size.
 * 
 * @param ptr Output pointer to allocated memory
 * @param size Size to allocate in bytes
 * @return cudaError_t from allocation operation
 */
cudaError_t hostPoolMallocRaw(void** ptr, size_t size) {
    if (ptr == nullptr) return cudaErrorInvalidValue;
    if (size == 0) { *ptr = nullptr; return cudaSuccess; }
    
    // Lazy initialization of lock (fallback if hostPoolInit() not called)
    // Use double-checked locking pattern with proper synchronization
    if (g_hostPoolLock == nullptr) {
        CmiNodeLock& initLock = getInitLock();
        CmiLock(initLock);
        // Double-check after acquiring lock
        if (g_hostPoolLock == nullptr) {
            g_hostPoolLock = CmiCreateLock();
        }
        CmiUnlock(initLock);
    }
    
    size_t bucket = hostBucketSize(size);
    
    // Try to reuse from pool
    CmiLock(g_hostPoolLock);
    auto it = g_hostFreePool.find(bucket);
    if (it != g_hostFreePool.end() && !it->second.empty()) {
        void* p = it->second.back();
        it->second.pop_back();
        g_hostPointerSizes[p] = bucket;
        
        // Update analytics
        g_poolAnalytics.totalHits++;
        g_poolAnalytics.timestepHits++;
        g_poolAnalytics.totalAllocations++;
        g_poolAnalytics.timestepAllocations++;
        g_poolAnalytics.bucketTotalAllocations[bucket]++;
        g_poolAnalytics.currentMemory += bucket;
        if (g_poolAnalytics.currentMemory > g_poolAnalytics.peakMemory) {
            g_poolAnalytics.peakMemory = g_poolAnalytics.currentMemory;
        }
        
        // Update usage tracking for capacity-based trimming
        g_bucketUsage[bucket].lastUsedTimestep = g_timestepCount;
        g_bucketUsage[bucket].recentHits++;
        
        CmiUnlock(g_hostPoolLock);
        *ptr = p;
        return cudaSuccess;
    }
    CmiUnlock(g_hostPoolLock);
    
    // Allocate new pinned block at bucket size (do this unlocked - it's slow)
    double startTime = CkWallTimer();
    void* newPtr = nullptr;
    cudaError_t res = cudaMallocHost(&newPtr, bucket);
    double allocTime = CkWallTimer() - startTime;
    
    if (res != cudaSuccess) return res;
    
    // Re-acquire lock for all map and analytics updates
    // This ensures thread-safety while keeping the slow cudaMallocHost unlocked
    CmiLock(g_hostPoolLock);
    
    // Track miss for adaptive refill
    g_bucketMissCount[bucket]++;
    
    // Update analytics
    g_poolAnalytics.totalMisses++;
    g_poolAnalytics.timestepMisses++;
    g_poolAnalytics.totalAllocations++;
    g_poolAnalytics.timestepAllocations++;
    g_poolAnalytics.bucketTotalAllocations[bucket]++;
    g_poolAnalytics.bucketTotalMisses[bucket]++;  // Track cumulative misses per bucket
    g_poolAnalytics.totalAllocTime += allocTime;
    g_poolAnalytics.totalAllocated += bucket;
    g_poolAnalytics.currentMemory += bucket;
    if (g_poolAnalytics.currentMemory > g_poolAnalytics.peakMemory) {
        g_poolAnalytics.peakMemory = g_poolAnalytics.currentMemory;
    }
    
    // Register pointer in pool tracking
    g_hostPointerSizes[newPtr] = bucket;
    
    CmiUnlock(g_hostPoolLock);
    *ptr = newPtr;
    return cudaSuccess;
}

/**
 * @brief Free host memory back to pool (raw, internal)
 * 
 * Returns the memory block to the pool for reuse. If the pointer
 * is unknown or the pool is not initialized, falls back to cudaFreeHost.
 * 
 * @param ptr Pointer to free
 * @return cudaError_t from operation
 */
cudaError_t hostPoolFreeRaw(void* ptr) {
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
    
    // Update analytics
    g_poolAnalytics.totalFrees++;
    g_poolAnalytics.timestepFrees++;
    g_poolAnalytics.currentMemory -= bucket;
    
    CmiUnlock(g_hostPoolLock);
    return cudaSuccess;
}

// ----------------------------------------------------------------------------
// Host Pool API with Logging - Pool-based Allocation
// ----------------------------------------------------------------------------

/**
 * @brief Host memory allocation from pool with logging
 * 
 * Allocates from the host pool and logs the event if logging is enabled.
 * Use the hostPoolMalloc macro instead of calling this directly.
 * 
 * @param ptr Output pointer to allocated memory
 * @param size Size to allocate in bytes
 * @param tag Pointer identifier tag
 * @param functionTag Function context tag
 * @param file Source file
 * @param line Line number
 * @return cudaError_t from allocation operation
 */
cudaError_t hostPoolMallocImpl(void** ptr, size_t size,
                               const char* tag, const char* functionTag,
                               const char* file, int line) {
    cudaError_t result = hostPoolMallocRaw(ptr, size);
    logHostAlloc(result, ptr, size, CkWallTimer(), file, line, tag, functionTag);
    return result;
}

/**
 * @brief Host memory deallocation to pool with logging
 * 
 * Frees to the pool and logs the event if logging is enabled.
 * Use the hostPoolFree macro instead of calling this directly.
 * 
 * @param ptr Pointer to free
 * @param tag Pointer identifier tag
 * @param functionTag Function context tag
 * @param file Source file
 * @param line Line number
 * @return cudaError_t from free operation
 */
cudaError_t hostPoolFreeImpl(void* ptr,
                             const char* tag, const char* functionTag,
                             const char* file, int line) {
    if (ptr == NULL) {
        logHostFree(cudaSuccess, ptr, CkWallTimer(), file, line, tag, functionTag);
        return cudaSuccess;
    }
    
    cudaError_t result = hostPoolFreeRaw(ptr);
    logHostFree(result, ptr, CkWallTimer(), file, line, tag, functionTag);
    return result;
}

// ----------------------------------------------------------------------------
// Host API with Logging - Direct Allocation (Bypass Pool)
// ----------------------------------------------------------------------------

/**
 * @brief Direct host memory allocation with logging (bypasses pool)
 * 
 * Allocates directly using cudaMallocHost and logs the event if logging is enabled.
 * Use the hostMalloc macro instead of calling this directly.
 * 
 * @param ptr Output pointer to allocated memory
 * @param size Size to allocate in bytes
 * @param tag Pointer identifier tag
 * @param functionTag Function context tag
 * @param file Source file
 * @param line Line number
 * @return cudaError_t from allocation operation
 */
cudaError_t hostMallocImpl(void** ptr, size_t size,
                           const char* tag, const char* functionTag,
                           const char* file, int line) {
    if (ptr == nullptr) {
        return cudaErrorInvalidValue;
    }
    
    if (size == 0) {
        *ptr = nullptr;
        return cudaSuccess;
    }
    
    // Use cudaMallocHost directly
    cudaError_t result = cudaMallocHost(ptr, size);
    logHostAlloc(result, ptr, size, CkWallTimer(), file, line, tag, functionTag);
    return result;
}

/**
 * @brief Direct host memory deallocation with logging (bypasses pool)
 * 
 * Frees directly using cudaFreeHost and logs the event if logging is enabled.
 * Use the hostFree macro instead of calling this directly.
 * 
 * @param ptr Pointer to free
 * @param tag Pointer identifier tag
 * @param functionTag Function context tag
 * @param file Source file
 * @param line Line number
 * @return cudaError_t from free operation
 */
cudaError_t hostFreeImpl(void* ptr,
                         const char* tag, const char* functionTag,
                         const char* file, int line) {
    if (ptr == NULL) {
        logHostFree(cudaSuccess, ptr, CkWallTimer(), file, line, tag, functionTag);
        return cudaSuccess;
    }
    
    // Direct free bypassing pool
    cudaError_t result = cudaFreeHost(ptr);
    logHostFree(result, ptr, CkWallTimer(), file, line, tag, functionTag);
    return result;
}

// ----------------------------------------------------------------------------
// Host Pool Maintenance & Trimming
// ----------------------------------------------------------------------------

/**
 * @brief Trim host pool to target capacity
 * 
 * Frees excess blocks from the pool to the OS, prioritizing large buckets
 * and respecting usage-based protection for recently-used buckets.
 * 
 * @param targetCapacityGB Target maximum pool capacity in GB
 * @param minCapacityPerBucketMB Minimum capacity to preserve per active bucket in MB
 * @return Total bytes freed to the OS
 */
size_t hostPoolTrim(double targetCapacityGB, size_t minCapacityPerBucketMB) {
    if (g_hostPoolLock == nullptr) {
        return 0;
    }
    
    size_t targetCapacity = (size_t)(targetCapacityGB * 1024.0 * 1024.0 * 1024.0);
    size_t minCapacityPerBucket = minCapacityPerBucketMB * 1024 * 1024;
    
    CmiLock(g_hostPoolLock);
    
    // Calculate current free pool capacity
    size_t currentCapacity = 0;
    for (const auto &kv : g_hostFreePool) {
        currentCapacity += kv.first * kv.second.size();
    }
    
    // If under target, nothing to trim
    if (currentCapacity <= targetCapacity) {
        CmiUnlock(g_hostPoolLock);
        return 0;
    }
    
    // Build list of (bucket, trimmableBlocks, bucketSize) sorted by bucketSize descending
    struct TrimCandidate {
        size_t bucket;
        int trimmableBlocks;
        size_t bytesPerBlock;
    };
    std::vector<TrimCandidate> candidates;
    
    for (auto &kv : g_hostFreePool) {
        size_t bucket = kv.first;
        int currentBlocks = kv.second.size();
        
        // Calculate protection for this bucket
        int protectedBlocks = 0;
        auto usageIt = g_bucketUsage.find(bucket);
        int recentHits = 0;
        if (usageIt != g_bucketUsage.end() && usageIt->second.recentHits > 0) {
            recentHits = usageIt->second.recentHits;
            // Protect based on recent usage, but limited by capacity
            int capacityLimit = (int)(minCapacityPerBucket / bucket);
            protectedBlocks = std::min(usageIt->second.recentHits, capacityLimit);
        }
        
        int trimmable = std::max(0, currentBlocks - protectedBlocks);
        
        if (trimmable > 0) {
            candidates.push_back({bucket, trimmable, bucket});
        }
    }
    
    // Sort by bucket size descending (trim largest first)
    std::sort(candidates.begin(), candidates.end(), 
              [](const TrimCandidate& a, const TrimCandidate& b) {
                  return a.bytesPerBlock > b.bytesPerBlock;
              });
    
    // Greedily trim blocks from largest buckets until target reached
    // Collect pointers to free while holding lock, then free them unlocked
    struct BlockToFree {
        void* ptr;
        size_t bucket;
    };
    std::vector<BlockToFree> blocksToFree;
    size_t totalFreed = 0;
    
    for (auto &candidate : candidates) {
        if (currentCapacity <= targetCapacity) {
            break;
        }
        
        auto &freeBlocks = g_hostFreePool[candidate.bucket];
        int blocksToFreeCount = std::min(candidate.trimmableBlocks, 
                                         (int)((currentCapacity - targetCapacity) / candidate.bucket) + 1);
        
        for (int i = 0; i < blocksToFreeCount && !freeBlocks.empty(); i++) {
            void* ptr = freeBlocks.back();
            freeBlocks.pop_back();
            g_hostPointerSizes.erase(ptr);
            blocksToFree.push_back({ptr, candidate.bucket});
            
            totalFreed += candidate.bucket;
            currentCapacity -= candidate.bucket;
        }
    }
    
    // Release lock before slow cudaFreeHost operations
    CmiUnlock(g_hostPoolLock);
    
    // Free all blocks unlocked (slow operation)
    for (const auto& block : blocksToFree) {
        cudaFreeHost(block.ptr);
    }
    
    // Re-acquire lock for analytics updates
    CmiLock(g_hostPoolLock);
    int blocksTrimmed = blocksToFree.size();
    for (const auto& block : blocksToFree) {
        g_poolAnalytics.totalFreedToOS++;
        g_poolAnalytics.timestepFreedToOS++;
        g_poolAnalytics.totalFreed += block.bucket;
    }
    CmiUnlock(g_hostPoolLock);
    
    return totalFreed;
}

/**
 * @brief Adaptive refill of hot buckets between timesteps
 * 
 * Refills buckets that had misses during the previous timestep to maintain
 * non-blocking behavior. Only refills if misses occurred and capacity is low.
 */
void hostPoolAdaptiveRefill() {
    if (g_hostPoolLock == nullptr) return;
    
    CmiLock(g_hostPoolLock);
    
    // Use shared warmup configuration for refill targets
    // Collect refill requirements while holding lock
    struct RefillRequest {
        size_t bucket;
        int count;
    };
    std::vector<RefillRequest> refillRequests;
    
    for (int i = 0; i < g_numWarmupSizes; i++) {
        size_t bucket = g_warmupConfig[i].size;
        auto missIt = g_bucketMissCount.find(bucket);
        
        // Only refill if we had misses and current free count is low
        if (missIt != g_bucketMissCount.end() && missIt->second > 0) {
            auto freeIt = g_hostFreePool.find(bucket);
            int currentFree = (freeIt != g_hostFreePool.end()) ? freeIt->second.size() : 0;
            
            // Refill up to 3 blocks if we had misses and have < 2 free
            int refillCount = std::min(3 - currentFree, missIt->second);
            if (refillCount > 0) {
                refillRequests.push_back({bucket, refillCount});
            }
        }
    }
    
    // Reset miss counts for next timestep
    g_bucketMissCount.clear();
    g_timestepCount++;
    
    // Decay recent hits to allow trimming of less-used buckets
    // Active buckets maintain protection through ongoing hits
    // Inactive buckets gradually lose protection and become trimmable
    for (auto &kv : g_bucketUsage) {
        if (kv.second.recentHits > 0) {
            kv.second.recentHits--;
        }
    }
    
    // Reset per-timestep analytics
    g_poolAnalytics.timestepAllocations = 0;
    g_poolAnalytics.timestepFrees = 0;
    g_poolAnalytics.timestepHits = 0;
    g_poolAnalytics.timestepMisses = 0;
    g_poolAnalytics.timestepFreedToOS = 0;
    
    CmiUnlock(g_hostPoolLock);
    
    // Allocate blocks unlocked (slow operation)
    struct AllocatedBlock {
        void* ptr;
        size_t bucket;
    };
    std::vector<AllocatedBlock> allocatedBlocks;
    int totalRefilled = 0;
    size_t totalRefillBytes = 0;
    
    for (const auto& req : refillRequests) {
        for (int j = 0; j < req.count; j++) {
            void* ptr = nullptr;
            cudaError_t res = cudaMallocHost(&ptr, req.bucket);
            if (res == cudaSuccess) {
                allocatedBlocks.push_back({ptr, req.bucket});
                totalRefilled++;
                totalRefillBytes += req.bucket;
            } else {
                break; // Stop if allocation fails
            }
        }
    }
    
    // Re-acquire lock to add allocated blocks to pool
    if (!allocatedBlocks.empty()) {
        CmiLock(g_hostPoolLock);
        for (const auto& block : allocatedBlocks) {
            g_hostFreePool[block.bucket].push_back(block.ptr);
            // Note: Free pool blocks are not tracked in g_hostPointerSizes
            // Only live (allocated) blocks are tracked there
        }
        CmiUnlock(g_hostPoolLock);
    }
    
    if (totalRefilled > 0) {
        CkPrintf("[HostPool] Refilled %d blocks (%.1f MB) for hot buckets\n", 
                 totalRefilled, totalRefillBytes / (1024.0 * 1024.0));
    }
}

// ----------------------------------------------------------------------------
// Host Pool Analytics & Reporting
// ----------------------------------------------------------------------------

/**
 * @brief Report host pool statistics with enhanced analytics
 * 
 * Prints current pool state including bucket count, free blocks, total memory,
 * hit rates, per-timestep statistics, and per-bucket distribution.
 * Also performs pool accounting validation and stores timestep snapshot.
 * 
 * @param prefix String prefix for output lines (for formatting)
 * @param targetCapacityGB Target pool capacity in GB (for trim reporting)
 */
void hostPoolReportStats(const char* prefix, double targetCapacityGB) {
    if (g_hostPoolLock == nullptr) {
        CkPrintf("%s Host pool not initialized\n", prefix);
        return;
    }
    
    CmiLock(g_hostPoolLock);
    
    size_t totalBlocks = 0;
    size_t totalBytes = 0;
    size_t numBuckets = g_hostFreePool.size();
    
    for (const auto &kv : g_hostFreePool) {
        size_t bucketSize = kv.first;
        size_t numBlocks = kv.second.size();
        totalBlocks += numBlocks;
        totalBytes += bucketSize * numBlocks;
    }
    
    size_t liveBlocks = g_hostPointerSizes.size();
    
    // Copy all analytics values while holding lock to avoid race conditions
    int totalAllocations = g_poolAnalytics.totalAllocations;
    int totalHits = g_poolAnalytics.totalHits;
    int totalMisses = g_poolAnalytics.totalMisses;
    int timestepAllocations = g_poolAnalytics.timestepAllocations;
    int timestepFrees = g_poolAnalytics.timestepFrees;
    int timestepHits = g_poolAnalytics.timestepHits;
    int timestepMisses = g_poolAnalytics.timestepMisses;
    int timestepFreedToOS = g_poolAnalytics.timestepFreedToOS;
    size_t peakMemory = g_poolAnalytics.peakMemory;
    size_t currentMemory = g_poolAnalytics.currentMemory;
    double totalAllocTime = g_poolAnalytics.totalAllocTime;
    size_t totalAllocated = g_poolAnalytics.totalAllocated;
    size_t totalFreed = g_poolAnalytics.totalFreed;
    size_t totalFreedToOS = g_poolAnalytics.totalFreedToOS;
    int timestepCount = g_timestepCount;
    
    CmiUnlock(g_hostPoolLock);
    
    // Calculate capacity usage
    double totalBytesGB = totalBytes / (1024.0 * 1024.0 * 1024.0);
    double capacityUsagePct = (targetCapacityGB > 0) ? (totalBytesGB / targetCapacityGB * 100.0) : 0.0;
    const char* capacityStatus = (totalBytesGB > targetCapacityGB) ? "OVER" : "OK";
    
    // Pool summary header
    CkPrintf("\n");
    CkPrintf("%s ═══════════════════════════════════════════════════════════════\n", prefix);
    CkPrintf("%s HOST POOL STATUS\n", prefix);
    CkPrintf("%s ═══════════════════════════════════════════════════════════════\n", prefix);
    
    // Capacity metrics
    CkPrintf("%s Free Pool Capacity:  %.2f GB / %.2f GB target (%.1f%% %s)\n",
             prefix, totalBytesGB, targetCapacityGB, capacityUsagePct, capacityStatus);
    CkPrintf("%s Pool Structure:      %zu buckets, %zu free blocks, %zu live blocks\n",
             prefix, numBuckets, totalBlocks, liveBlocks);
    
    // Hit rate and efficiency
    double hitRate = (totalAllocations > 0) ? 
        (100.0 * totalHits) / totalAllocations : 0.0;
    
    CkPrintf("%s Hit Rate:            %.1f%% (%d hits, %d misses)\n",
             prefix, hitRate, totalHits, totalMisses);
    
    // Memory usage
    CkPrintf("%s Memory Usage:        Peak %.2f GB, Current %.2f GB\n",
             prefix, 
             peakMemory / (1024.0 * 1024.0 * 1024.0),
             currentMemory / (1024.0 * 1024.0 * 1024.0));
    
    // This timestep activity
    CkPrintf("%s This Step:           %d allocs, %d frees, %d hits, %d misses, %d freed to OS\n",
             prefix, timestepAllocations, timestepFrees,
             timestepHits, timestepMisses, timestepFreedToOS);
    
    // Cumulative stats
    CkPrintf("%s Cumulative:          %d allocs, %.3f ms avg time, %.2f GB from OS\n",
             prefix, totalAllocations,
             (totalAllocTime * 1000.0) / std::max(1, totalAllocations),
             totalAllocated / (1024.0 * 1024.0 * 1024.0));
    
    // Pool accounting check
    // Total blocks owned by pool = (allocated from OS) - (freed to OS)
    // These blocks are either: (1) free in pool, or (2) currently checked out (live)
    size_t totalOwnedBlocks = totalMisses - totalFreedToOS;
    size_t accountedBlocks = totalBlocks + liveBlocks;
    
    if (totalOwnedBlocks != accountedBlocks) {
        size_t diff = (totalOwnedBlocks > accountedBlocks) ? 
                      (totalOwnedBlocks - accountedBlocks) : (accountedBlocks - totalOwnedBlocks);
        CkPrintf("%s WARNING: Pool accounting mismatch! Owned: %zu, Accounted: %zu (Free: %zu + Live: %zu), Diff: %zu\n",
                 prefix, totalOwnedBlocks, accountedBlocks, totalBlocks, liveBlocks, diff);
    }
    
    // Per-bucket diagnostics (sorted by bucket size)
    CkPrintf("%s Per-Bucket Distribution:\n", prefix);
    
    CmiLock(g_hostPoolLock);
    
    // Collect and sort buckets by size
    std::vector<size_t> bucketSizes;
    for (const auto &kv : g_hostFreePool) {
        bucketSizes.push_back(kv.first);
    }
    std::sort(bucketSizes.begin(), bucketSizes.end());
    
    // Print header
    CkPrintf("%s   %10s | %6s | %6s | %10s | %8s\n", prefix, 
             "Bucket", "Free", "Misses", "TotalAllocs", "Usage");
    CkPrintf("%s   %s\n", prefix, "----------------------------------------------------------------");
    
    for (size_t bucket : bucketSizes) {
        int freeCount = g_hostFreePool[bucket].size();
        int totalMisses = g_poolAnalytics.bucketTotalMisses.count(bucket) ? 
                         g_poolAnalytics.bucketTotalMisses[bucket] : 0;
        int totalAllocs = g_poolAnalytics.bucketTotalAllocations.count(bucket) ? 
                         g_poolAnalytics.bucketTotalAllocations[bucket] : 0;
        
        // Format bucket size nicely
        char bucketStr[16];
        if (bucket < 1024) {
            snprintf(bucketStr, sizeof(bucketStr), "%zu B", bucket);
        } else if (bucket < 1024 * 1024) {
            snprintf(bucketStr, sizeof(bucketStr), "%zu KB", bucket / 1024);
        } else {
            snprintf(bucketStr, sizeof(bucketStr), "%zu MB", bucket / (1024 * 1024));
        }
        
        // Calculate usage percentage: (blocks in use) / (total blocks allocated from OS) * 100
        // Total blocks from OS = totalMisses (each miss allocates a new block from OS)
        // Blocks in use = totalMisses - freeCount
        double usagePct = (totalMisses > 0) ? (100.0 * (totalMisses - freeCount) / totalMisses) : 0.0;
        
        CkPrintf("%s   %10s | %6d | %6d | %10d | %7.1f%%\n", prefix,
                 bucketStr, freeCount, totalMisses, totalAllocs, usagePct);
    }
    
    CkPrintf("%s   %s\n", prefix, "----------------------------------------------------------------");
    CkPrintf("%s ═══════════════════════════════════════════════════════════════\n", prefix);
    CkPrintf("\n");
    
    CmiUnlock(g_hostPoolLock);  // Release lock after per-bucket diagnostics
    
    // Store timestep snapshot for historical tracking (protected by lock)
    TimestepSnapshot snapshot;
    snapshot.timestep = timestepCount;
    snapshot.totalBuckets = numBuckets;
    snapshot.totalFreeBlocks = totalBlocks;
    snapshot.totalFreeBytes = totalBytes;
    snapshot.liveBlocks = liveBlocks;
    snapshot.totalAllocated = totalAllocated;
    snapshot.totalFreed = totalFreed;
    snapshot.timestamp = CkWallTimer();
    
    // Re-acquire lock to update history
    CmiLock(g_hostPoolLock);
    g_timestepHistory.push_back(snapshot);
    
    // Keep only last 100 timesteps to prevent memory growth
    if (g_timestepHistory.size() > 100) {
        g_timestepHistory.erase(g_timestepHistory.begin());
    }
    CmiUnlock(g_hostPoolLock);
}

/**
 * @brief Analyze warmup effectiveness by comparing targets vs. actual usage
 * 
 * Identifies over-provisioned and under-provisioned buckets based on
 * warmup targets vs. actual free counts and misses. Provides recommendations
 * for adjusting warmup parameters.
 */
void hostPoolAnalyzeWarmupEffectiveness() {
    if (g_hostPoolLock == nullptr) {
        CkPrintf("[PoolAnalysis] Host pool not initialized\n");
        return;
    }
    
    // Build warmup targets from shared configuration - always in sync
    std::unordered_map<size_t, int> warmupTargets;
    for (int i = 0; i < g_numWarmupSizes; i++) {
        warmupTargets[g_warmupConfig[i].size] = g_warmupConfig[i].count;
    }
    
    CkPrintf("\n=== WARMUP EFFECTIVENESS ANALYSIS ===\n");
    
    CmiLock(g_hostPoolLock);
    
    // Collect buckets sorted by size
    std::vector<size_t> bucketSizes;
    for (const auto &kv : g_hostFreePool) {
        bucketSizes.push_back(kv.first);
    }
    std::sort(bucketSizes.begin(), bucketSizes.end());
    
    CkPrintf("Bucket Analysis (Target vs. Actual):\n");
    CkPrintf("  %10s | %8s | %6s | %6s | %10s | Status\n", 
             "Bucket", "Target", "Free", "Misses", "TotalAllocs", "");
    CkPrintf("  --------------------------------------------------------------------\n");
    
    int overProvisionedCount = 0;
    int underProvisionedCount = 0;
    int wellProvisionedCount = 0;
    
    for (size_t bucket : bucketSizes) {
        int target = warmupTargets.count(bucket) ? warmupTargets[bucket] : 0;
        int freeCount = g_hostFreePool[bucket].size();
        int totalMisses = g_poolAnalytics.bucketTotalMisses.count(bucket) ? 
                         g_poolAnalytics.bucketTotalMisses[bucket] : 0;
        int totalAllocs = g_poolAnalytics.bucketTotalAllocations.count(bucket) ? 
                         g_poolAnalytics.bucketTotalAllocations[bucket] : 0;
        
        // Format bucket size
        char bucketStr[16];
        if (bucket < 1024) {
            snprintf(bucketStr, sizeof(bucketStr), "%zu B", bucket);
        } else if (bucket < 1024 * 1024) {
            snprintf(bucketStr, sizeof(bucketStr), "%zu KB", bucket / 1024);
        } else {
            snprintf(bucketStr, sizeof(bucketStr), "%zu MB", bucket / (1024 * 1024));
        }
        
        // Determine status
        const char* status = "";
        if (target == 0) {
            // Bucket not in warmup - likely dynamic allocation
            status = "DYNAMIC";
        } else if (freeCount > target * 1.5) {
            // Over-provisioned: way more free than target
            status = "OVER";
            overProvisionedCount++;
        } else if (totalMisses > target * 0.1) {
            // Under-provisioned: significant misses after warmup
            status = "UNDER";
            underProvisionedCount++;
        } else {
            // Well-provisioned: close to target, few misses
            status = "OK";
            wellProvisionedCount++;
        }
        
        CkPrintf("  %10s | %8d | %6d | %6d | %10d | %s\n",
                 bucketStr, target, freeCount, totalMisses, totalAllocs, status);
    }
    
    CmiUnlock(g_hostPoolLock);
    
    // Summary
    CkPrintf("\nSummary:\n");
    CkPrintf("  Well-provisioned buckets: %d\n", wellProvisionedCount);
    CkPrintf("  Over-provisioned buckets: %d (reduce warmup counts)\n", overProvisionedCount);
    CkPrintf("  Under-provisioned buckets: %d (increase warmup counts)\n", underProvisionedCount);
    
    // Recommendations
    if (overProvisionedCount > 0) {
        CkPrintf("\nRecommendation: Reduce warmup counts for OVER buckets to save initial memory.\n");
    }
    if (underProvisionedCount > 0) {
        CkPrintf("\nRecommendation: Increase warmup counts for UNDER buckets to reduce misses.\n");
    }
    if (wellProvisionedCount == bucketSizes.size()) {
        CkPrintf("\nExcellent! All buckets are well-provisioned.\n");
    }
    
    CkPrintf("=== END WARMUP ANALYSIS ===\n\n");
}

/**
 * @brief Analyze memory growth trends over time
 * 
 * Reports memory growth patterns and detects potential leaks by analyzing
 * historical timestep snapshots. Identifies concerning trends in bucket
 * proliferation, memory growth, and live block accumulation.
 */
void hostPoolAnalyzeGrowth() {
    if (g_hostPoolLock == nullptr) {
        CkPrintf("[PoolAnalysis] Host pool not initialized\n");
        return;
    }
    
    // Protect all access to g_timestepHistory with lock
    CmiLock(g_hostPoolLock);
    
    if (g_timestepHistory.size() < 2) {
        CmiUnlock(g_hostPoolLock);
        CkPrintf("[PoolAnalysis] Not enough data for growth analysis (need at least 2 timesteps)\n");
        return;
    }
    
    // Copy data we need while holding lock
    size_t historySize = g_timestepHistory.size();
    TimestepSnapshot first = g_timestepHistory.front();
    TimestepSnapshot last = g_timestepHistory.back();
    
    // Copy recent snapshots if available
    std::vector<TimestepSnapshot> recentSnapshots;
    if (historySize >= 5) {
        for (int i = historySize - 5; i < (int)historySize; i++) {
            recentSnapshots.push_back(g_timestepHistory[i]);
        }
    }
    
    CmiUnlock(g_hostPoolLock);
    
    // Now do calculations and printing unlocked (using copied data)
    CkPrintf("\n=== POOL GROWTH ANALYSIS ===\n");
    
    // Calculate growth rates
    double timeSpan = last.timestamp - first.timestamp;
    size_t bucketGrowth = last.totalBuckets - first.totalBuckets;
    size_t memoryGrowth = last.totalFreeBytes - first.totalFreeBytes;
    size_t liveGrowth = last.liveBlocks - first.liveBlocks;
    
    CkPrintf("Time span: %.1f seconds (%d timesteps)\n", timeSpan, last.timestep - first.timestep + 1);
    CkPrintf("Bucket growth: %zu buckets (%.1f per timestep)\n", 
             bucketGrowth, (double)bucketGrowth / (last.timestep - first.timestep + 1));
    CkPrintf("Memory growth: %.2f GB (%.3f GB per timestep)\n", 
             memoryGrowth / (1024.0 * 1024.0 * 1024.0),
             memoryGrowth / (1024.0 * 1024.0 * 1024.0) / (last.timestep - first.timestep + 1));
    CkPrintf("Live block growth: %zu blocks (%.1f per timestep)\n",
             liveGrowth, (double)liveGrowth / (last.timestep - first.timestep + 1));
    
    // Detect concerning patterns
    if (bucketGrowth > 10) {
        CkPrintf("WARNING: High bucket proliferation detected! Consider adjusting bucket granularity.\n");
    }
    
    if (memoryGrowth > 1024 * 1024 * 1024) { // > 1GB
        CkPrintf("WARNING: Significant memory growth detected! Check for memory leaks.\n");
    }
    
    if (liveGrowth > 0) {
        CkPrintf("WARNING: Live block count increasing - potential memory leak!\n");
    }
    
    // Show recent trend (last 5 timesteps)
    if (!recentSnapshots.empty()) {
        CkPrintf("\nRecent trend (last 5 timesteps):\n");
        for (const auto& snap : recentSnapshots) {
            CkPrintf("  Step %d: %zu buckets, %.2f GB free, %zu live\n",
                     snap.timestep, snap.totalBuckets, 
                     snap.totalFreeBytes / (1024.0 * 1024.0 * 1024.0), snap.liveBlocks);
        }
    }
    
    CkPrintf("=== END ANALYSIS ===\n\n");
}

// ============================================================================
// MEMORY LOGGING INFRASTRUCTURE
// ============================================================================

// ----------------------------------------------------------------------------
// DataManager Log Initialization
// ----------------------------------------------------------------------------

/**
 * @brief Initialize GPU memory log on DataManager
 * 
 * Called collectively across all DataManagers to set the log filename
 * and enable/disable GPU memory logging.
 * 
 * @param _fileName Output filename for GPU memory log
 * @param bGpuMemLoggerFlag Flag to enable (1) or disable (0) GPU logging
 * @param cb Callback to contribute when initialization completes
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
 * 
 * Called collectively across all DataManagers to set the log filename
 * and enable/disable CPU memory logging.
 * 
 * @param _fileName Output filename for CPU memory log
 * @param bCpuMemLoggerFlag Flag to enable (1) or disable (0) CPU logging
 * @param cb Callback to contribute when initialization completes
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

// ----------------------------------------------------------------------------
// Main Chare Log Initialization
// ----------------------------------------------------------------------------

/**
 * @brief Initialize GPU memory log file on PE 0
 * 
 * Called from Main chare to create the GPU memory log file and write header.
 * Sends configuration to all DataManagers before creating the file.
 */
void Main::initMemLog() {
    std::string memLogFile = "gpumemlog.out";
    
    // Send filename and flag to all DataManagers and wait
    dMProxy.initMemLog(memLogFile, param.bGpuMemLogger, CkCallbackResumeThread());
    
    // PE 0 creates the file and writes header
    if (CkMyPe() == 0) {
        FILE* fpLog = CmiFopen(memLogFile.c_str(), "w");
        fprintf(fpLog, "# ChaNGa Memory Log v1.1\n");
        fprintf(fpLog, "# NodeID OpType Size Address Timestamp File:Line PointerID FunctionTag\n");
        int close_err = CmiFclose(fpLog);
        if (close_err != 0) {
            CkPrintf("WARNING: PE 0 failed to close gpumemlog file: %s (Error %d)\n", memLogFile.c_str(), close_err);
        }
    }
}

/**
 * @brief Initialize CPU memory log file on PE 0
 * 
 * Called from Main chare to create the CPU memory log file and write header.
 * Sends configuration to all DataManagers before creating the file.
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

// ----------------------------------------------------------------------------
// DataManager Log Flushing
// ----------------------------------------------------------------------------

/**
 * @brief Flush GPU memlog sequentially across nodes
 * 
 * Sequentially flushes GPU memory logs from each node to avoid file corruption.
 * Each node triggers the next node in sequence before invoking the callback.
 * 
 * @param cb Callback to invoke after all nodes have flushed
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
 * 
 * Sequentially flushes CPU memory logs from each node to avoid file corruption.
 * Each node triggers the next node in sequence before invoking the callback.
 * 
 * @param cb Callback to invoke after all nodes have flushed
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

// ----------------------------------------------------------------------------
// MemLog Implementation
// ----------------------------------------------------------------------------

/**
 * @brief Write buffered memory log events to file
 * 
 * Flushes all buffered memory log events to the log file in append mode.
 * Converts operation types to human-readable strings and formats each event.
 * Clears the event buffer after writing.
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

