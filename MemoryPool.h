#ifndef MEMORY_POOL_H
#define MEMORY_POOL_H

#ifdef CUDA
#include <cuda_runtime.h>

// ============================================================================
// GPU MEMORY POOL
// ============================================================================

// ----------------------------------------------------------------------------
// GPU Pool Initialization & Warmup
// ----------------------------------------------------------------------------

/**
 * @brief Initialize the GPU memory pool
 * 
 * Configures the default CUDA memory pool with a 5GB retention limit
 * and pre-warms the pool with common allocation sizes.
 */
void gpuPoolInit();

/**
 * @brief Pre-warm the GPU memory pool with common allocation sizes
 * 
 * Allocates and immediately frees blocks to populate the pool's free-list.
 * 
 * @param sizes Array of allocation sizes in bytes
 * @param counts Array of how many allocations per size
 * @param numSizes Number of different sizes to pre-allocate
 * @param stream CUDA stream to use for allocations (use 0 for default)
 * @return cudaError_t from allocation operations
 */
cudaError_t gpuPoolWarmup(const size_t* sizes, const int* counts, int numSizes, cudaStream_t stream);

// ----------------------------------------------------------------------------
// GPU Pool Raw Operations (Internal - use macros for public API)
// ----------------------------------------------------------------------------

/**
 * @brief Allocate GPU memory from the pool (raw, internal, stream-aware)
 * 
 * @param devPtr Output pointer to allocated device memory
 * @param size Size of memory to allocate in bytes
 * @param stream CUDA stream for ordering (NULL for synchronous)
 * @return cudaError_t result from allocation
 */
cudaError_t gpuPoolMallocRaw(void** devPtr, size_t size, cudaStream_t stream);

/**
 * @brief Free GPU memory back to the pool (raw, internal, stream-aware)
 * 
 * @param ptr Device pointer to free
 * @param stream CUDA stream for ordering (NULL for synchronous)
 * @return cudaError_t result from free operation
 */
cudaError_t gpuPoolFreeRaw(void* ptr, cudaStream_t stream);

// ----------------------------------------------------------------------------
// GPU Pool API with Logging - Pool-based Allocation
// ----------------------------------------------------------------------------

/**
 * @brief GPU memory allocation from pool with logging (use gpuPoolMalloc macro)
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
                              const char* file, int line);

/**
 * @brief GPU memory deallocation to pool with logging (use gpuPoolFree macro)
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
                            const char* file, int line);

// ----------------------------------------------------------------------------
// GPU API with Logging - Direct Allocation (Bypass Pool)
// ----------------------------------------------------------------------------

/**
 * @brief Direct GPU memory allocation with logging (use gpuMalloc macro)
 * 
 * Bypasses the pool and allocates directly using cudaMalloc.
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
                          const char* file, int line);

/**
 * @brief Direct GPU memory deallocation with logging (use gpuFree macro)
 * 
 * Bypasses the pool and frees directly using cudaFree.
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
                        const char* file, int line);

// ----------------------------------------------------------------------------
// GPU Pool Type-Safe Wrappers (Templates)
// ----------------------------------------------------------------------------

/**
 * @brief Type-safe wrapper for GPU pool allocation
 */
template <typename T>
inline cudaError_t gpuPoolMallocTyped(T** ptr, size_t size, cudaStream_t stream, 
                                      const char* pointerIdTag, const char* functionTag, 
                                      const char* file, int line) {
    return gpuPoolMallocImpl(reinterpret_cast<void**>(ptr), size, stream, 
                             pointerIdTag, functionTag, file, line);
}

/**
 * @brief Type-safe wrapper for direct GPU allocation (bypass pool)
 */
template <typename T>
inline cudaError_t gpuMallocTyped(T** ptr, size_t size, cudaStream_t stream, 
                                  const char* pointerIdTag, const char* functionTag, 
                                  const char* file, int line) {
    return gpuMallocImpl(reinterpret_cast<void**>(ptr), size, stream, 
                         pointerIdTag, functionTag, file, line);
}

// ----------------------------------------------------------------------------
// GPU Pool Public API Macros
// ----------------------------------------------------------------------------

// Pool-based allocation
#define gpuPoolMalloc(ptr, size, stream, funcTag) \
    gpuPoolMallocTyped(ptr, size, stream, #ptr, funcTag, __FILE__, __LINE__)

#define gpuPoolFree(ptr, stream, funcTag) \
    gpuPoolFreeImpl(ptr, stream, #ptr, funcTag, __FILE__, __LINE__)

// Direct allocation (bypasses pool)
#define gpuMalloc(ptr, size, stream, funcTag) \
    gpuMallocTyped(ptr, size, stream, #ptr, funcTag, __FILE__, __LINE__)

#define gpuFree(ptr, stream, funcTag) \
    gpuFreeImpl(ptr, stream, #ptr, funcTag, __FILE__, __LINE__)

// ----------------------------------------------------------------------------
// GPU Pool Cleanup
// ----------------------------------------------------------------------------

/**
 * @brief Destroy the GPU memory pool and cleanup
 */
void gpuPoolDestroy();

// ============================================================================
// Host Memory Pool (Pinned Memory)
// ============================================================================

// ----------------------------------------------------------------------------
// Host Pool Initialization & Warmup
// ----------------------------------------------------------------------------

/**
 * @brief Initialize the host memory pool
 * 
 * Creates the pool lock and pre-warms the pool with common allocation sizes
 * ranging from 256B to 64MB across 19 bucket sizes.
 */
void hostPoolInit();

// ----------------------------------------------------------------------------
// Host Pool Raw Operations (Internal - use macros for public API)
// ----------------------------------------------------------------------------

/**
 * @brief Allocate host memory from pool (raw, internal)
 * 
 * @param ptr Output pointer to allocated memory
 * @param size Size to allocate in bytes
 * @return cudaError_t from allocation operation
 */
cudaError_t hostPoolMallocRaw(void** ptr, size_t size);

/**
 * @brief Free host memory back to pool (raw, internal)
 * 
 * @param ptr Pointer to free
 * @return cudaError_t from operation
 */
cudaError_t hostPoolFreeRaw(void* ptr);

// ----------------------------------------------------------------------------
// Host Pool Maintenance & Trimming
// ----------------------------------------------------------------------------

/**
 * @brief Get warmup target for a given bucket size
 * 
 * @param bucket Bucket size in bytes
 * @return Target number of blocks (0 = no target, use dynamic trimming)
 */
int hostPoolGetWarmupTarget(size_t bucket);

/**
 * @brief Trim host pool to target capacity
 * 
 * Frees excess blocks from the pool to the OS, prioritizing large buckets.
 * 
 * @param targetCapacityGB Target maximum pool capacity in GB
 * @param minCapacityPerBucketMB Minimum capacity to preserve per active bucket in MB
 * @return Total bytes freed to the OS
 */
size_t hostPoolTrim(double targetCapacityGB, size_t minCapacityPerBucketMB);

/**
 * @brief Adaptive refill of hot buckets between timesteps
 * 
 * Refills buckets that had misses during the previous timestep.
 */
void hostPoolAdaptiveRefill();

// ----------------------------------------------------------------------------
// Host Pool Analytics & Reporting
// ----------------------------------------------------------------------------

/**
 * @brief Report host pool statistics with enhanced analytics
 * 
 * @param prefix String prefix for output lines (for formatting)
 * @param targetCapacityGB Target pool capacity in GB (for trim reporting)
 */
void hostPoolReportStats(const char* prefix, double targetCapacityGB);

/**
 * @brief Analyze warmup effectiveness by comparing targets vs. actual usage
 */
void hostPoolAnalyzeWarmupEffectiveness();

/**
 * @brief Analyze memory growth trends over time
 */
void hostPoolAnalyzeGrowth();

// ----------------------------------------------------------------------------
// Host Pool API with Logging - Pool-based Allocation
// ----------------------------------------------------------------------------

/**
 * @brief Host memory allocation from pool with logging (use hostPoolMalloc macro)
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
                               const char* file, int line);

/**
 * @brief Host memory deallocation to pool with logging (use hostPoolFree macro)
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
                             const char* file, int line);

// ----------------------------------------------------------------------------
// Host API with Logging - Direct Allocation (Bypass Pool)
// ----------------------------------------------------------------------------

/**
 * @brief Direct host memory allocation with logging (use hostMalloc macro)
 * 
 * Bypasses the pool and allocates directly using cudaMallocHost.
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
                           const char* file, int line);

/**
 * @brief Direct host memory deallocation with logging (use hostFree macro)
 * 
 * Bypasses the pool and frees directly using cudaFreeHost.
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
                         const char* file, int line);

// ----------------------------------------------------------------------------
// Host Pool Type-Safe Wrappers (Templates)
// ----------------------------------------------------------------------------

/**
 * @brief Type-safe wrapper for host pool allocation
 */
template <typename T>
inline cudaError_t hostPoolMallocTyped(T** ptr, size_t size,
                                       const char* pointerIdTag, const char* functionTag,
                                       const char* file, int line) {
    return hostPoolMallocImpl(reinterpret_cast<void**>(ptr), size,
                              pointerIdTag, functionTag, file, line);
}

/**
 * @brief Type-safe wrapper for direct host allocation (bypass pool)
 */
template <typename T>
inline cudaError_t hostMallocTyped(T** ptr, size_t size,
                                   const char* pointerIdTag, const char* functionTag,
                                   const char* file, int line) {
    return hostMallocImpl(reinterpret_cast<void**>(ptr), size,
                          pointerIdTag, functionTag, file, line);
}

// ----------------------------------------------------------------------------
// Host Pool Public API Macros
// ----------------------------------------------------------------------------

// Pool-based allocation
#define hostPoolMalloc(ptr, size, funcTag) \
    hostPoolMallocTyped(ptr, size, #ptr, funcTag, __FILE__, __LINE__)

#define hostPoolFree(ptr, funcTag) \
    hostPoolFreeImpl(ptr, #ptr, funcTag, __FILE__, __LINE__)

// Direct allocation (bypasses pool)
#define hostMalloc(ptr, size, funcTag) \
    hostMallocTyped(ptr, size, #ptr, funcTag, __FILE__, __LINE__)

#define hostFree(ptr, funcTag) \
    hostFreeImpl(ptr, #ptr, funcTag, __FILE__, __LINE__)

// ----------------------------------------------------------------------------
// Host Pool Cleanup
// ----------------------------------------------------------------------------

/**
 * @brief Cleanup host memory pool
 * 
 * Frees all blocks from the pool back to the OS and clears internal data structures.
 * 
 * @return cudaSuccess after cleanup completes
 */
cudaError_t hostPoolCleanup();

#endif // CUDA

#endif // MEMORY_POOL_H
