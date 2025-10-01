#ifndef MEMORY_POOL_H
#define MEMORY_POOL_H

#ifdef CUDA
#include <cuda_runtime.h>

/**
 * @brief Initialize the GPU memory pool
 * Configures the default CUDA memory pool to retain freed memory
 */
void poolInit();

/**
 * @brief Initialize GPU pool with automatic warmup for typical ChaNGa workloads
 * Combines poolInit() with pre-warming based on profiled allocation patterns
 */
void poolInitWithWarmup();

/**
 * @brief Pre-warm the GPU memory pool with common allocation sizes
 * Allocates and immediately frees blocks to populate the pool's free-list
 * @param sizes Array of allocation sizes in bytes
 * @param counts Array of how many allocations per size
 * @param numSizes Number of different sizes to pre-allocate
 * @param stream CUDA stream to use for allocations (use 0 for default)
 */
cudaError_t poolWarmup(const size_t* sizes, const int* counts, int numSizes, cudaStream_t stream);

/**
 * @brief Allocate GPU memory from the pool (stream-aware)
 * @param devPtr Output pointer to the allocated device memory
 * @param size The size of memory to allocate in bytes
 * @param stream The CUDA stream that will use this memory
 * @return cudaError_t result from the allocation operation
 */
cudaError_t poolMalloc(void** devPtr, size_t size, cudaStream_t stream);

/**
 * @brief Free GPU memory back to the pool (stream-aware)
 * @param ptr The device pointer to return to the pool
 * @param stream The CUDA stream on which the memory was last used
 * @return cudaError_t result, typically cudaSuccess
 */
cudaError_t poolFree(void* ptr, cudaStream_t stream);

/**
 * @brief Destroy the memory pool and cleanup host allocations
 */
void poolDestroy();

// ============================================================================
// Host Memory Pool (Pinned Memory)
// ============================================================================

/**
 * @brief Initialize host memory pool (no-op, lazy initialization)
 */
cudaError_t hostMemoryPoolInit(size_t initialSize);

/**
 * @brief Warm up host pool with typical ChaNGa allocation patterns
 * Pre-allocates and frees common sizes to populate the pool
 */
void hostPoolInitWithWarmup();

/**
 * @brief Allocate pinned host memory from the pool
 * Uses bucketed free-list for reuse (64KB granularity)
 */
cudaError_t hostMallocPool(void** ptr, size_t size);

/**
 * @brief Stream-aware host allocation (currently ignores stream)
 */
cudaError_t hostMallocPoolStream(void** ptr, size_t size, cudaStream_t stream);

/**
 * @brief Free pinned host memory back to the pool
 */
cudaError_t hostFreePool(void* ptr);

/**
 * @brief Stream-aware host free (currently ignores stream)
 */
cudaError_t hostFreePoolStream(void* ptr, cudaStream_t stream);

/**
 * @brief Cleanup all host pool allocations
 */
cudaError_t hostMemoryPoolCleanup();

// ============================================================================
// GPU Memory Allocation with Logging
// ============================================================================

/**
 * @brief GPU memory allocation with optional logging
 */
cudaError_t gpuMalloc(void** devPtr, size_t size, cudaStream_t stream, 
                      const char* tag, const char* functionTag, 
                      const char* file, int line);

/**
 * @brief GPU memory deallocation with optional logging
 */
cudaError_t gpuFree(void* devPtr, cudaStream_t stream, 
                    const char* tag, const char* functionTag, 
                    const char* file, int line);

// Templated wrapper for type safety
template <typename T>
inline cudaError_t gpuMallocTyped(T** ptr, size_t size, cudaStream_t stream, 
                                  const char* pointerIdTag, const char* functionTag, 
                                  const char* file, int line) {
    return gpuMalloc(reinterpret_cast<void**>(ptr), size, stream, 
                     pointerIdTag, functionTag, file, line);
}

// Convenience macros
#define gpuMallocHelper(ptr, size, stream, funcTag) \
    gpuMallocTyped(ptr, size, stream, #ptr, funcTag, __FILE__, __LINE__)

#define gpuFreeHelper(ptr, stream, funcTag) \
    gpuFree(ptr, stream, #ptr, funcTag, __FILE__, __LINE__)

#endif // CUDA

#endif // MEMORY_POOL_H

