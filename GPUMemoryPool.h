#ifndef GPU_MEMORY_POOL_H
#define GPU_MEMORY_POOL_H

#include <cuda_runtime.h>
#include <stddef.h> // For size_t

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Wrapper for cudaMalloc that prints allocation size and context.
 * @param devPtr Pointer to the allocated device memory.
 * @param size Size of memory to allocate in bytes.
 * @param tag String representing the variable name being allocated.
 * @param functionTag String representing the calling function/context.
 * @param file Source file where the allocation occurs.
 * @param line Line number where the allocation occurs.
 * @return cudaError_t result from cudaMalloc.
 */
cudaError_t gpuMallocTracked(void** devPtr, size_t size, const char* tag, const char* functionTag, const char* file, int line);

/**
 * @brief Wrapper for cudaFree that prints a message and context.
 * @param devPtr Pointer to the device memory to free.
 * @param tag String representing the variable name being freed.
 * @param functionTag String representing the calling function/context.
 * @param file Source file where the free occurs.
 * @param line Line number where the free occurs.
 * @return cudaError_t result from cudaFree.
 */
cudaError_t gpuFreeTracked(void* devPtr, const char* tag, const char* functionTag, const char* file, int line);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
// Templated inline wrapper to handle type casting safely and pass context
template <typename T>
inline cudaError_t gpuMallocTyped(T** ptr, size_t size, const char* pointerIdTag, const char* functionTag, const char* file, int line) {
    // Call the C-linkage function
    return gpuMallocTracked(reinterpret_cast<void**>(ptr), size, pointerIdTag, functionTag, file, line);
}

// Macro helper to call the typed wrapper, capturing context automatically
#define gpuMallocHelper(ptr, size, funcTag) gpuMallocTyped(ptr, size, #ptr, funcTag, __FILE__, __LINE__)

// Macro helper for gpuFree, capturing context automatically
#define gpuFreeHelper(ptr, funcTag) gpuFreeTracked(ptr, #ptr, funcTag, __FILE__, __LINE__)
#endif // __cplusplus

#endif // GPU_MEMORY_POOL_H 