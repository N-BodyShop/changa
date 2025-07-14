
cudaError_t gpuMallocTracked(void** devPtr, size_t size, const char* tag, const char* functionTag, const char* file, int line);
cudaError_t gpuFreeTracked(void* devPtr, const char* tag, const char* functionTag, const char* file, int line);


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
