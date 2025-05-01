#ifndef MEMLOG_HINCLUDED
#define MEMLOG_HINCLUDED

#include <vector>
#include <string>
#include <cstdint> // For uintptr_t
#include "pup_stl.h" // Include if PUP support is needed later, for now keep std::vector/string

// Enum to represent memory operation types
enum MemLogOpType {
    MEMLOG_ALLOC,
    MEMLOG_FREE,
    MEMLOG_ALLOC_FAIL,
    MEMLOG_FREE_FAIL,
    MEMLOG_FREE_SKIP // For attempting to free NULL
};

/** @brief Holds information about a single CUDA memory event */
class MemLogEvent {
public:
    int nodeId;               // ID of the node where the event occurred (CkMyNode())
    MemLogOpType opType;      // Type of memory operation
    size_t size;              // Size of memory block (for alloc/alloc_fail)
    uintptr_t address;        // Memory address (device pointer)
    double timestamp;         // Timestamp of the event (e.g., CkWallTimer())
    std::string location;     // Code location ("file:line")
    std::string tag;          // User-provided tag for the operation

    // Default constructor
    MemLogEvent() : nodeId(-1), opType(MEMLOG_ALLOC), size(0), address(0), timestamp(0.0) {}

    // Constructor
    MemLogEvent(int _nodeId, MemLogOpType _opType, size_t _size, uintptr_t _address, double _timestamp, const char* _file, int _line, const char* _tag)
        : nodeId(_nodeId), opType(_opType), size(_size), address(_address), timestamp(_timestamp), tag(_tag ? _tag : "") {
        // Combine file and line into location string
        location = std::string(_file ? _file : "unknown") + ":" + std::to_string(_line);
    }

    // PUP routine (optional, uncomment if migration needed)
    /*
    void pup(PUP::er& p) {
        p | opType; // Need to handle enum pup
        p | size;
        p | address;
        p | timestamp;
        p | location;
        p | tag;
    }
    */
};

/** @brief Manages buffering and flushing of memory log events */
class MemLog {
public:
    std::vector<MemLogEvent> meTab; // Buffer for memory events
    std::string fileName;          // Output log file name

    MemLog() : fileName(".memlog") {} // Default filename

    // Method to flush buffered events to the file (implementation in .cpp)
    void flush();

    // Method to write metadata/header to the log file (implementation in .cpp)
    static void logMetaData(std::ofstream &ofsLog); // Static as it doesn't depend on instance data

    // PUP routine (optional, uncomment if migration needed)
    // Need to inherit from PUP::able if pup routine is uncommented
    /*
    PUPable_decl(MemLog);
    MemLog(CkMigrateMessage *m) : PUP::able(m) {} // PUP constructor
    void pup(PUP::er& p) {
        PUP::able::pup(p);
        p | meTab;
        p | fileName;
    }
    */
};

#endif // MEMLOG_HINCLUDED 