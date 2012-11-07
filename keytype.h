

#ifdef BIGKEYS
# if CMK_HAS_INT16
typedef CmiUInt16 KeyType;
# else
#error  "128-bit integer not supported."
# endif
#else
typedef CmiUInt8  KeyType;
#endif

