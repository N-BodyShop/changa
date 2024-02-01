#ifndef RAND_HINCLUDED
#define RAND_HINCLUDED

/** @brief basic random number generator
 *  This is implemented as an object so that it can easily be used in
 *  threads: one Rand instance per thread.
 */
class Rand 
{
 private:
    u_int64_t u,v,w;            /* State holders */
 public:
    Rand(u_int64_t j = 0            /* Seed */
         ) : v(4101842887655102017LL), w(1) {
        u = j ^ v; int64();
        v = u; int64();
        w = v; int64();
    }
    /** @brief return 64 bit random number
     */
    inline u_int64_t int64() {
        u = u * 2862933555777941757LL + 7046029254386353087LL;
        v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
        w = 4294957665U*(w & 0xffffffff) + (w >> 32);
        u_int64_t x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
        return (x + v) ^ w;
    }
    /** @brief return random double 0<x<1
     */
    inline double dbl() {
        const double iRANDMAX = 5.42101086242752217E-20;
        return iRANDMAX * int64();
    }
    /** @brief return 32 bit random number
     */
    inline u_int32_t int32() { return (u_int32_t)int64(); }
#ifdef __CHARMC__
    inline void pup(PUP::er &p) {
        p|u;
        p|v;
        p|w;
    }
#endif
};

#endif // RAND_HINCLUDED
