#ifndef __SSEDEFS_H__
#define __SSEDEFS_H__

#include "cosmoType.h"

// The lane count is whatever native_simd gives for the compiler's target
// flags (-msse2 baseline, -mavx, -march=native), so it is not tied to an
// arch name here.
#include <experimental/simd>
namespace stdx = std::experimental;
typedef stdx::native_simd<cosmoType> SSEcosmoType;
typedef stdx::native_simd_mask<cosmoType> SSEMask;
static constexpr int SSE_VECTOR_WIDTH = static_cast<int>(SSEcosmoType::size());
#define FORCE_INPUT_LIST_PAD (SSE_VECTOR_WIDTH - 1)
#define SSELoad(where, arr, idx, field) \
  where(&arr[idx]field, stdx::element_aligned)
#define SSEStore(what, arr, idx, field) \
  what.copy_to(&arr[idx]field, stdx::element_aligned)

#endif
