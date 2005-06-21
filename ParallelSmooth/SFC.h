/** @file SFC.h
 Structures, defines, functions relating to space-filling curves,
 as used to build trees for particle data.
 */

#ifndef SFC_H
#define SFC_H

#include <iostream>

#include "Vector3D.h"
#include "OrientedBox.h"

namespace SFC {

typedef unsigned long long Key;

void printFloatBits(float f, std::ostream& os);
void printKeyBits(Key k, std::ostream& os);
void printIntBits(int k, std::ostream& os);

/** The very first possible key a particle can take on. */
const Key firstPossibleKey = static_cast<Key>(0);
/** The very last possible key a particle can take on. */
const Key lastPossibleKey = ~(static_cast<Key>(1) << 63);

/// Out of three floats, make a Morton order (z-order) space-filling curve key
// Cannot make this function inline, or g++ with -O2 or higher generates incorrect code
Key makeKey(Vector3D<float> v);

template <typename T, typename T2>
Key generateKey(const Vector3D<T>& v, const OrientedBox<T2>& boundingBox) {
	return makeKey((v - boundingBox.lesser_corner) / (boundingBox.greater_corner - boundingBox.lesser_corner) + Vector3D<float>(1, 1, 1));
}

template <typename T>
OrientedBox<T> cutBoxLeft(const OrientedBox<T>& box, const int axis) {
	OrientedBox<T> newBox = box;
	switch(axis) {
		case 0: //cut x
			newBox.greater_corner.x = (box.greater_corner.x + box.lesser_corner.x) / 2.0;
			break;
		case 1: //cut y
			newBox.greater_corner.y = (box.greater_corner.y + box.lesser_corner.y) / 2.0;				
			break;
		case 2: //cut z
			newBox.greater_corner.z = (box.greater_corner.z + box.lesser_corner.z) / 2.0;
			break;
	}
	return newBox;
}

template <typename T>
OrientedBox<T> cutBoxRight(const OrientedBox<T>& box, const int axis) {
	OrientedBox<T> newBox = box;
	switch(axis) {
		case 0: //cut x
			newBox.lesser_corner.x = (box.greater_corner.x + box.lesser_corner.x) / 2.0;
			break;
		case 1: //cut y
			newBox.lesser_corner.y = (box.greater_corner.y + box.lesser_corner.y) / 2.0;				
			break;
		case 2: //cut z
			newBox.lesser_corner.z = (box.greater_corner.z + box.lesser_corner.z) / 2.0;
			break;
	}
	return newBox;
}

/** Increment the floating point number until the last two bits
 of the mantissa are zero.
 */
float bumpLastTwoBits(float f, float direction);
void bumpBox(OrientedBox<float>& b, float direction);
void cubize(OrientedBox<float>& b);

} //close namespace SFC

#endif //SFC_H
