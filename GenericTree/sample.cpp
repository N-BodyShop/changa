//sample.cpp

#include <iostream>
#include <typeinfo>

#include "Vector3D.h"
#include "GravityTreeNode.h"
#include "Sphere.h"
#include "Space.h"

using namespace Tree;
using namespace std;

/// A geometrical factor used in doing multipole acceptance checks.
double opening_geometry_factor;
/// The opening angle for multipole acceptance checks.
double theta;

/// A simple particle definition
struct Particle {
	double mass;
	Vector3D<double> position;
};

/// An array of particles
Particle* myParticles;

/** Walking a tree to calculate gravity at a position.  In general this
 this function first checks if the multipole expansion is acceptable for
 this node.  If so, it calculates the acceleration due to this node,
 using the multipole information.  If the multipole is not accepted,
 we do a switch on node type, taking different action for different
 types of node.  For valid children, we recurse this function. */
void walkTree(GravityTreeNode* node, const Vector3D<double>& position, Vector3D<double>& acceleration) {
	Vector3D<double> r = node->moments.cm / node->moments.totalMass - position;
	double rsq = r.lengthSquared();
	double r_open = opening_geometry_factor * node->moments.radius / theta;
	if(rsq > r_open * r_open) {
		acceleration += node->moments.totalMass * r / rsq / sqrt(rsq);
	} else if(node->getType() == Bucket) {
		for(u_int64_t i = node->beginParticle; i < node->endParticle; ++i) {
			r = myParticles[i].position - position;
			rsq = r.lengthSquared();
			if(rsq != 0)
				acceleration += myParticles[i].mass * r / rsq / sqrt(rsq);
		}
	} else if(node->getType() == NonLocal) {
		/// What to do here depends on the implementation of the tree
		//unfilledRequests[mySerialNumber].numAdditionalRequests++;
		//streamingProxy[node->remoteIndex].fillRequestTree(req);
	} else {
		GenericTreeNode** childrenIterator = node->getChildren();
		for(unsigned int i = 0; i < node->numChildren(); ++i) {
			if(childrenIterator[i])
				walkTree(dynamic_cast<GravityTreeNode *>(childrenIterator[i]), position, acceleration);
		}
	}
}

/** Walking a tree to find particles in a sphere, for example like
 in an SPH calculation.  In general this this function first checks 
 for intersection of the node and the sphere of interest.  If they
 intersect, it switches on node type.  If a bucket, it checks for
 containmnent of the particles themselves.  For valid children, 
 we recurse this function. */
void walkTree(GenericTreeNode* node, const Sphere<double>& s) {
	if(Space::intersect(node->boundingBox, s)) {
		if(node->getType() == Bucket) {
			for(u_int64_t i = node->beginParticle; i < node->endParticle; ++i) {
				if(Space::contains(s, myParticles[i].position)) {
					//found a particle in the sphere, do something with it
				}
			}
		} else if(node->getType() == NonLocal) {
			/// What to do here depends on the implementation of the tree
			//unfilledRequests[mySerialNumber].numAdditionalRequests++;
			//streamingProxy[node->remoteIndex].fillRequestTree(req);
		} else {
			GenericTreeNode** childrenIterator = node->getChildren();
			for(unsigned int i = 0; i < node->numChildren(); ++i) {
				if(childrenIterator[i])
					walkTree(childrenIterator[i], s);
			}
		}
	}
}

int main() {
	
	cout << "Typeof GenericTreeNode: " << typeid(GenericTreeNode).name() << endl;
	cout << "Typeof BinaryTreeNode: " << typeid(BinaryTreeNode).name() << endl;
	cout << "Typeof GravityTreeNode: " << typeid(GravityTreeNode).name() << endl;
	
	BinaryTreeNode* bn = new BinaryTreeNode;	
	bn->rightChild = new BinaryTreeNode;
	bn->rightChild->parent = bn;
	bn->leftChild = new BinaryTreeNode;
	bn->leftChild->parent = bn;
	bn->leftChild->rightChild = new BinaryTreeNode;
	bn->leftChild->rightChild->parent = bn->leftChild;
	
	cout << "bn is " << typeid(*bn).name() << endl;
	cout << "bn is at: " << bn << endl;
	GenericTreeNode** iter = bn->getChildren();
	for(int i = 0; i < bn->numChildren(); ++i) {
		cout << "Child " << i << " is at " << iter[i] << endl;
		cout << "Typeof child: " << typeid(*(iter[i])).name() << endl;
		cout << "Child's parent is " << iter[i]->parent << endl;
		cout << "Typeof parent: " << typeid(*(iter[i]->parent)).name() << endl;
	}
	cout << "bn's left child : " << bn->leftChild << endl;
	cout << "bn's left child's parent : " << bn->leftChild->parent << endl;
	cout << "bn's right child: " << bn->rightChild << endl;
	cout << "bn's right child's parent: " << bn->rightChild->parent << endl;
	
	delete bn;
	
	cerr << "Done." << endl;
}

		
