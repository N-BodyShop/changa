/** @file SPH_Requests.h
 @author Graeme Lufkin (gwl@u.washington.edu)
 */

#ifndef SPH_REQUESTS_H
#define SPH_REQUESTS_H

#include <vector>
#include <queue>

#include "SPH_Kernel.h"
#include "RequestResponse.h"

#include "pup_stl.h"

using std::vector;
using std::priority_queue;

class BallCountResponse : public Response {
	int count;
	friend class BallCountRequest;
public:
		
	BallCountResponse(int serialNumber) : Response(serialNumber), count(0) { }
	
	BallCountResponse(CkMigrateMessage* m) : Response(m) { }
	
	virtual void merge(const Response* r) {
		count += dynamic_cast<const BallCountResponse *>(r)->count;
	}
		
	void updateInitialParticle(FullParticle& p) {
		p.neighborCount = count;
	}
	
	virtual void pup(PUP::er& p) {
		Response::pup(p);
		p(count);
	}
	PUPable_decl(BallCountResponse);
};

/** Count the number of particles in a ball. */
class BallCountRequest : public TreeRequestRequiringResponse {
public:

	BallCountRequest(const Vector3D<double>& pos, const double r, const int tp, const int sn) : TreeRequestRequiringResponse(pos, r, tp, sn) { }

	BallCountRequest(CkMigrateMessage* m) : TreeRequestRequiringResponse(m) { }

	virtual Response* createResponse(int serialNumber) {
		return new BallCountResponse(serialNumber);
	}
	
	virtual void makeContribution(Response* resp, FullParticle* q, const Vector3D<double>& offset) {
		dynamic_cast<BallCountResponse *>(resp)->count++;
	}
	
	PUPable_decl(BallCountRequest);
};


class NeighborSearchResponse : public Response {
	//properties of the initiating particle needed for the calculation
	const unsigned int numNeighbors;
	//results we're building to give back to the initiating particle
	vector<double> radii_vector;
	priority_queue<double> radii;
	friend class NeighborSearchRequest;

public:
	
	NeighborSearchResponse(const int sn, const unsigned int n) : Response(sn), numNeighbors(n) { }
	
	NeighborSearchResponse(CkMigrateMessage* m) : Response(m) { }

	virtual void merge(const Response* r) {
		const NeighborSearchResponse* resp = dynamic_cast<const NeighborSearchResponse *>(r);
		for(vector<double>::const_iterator iter = resp->radii_vector.begin(); iter != resp->radii_vector.end(); ++iter)
			radii.push(*iter);
		while(radii.size() > numNeighbors)
			radii.pop();
	}
	
	virtual void prepareForSending() {
		radii_vector.reserve(radii.size());
		while(!radii.empty()) {
			radii_vector.push_back(radii.top());
			radii.pop();
		}
	}
		
	virtual void updateInitialParticle(FullParticle& p) {
		p.smoothingRadius = radii.top() / 2;
	}
	
	virtual void pup(PUP::er& p) {
		Response::pup(p);
		p | radii_vector;
	}
	PUPable_decl(NeighborSearchResponse);
};

/** Given an upper bound, find the radius of a ball that contains a specified number of particles. */
class NeighborSearchRequest : public TreeRequestRequiringResponse {
	unsigned int numNeighbors;
public:

	NeighborSearchRequest(const Vector3D<double>& pos, const double r, const int tp, const int sn, const unsigned int n) : TreeRequestRequiringResponse(pos, r, tp, sn), numNeighbors(n) { }

	NeighborSearchRequest(CkMigrateMessage* m) : TreeRequestRequiringResponse(m) { }

	virtual Response* createResponse(int serialNumber) {
		return new NeighborSearchResponse(serialNumber, numNeighbors);
	}
	
	virtual void makeContribution(Response* resp, FullParticle* q, const Vector3D<double>& offset) {
		NeighborSearchResponse* nsresp = dynamic_cast<NeighborSearchResponse *>(resp);
		nsresp->radii.push(offset.length());
		//remove any extra radii
		while(nsresp->radii.size() > numNeighbors)
			nsresp->radii.pop();
		//possibly shrink the ball
		if(nsresp->radii.size() == numNeighbors)
			s.radius = nsresp->radii.top();
	}

	virtual void pup(PUP::er& p) {
		TreeRequestRequiringResponse::pup(p);
		p(numNeighbors);
	}
	PUPable_decl(NeighborSearchRequest);
};

class DensityResponse : public Response {
	//results we're building to give back to the initiating particle
	double density;
	friend class DensityRequest;
public:
	
	DensityResponse(const int sn) : Response(sn), density(0) { }
	
	DensityResponse(CkMigrateMessage* m) : Response(m) { }

	virtual void merge(const Response* r) {
		density += dynamic_cast<const DensityResponse *>(r)->density;
	}

	virtual void updateInitialParticle(FullParticle& p) {
		p.density += density;
	}
	
	virtual void pup(PUP::er& p) {
		Response::pup(p);
		p(density);
	}
	PUPable_decl(DensityResponse);
};

/** Calculate the density in a ball. */
class DensityRequest : public TreeRequestRequiringResponse {
	double mass;
public:

	DensityRequest(const Vector3D<double>& pos, const double r, const int tp, const int sn, const double m) : TreeRequestRequiringResponse(pos, r, tp, sn), mass(m) { }

	DensityRequest(CkMigrateMessage* m) : TreeRequestRequiringResponse(m) { }

	virtual Response* createResponse(int serialNumber) {
		return new DensityResponse(serialNumber);
	}
	
	virtual void makeContribution(Response* resp, FullParticle* q, const Vector3D<double>& offset) {
		double factor = 0.5 * kernelEvaluate(offset.length(), s.radius / 2);
		//gather part
		dynamic_cast<DensityResponse *>(resp)->density += factor * q->mass;
		//scatter part
		q->density += factor * mass;		
	}
	
	virtual void pup(PUP::er& p) {
		TreeRequestRequiringResponse::pup(p);
		p(mass);
	}
	PUPable_decl(DensityRequest);
};

class VelDivCurlResponse : public Response {
	//results we're building to give back to the initiating particle
	double divv;
	Vector3D<double> curlv;
	Vector3D<double> meanVelocity;
	friend class VelDivCurlRequest;
public:
	
	VelDivCurlResponse(const int sn) : Response(sn), divv(0), curlv(0, 0, 0), meanVelocity(0, 0, 0) { }
	
	VelDivCurlResponse(CkMigrateMessage* m) : Response(m) { }

	virtual void merge(const Response* r) {
		divv += dynamic_cast<const VelDivCurlResponse *>(r)->divv;
		curlv += dynamic_cast<const VelDivCurlResponse *>(r)->curlv;
		meanVelocity += dynamic_cast<const VelDivCurlResponse *>(r)->meanVelocity;
	}
		
	virtual void updateInitialParticle(FullParticle& p) {
		p.divv += divv;
		p.curlv += curlv;
		p.meanVelocity += meanVelocity;
	}
	
	virtual void pup(PUP::er& p) {
		Response::pup(p);
		p(divv);
		p | curlv;
		p | meanVelocity;
	}
	PUPable_decl(VelDivCurlResponse);
};

/** Calculate the divergence, curl and mean of the velocity in a ball. */
class VelDivCurlRequest : public TreeRequestRequiringResponse {
	double mass;
	double density;
	Vector3D<double> velocity;
public:

	VelDivCurlRequest(const Vector3D<double>& pos, const double r, const int tp, const int sn, const double m, const double d, const Vector3D<double> vel) : TreeRequestRequiringResponse(pos, r, tp, sn), mass(m), density(d), velocity(vel) { }

	VelDivCurlRequest(CkMigrateMessage* m) : TreeRequestRequiringResponse(m) { }

	virtual Response* createResponse(int serialNumber) {
		return new VelDivCurlResponse(serialNumber);
	}
	
	virtual void makeContribution(Response* resp, FullParticle* q, const Vector3D<double>& offset) {
		VelDivCurlResponse* p = dynamic_cast<VelDivCurlResponse *>(resp);
		double common = 0.5 * kernelEvaluateGradient(offset.length(), s.radius / 2);
		double divfactor = common * dot(velocity - q->velocity, offset);
		//gather part
		p->divv += divfactor * q->mass / density;
		//scatter part
		q->divv += divfactor * mass / q->density;
	
		Vector3D<double> curlfactor = common * cross(velocity - q->velocity, offset);
		//gather part
		p->curlv += curlfactor * q->mass / density;
		//scatter part
		q->curlv += curlfactor * mass / q->density;
	
		double factor = 0.5 * kernelEvaluate(offset.length(), s.radius / 2);
		//gather part
		p->meanVelocity += factor * q->mass / q->density * q->velocity;
		//scatter part
		q->meanVelocity += factor * mass / density * velocity;
	}
	
	virtual void pup(PUP::er& p) {
		TreeRequestRequiringResponse::pup(p);
		p(mass);
		p(density);
		p | velocity;
	}
	PUPable_decl(VelDivCurlRequest);
};

class VelDispResponse : public Response {
	//results we're building to give back to the initiating particle
	double velDispSq;
	friend class VelDispRequest;
public:
	
	VelDispResponse(const int sn) : Response(sn), velDispSq(0) { }
	
	VelDispResponse(CkMigrateMessage* m) : Response(m) { }

	virtual void merge(const Response* r) {
		velDispSq += dynamic_cast<const VelDispResponse *>(r)->velDispSq;
	}
		
	virtual void updateInitialParticle(FullParticle& p) {
		p.velDispSq += velDispSq;
	}
	
	virtual void pup(PUP::er& p) {
		Response::pup(p);
		p(velDispSq);
	}
	PUPable_decl(VelDispResponse);
};

/** Calculate the dispersion of the velocity in a ball. */
class VelDispRequest : public TreeRequestRequiringResponse {
	double mass;
	double density;
	Vector3D<double> velocity;
	Vector3D<double> meanVelocity;
	double divv;
public:

	VelDispRequest(const Vector3D<double>& pos, const double r, const int tp, const int sn, const double m, const double d, const Vector3D<double> vel, const Vector3D<double> mVel, const double div) : TreeRequestRequiringResponse(pos, r, tp, sn), mass(m), density(d), velocity(vel), meanVelocity(mVel), divv(div) { }

	VelDispRequest(CkMigrateMessage* m) : TreeRequestRequiringResponse(m) { }

	virtual Response* createResponse(int serialNumber) {
		return new VelDispResponse(serialNumber);
	}
	
	virtual void makeContribution(Response* resp, FullParticle* q, const Vector3D<double>& offset) {	
		double factor = 0.5 * kernelEvaluate(offset.length(), s.radius / 2);
		//gather part
		dynamic_cast<VelDispResponse *>(resp)->velDispSq += factor * q->mass / q->density * (q->velocity - meanVelocity - divv * offset).lengthSquared();
		//scatter part
		q->velDispSq += factor * mass / density * (q->velocity - meanVelocity - divv * offset).lengthSquared();
	}
	
	virtual void pup(PUP::er& p) {
		TreeRequestRequiringResponse::pup(p);
		p(mass);
		p(density);
		p | velocity;
		p | meanVelocity;
		p(divv);
	}
	PUPable_decl(VelDispRequest);
};

#endif //SPH_REQUESTS_H
