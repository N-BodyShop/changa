/** @file RequestResponse.h
 @author Graeme Lufkin (gwl@u.washington.edu)
 */

#ifndef REQUESTRESPONSE_H
#define REQUESTRESPONSE_H

#include "pup.h"

#include "Vector3D.h"
#include "Sphere.h"
#include "SPH_Kernel.h"

#include "Particle.h"

//Declarations to include in an abstract PUP::able's body.
//  Abstract PUP::ables do not need def or reg.
#define PUPable_abstract(className) \
public:\
    friend inline void operator|(PUP::er &p,className &a) {a.pup(p);}\
    friend inline void operator|(PUP::er &p,className* &a) {\
        PUP::able *pa=a;  p(&pa);  a=(className *)pa;\
    }

/** A base class representing a response built to answer a tree request. */
class Response : public PUP::able {
public:
	
	/// The serial number used to file this response
	int serialNumber;
	
	Response(int sn = 0) : serialNumber(sn) { }
	
	Response(CkMigrateMessage* m) : PUP::able(m) { }
	
	virtual ~Response() { };

	virtual void merge(const Response* r) = 0;

	virtual void prepareForSending() { }
	
	virtual void updateInitialParticle(FullParticle& p) { }
	
	virtual void pup(PUP::er& p) {
		PUP::able::pup(p);
		p(serialNumber);
	}
	PUPable_abstract(Response);
};

/** A base class representing a request made of the tree. */
class TreeRequest : public PUP::able {
	static const double spline_PI = 3.14159265358979323846;
public:

	/** The ball to search in.  The radius of this sphere will be twice the smoothing radius
	 of the particle it belongs to. */
	Sphere<double> s;
	
	Key startingNode;

	TreeRequest() { }

	TreeRequest(const Vector3D<double>& pos, const double r) : s(pos, r) { }
	
	TreeRequest(CkMigrateMessage* m) : PUP::able(m) { }
	
	virtual ~TreeRequest() { };
	
	virtual bool isResponseRequired() = 0;
	
	virtual void makeContribution(Response* resp, FullParticle* p, const Vector3D<double>& offset) = 0;

	inline double kernelEvaluate(double r, double h) const {
		double q = r / h;
		if(q < 1)
			return (1 - 1.5 * q * q + 0.75 * q * q * q) / spline_PI / h / h / h;
		else
			return 0.25 * (2 - q) * (2 - q) * (2 - q) / spline_PI / h / h / h;
	}
	
	inline double kernelEvaluateGradient(double r, double h) const {
		double q = r / h;
		if(q < 1)
			return (0.75 * q - 1) * 3 / spline_PI / h / h / h / h / h;
		else
			return (-0.25 * q - 1 / q + 1) * 3 / spline_PI / h / h / h / h / h;
	}

	virtual void pup(PUP::er& p) {
		PUP::able::pup(p);
		p(s.origin.x);
		p(s.origin.y);
		p(s.origin.z);
		p(s.radius);
		p(startingNode);
	}
	PUPable_abstract(TreeRequest);
};

class TreeRequestRequiringResponse : public TreeRequest {
public:

	/// The array index of the TreePiece which made this request, and expects a response
	int requestingTreePiece;
	int requestingSerialNumber;

	TreeRequestRequiringResponse(const Vector3D<double>& pos, const double r, const int tp, const int sn) : TreeRequest(pos, r), requestingTreePiece(tp), requestingSerialNumber(sn) { }

	TreeRequestRequiringResponse(CkMigrateMessage* m) : TreeRequest(m) { }

	bool isResponseRequired() {
		return true;
	}
	
	virtual Response* createResponse(int serialNumber) = 0;
	
	virtual void pup(PUP::er& p) {
		TreeRequest::pup(p);
		p(requestingTreePiece);
		p(requestingSerialNumber);
	}
	PUPable_abstract(TreeRequestRequiringResponse);
};

/** The information stored by a TreePiece when it is waiting for responses from itself and 
 any other TreePieces it may ask for responses from. */
class PendingResponse {
public:
	
	int requestingTreePiece;
	int requestingSerialNumber;
	int numResponsesPending;
	Response* resp;
	
	PendingResponse() { }
	
	PendingResponse(const TreeRequestRequiringResponse& req, Response* r) : requestingTreePiece(req.requestingTreePiece), requestingSerialNumber(req.requestingSerialNumber), numResponsesPending(1), resp(r) { }
};

typedef TreeRequest* TreeRequest_Pointer;
typedef Response* Response_Pointer;

#endif //REQUESTRESPONSE_H
