/** @file Smooth_TreePiece.cpp
 @author Graeme Lufkin (gwl@u.washington.edu)
*/

#include <fstream>
#include <algorithm>

#include "Tree.h"

#include "Smooth_TreePiece.h"
#include "SPH_Requests.h"

using std::string;
using std::priority_queue;

void Smooth_TreePiece::findSmoothingRadius(const int n, const CkCallback& cb) {
	numNeighbors = n;
	numComplete = 0;
	totalNumPending = numParticles;
	callback = cb;
	
	if(numParticles <= numNeighbors) {
		cerr << thisIndex << ": Smooth_TreePiece: Fatal: I don't have enough particles to do initial radius estimate" << endl;
		cb.send(0);
		return;
	}
	
	NeighborSearchRequest* req;
	NeighborSearchResponse* resp;
	
	double r;
	double maxr = 0;
	int lowIndex;

	//find the upper bound and create a request using it
	for(int i = 1; i < numParticles + 1; ++i) {
		if(i - numNeighbors / 2 < 1)
			lowIndex = 1;
		else if(i + (numNeighbors + 1) / 2 > numParticles)
			lowIndex = numParticles - numNeighbors;
		else
			lowIndex = i - numNeighbors / 2;
		maxr = 0;
		for(int j = lowIndex; j < lowIndex + numNeighbors; j++) {
			if(j == i)
				continue;
			r = pbc.distance(myParticles[j].position, myParticles[i].position);
			if(maxr == 0 || r > maxr)
				maxr = r;
		}
		
		//the original serial number is the particle number
		req = new NeighborSearchRequest(myParticles[i].position, maxr, thisIndex, i, numNeighbors);
		req->startingNode = root->lookupKey();
		
		resp = dynamic_cast<NeighborSearchResponse *>(req->createResponse(mySerialNumber));
		//file the pointer to Response using my serial number
		pendingResponses.insert(make_pair(mySerialNumber, PendingResponse(*req, resp)));
		req->requestingSerialNumber = mySerialNumber;
		mySerialNumber++;

		//traverse tree, handling criterion, building self-response and sending off requests
		localOnlyTraverse(root, req, resp);
		if(pbc.xPeriod && resp->radii.top() > pbc.xPeriod / 2)
			cerr << thisIndex << ": Smooth_TreePiece: Not good, initial radius guess is larger than the x-axis period: " << maxr << endl;
		if(pbc.yPeriod && resp->radii.top() > pbc.yPeriod / 2)
			cerr << thisIndex << ": Smooth_TreePiece: Not good, initial radius guess is larger than the y-axis period: " << maxr << endl;
		if(pbc.zPeriod && resp->radii.top() > pbc.zPeriod / 2)
			cerr << thisIndex << ": Smooth_TreePiece: Not good, initial radius guess is larger than the z-axis period: " << maxr << endl;
		remoteOnlyTraverse(root, req, resp);
		
		//send yourself the response you've built
		receiveResponse(resp);
		
		delete req;
	}
}

void Smooth_TreePiece::performSmoothOperation(const SmoothOperation op, const CkCallback& cb) {
	numComplete = 0;
	callback = cb;
	totalNumPending = numParticles;
	TreeRequest* req;
	for(int i = 1; i < numParticles + 1; ++i) {
		switch(op) {
			case Density:
				req = new DensityRequest(myParticles[i].position, 2 * myParticles[i].smoothingRadius, thisIndex, i, myParticles[i].mass);
				break;
			case VelocityDivCurl:
				req = new VelDivCurlRequest(myParticles[i].position, 2 * myParticles[i].smoothingRadius, thisIndex, i, myParticles[i].mass, myParticles[i].density, myParticles[i].velocity);
				break;
			case BallCount:
				req = new BallCountRequest(myParticles[i].position, 2 * myParticles[i].smoothingRadius, thisIndex, i);
				break;
			case VelocityDispersion:
				req = new VelDispRequest(myParticles[i].position, 2 * myParticles[i].smoothingRadius, thisIndex, i, myParticles[i].mass, myParticles[i].density, myParticles[i].velocity, myParticles[i].meanVelocity, myParticles[i].divv);
				break;
			case User:
				//req = dm->createUserTreeRequest(myParticles[i], thisIndex, i);
				break;
			default:
				cerr << "Aah!" << endl;
		}
		req->startingNode = root->lookupKey();
		handleTreeRequest(req);
	}
}

void Smooth_TreePiece::densityCutOperation(const SmoothOperation op, double minDensity, const CkCallback& cb) {
	numComplete = 0;
	callback = cb;
	minDensity = pow(10.0, minDensity);
	TreeRequest* req;
	totalNumPending = 0;
	for(int i = 1; i < numParticles + 1; ++i) {
		if(myParticles[i].density >= minDensity)
			totalNumPending++;
	}
	if(totalNumPending == 0) {
		contribute(0, 0, CkReduction::concat, callback);
		return;
	}
	for(int i = 1; i < numParticles + 1; ++i) {
		if(myParticles[i].density >= minDensity) {
			switch(op) {
				case Density:
					req = new DensityRequest(myParticles[i].position, 2 * myParticles[i].smoothingRadius, thisIndex, i, myParticles[i].mass);
					break;
				case VelocityDivCurl:
					myParticles[i].divv = 0;
					myParticles[i].curlv = myParticles[i].meanVelocity = Vector3D<double>(0, 0, 0);
					req = new VelDivCurlRequest(myParticles[i].position, 2 * myParticles[i].smoothingRadius, thisIndex, i, myParticles[i].mass, myParticles[i].density, myParticles[i].velocity);
					break;
				case BallCount:
					req = new BallCountRequest(myParticles[i].position, 2 * myParticles[i].smoothingRadius, thisIndex, i);
					break;
				case VelocityDispersion:
					req = new VelDispRequest(myParticles[i].position, 2 * myParticles[i].smoothingRadius, thisIndex, i, myParticles[i].mass, myParticles[i].density, myParticles[i].velocity, myParticles[i].meanVelocity, myParticles[i].divv);
					break;
				case User:
					//req = dm->createUserTreeRequest(myParticles[i], thisIndex, i);
					break;
				default:
					cerr << "Aah!" << endl;
			}
			req->startingNode = root->lookupKey();
			handleTreeRequest(req);
		}
	}
}

/// Receive a request to perform a calculation from someone.
void Smooth_TreePiece::handleTreeRequest(TreeRequest* req) {
	//look up starting node
	map<Key, TreeNode *>::iterator nodeIter = nodeLookup.find(req->startingNode);
	if(nodeIter == nodeLookup.end()) {	
		//I don't actually have this node, find who I think does
		nodeIter = nodeLookup.find(req->startingNode >> 1);
		if(nodeIter == nodeLookup.end()) {
			cerr << thisIndex << ": Smooth_TreePiece: Not Good: Got a request for a node I seriously don't have!" << endl;
		} else if(nodeIter->second->getType() == NonLocal) {
			cerr << "How does this sort of thing happen?" << endl;
			TreeRequestRequiringResponse* myReq = dynamic_cast<TreeRequestRequiringResponse *>(req);
			cerr << thisIndex << ": Chare " << myReq->requestingTreePiece << " thinks I have node ";
			for(int i = 0; i < 63; i++)
				cerr << (myReq->startingNode & (static_cast<Key>(1) << (62 - i)) ? 1 : 0);
			cerr << endl;
		} else { //I do have the node's parent!
			//the requested node is null; send empty response, letting the requestor know that I don't have anything to contribute
			if(req->isResponseRequired()) {
				TreeRequestRequiringResponse* myReq = dynamic_cast<TreeRequestRequiringResponse *>(req);
				smoothTreePieces[myReq->requestingTreePiece].receiveResponse(myReq->createResponse(myReq->requestingSerialNumber));
			}
		}
		delete req;
		return;
	}
	
	if(nodeIter->second->getType() == NonLocal) {
		//requested node isn't really me, although something thought so
		//forward the request to the right piece
		smoothTreePieces[nodeIter->second->chareID].handleTreeRequest(req);
		if(verbosity) {
			cout << "Type of request: " << typeid(req).name() << endl;
			if(req->isResponseRequired())
				cout << thisIndex << ": Node " << req->startingNode << " from " << dynamic_cast<TreeRequestRequiringResponse *>(req)->requestingTreePiece << " is not really me, passing the buck to " << nodeIter->second->chareID << endl;
			else
				cout << thisIndex << ": Node " << req->startingNode << " is not really me, passing the buck to " << nodeIter->second->chareID << endl;
		}
		delete req;
		return;
	}
	
	Response* resp = 0;
	if(req->isResponseRequired()) {
		TreeRequestRequiringResponse* myReq = dynamic_cast<TreeRequestRequiringResponse *>(req);
		resp = myReq->createResponse(mySerialNumber);
		//file the pointer to Response using my serial number
		pendingResponses.insert(make_pair(mySerialNumber, PendingResponse(*myReq, resp)));
		//make request point back to me
		myReq->requestingTreePiece = thisIndex;
		myReq->requestingSerialNumber = mySerialNumber;
	}
	
	//traverse tree, handling criterion, building self-response and sending off requests
	preOrderTraverse(nodeIter->second, req, resp);
	//localFirstTraverse(nodeIter->second, req, resp);
	//localOnlyTraverse(nodeIter->second, req, resp);
	//remoteOnlyTraverse(nodeIter->second, req, resp);

	//send yourself the response you've built
	if(req->isResponseRequired())
		receiveResponse(resp);
	
	mySerialNumber++;
	delete req;
}

/// Traverse my tree in an attempt to handle a request for calculation
void Smooth_TreePiece::preOrderTraverse(TreeNode* node, TreeRequest* req, Response* resp) {
	if(node->box.intersects(pbc, req->s)) {
		if(node->isBucket()) {
			for(FullParticle* p = node->beginBucket; p != node->endBucket; ++p) {
				//calculate offset, pass to response->receiveContribution if p is in the sphere
				Vector3D<double> offset = pbc.offset(p->position, req->s.origin);
				if(offset.length() <= req->s.radius)
					req->makeContribution(resp, p, offset);
			}
		} else if(node->isNonLocal()) {
			if(req->isResponseRequired())
				pendingResponses[resp->serialNumber].numResponsesPending++;
			req->startingNode = node->lookupKey();
			smoothTreePieces[node->chareID].handleTreeRequest(req);
		} else {
			if(node->leftChild != 0) //visit left child
				preOrderTraverse(node->leftChild, req, resp);
			if(node->rightChild != 0) //visit right child
				preOrderTraverse(node->rightChild, req, resp);
		}
	}
}

/// Walk the tree, only interacting with local nodes
void Smooth_TreePiece::localOnlyTraverse(TreeNode* node, TreeRequest* req, Response* resp) {
	if(node->getType() != NonLocal && node->box.intersects(pbc, req->s)) {
		if(node->isBucket()) {
			for(FullParticle* p = node->beginBucket; p != node->endBucket; ++p) {
				//calculate offset, pass to response->receiveContribution if p is in the sphere
				Vector3D<double> offset = pbc.offset(p->position, req->s.origin);
				if(offset.length() <= req->s.radius)
					req->makeContribution(resp, p, offset);
			}
		} else {
			if(node->leftChild != 0) //visit left child
				localOnlyTraverse(node->leftChild, req, resp);
			if(node->rightChild != 0) //visit right child
				localOnlyTraverse(node->rightChild, req, resp);
		}
	}
}

/// Walk the tree, only considering remote nodes, and sending the request to them
void Smooth_TreePiece::remoteOnlyTraverse(TreeNode* node, TreeRequest* req, Response* resp) {
	if(!node->isInternal() && node->box.intersects(pbc, req->s)) {
		if(node->getType() == NonLocal) {
			if(req->isResponseRequired())
				pendingResponses[resp->serialNumber].numResponsesPending++;
			req->startingNode = node->lookupKey();
			smoothTreePieces[node->chareID].handleTreeRequest(req);
		} else {
			if(node->leftChild != 0) //visit left child
				remoteOnlyTraverse(node->leftChild, req, resp);
			if(node->rightChild != 0) //visit right child
				remoteOnlyTraverse(node->rightChild, req, resp);
		}
	}
}

void Smooth_TreePiece::localFirstTraverse(TreeNode* node, TreeRequest* req, Response* resp) {
	if(node->box.intersects(pbc, req->s)) {
		if(node->isBucket()) {
			Vector3D<double> offset;
			for(FullParticle* p = node->beginBucket; p != node->endBucket; ++p) {
				//calculate offset, pass to response->receiveContribution if p is in the sphere
				offset = pbc.offset(p->position, req->s.origin);
				if(offset.length() <= req->s.radius)
					req->makeContribution(resp, p, offset);
			}
		} else if(node->isNonLocal()) {
			if(req->isResponseRequired())
				pendingResponses[resp->serialNumber].numResponsesPending++;
			req->startingNode = node->lookupKey();
			smoothTreePieces[node->chareID].handleTreeRequest(req);
		} else if(node->getType() == Boundary) {
			if(node->leftChild == 0)
				localFirstTraverse(node->rightChild, req, resp);
			else if(node->rightChild == 0)
				localFirstTraverse(node->leftChild, req, resp);
			else if(node->leftChild->isNonLocal()) {
				localFirstTraverse(node->rightChild, req, resp);
				localFirstTraverse(node->leftChild, req, resp);
			} else {
				localFirstTraverse(node->leftChild, req, resp);
				localFirstTraverse(node->rightChild, req, resp);
			}
		} else {
			if(node->leftChild != 0) //visit left child
				localFirstTraverse(node->leftChild, req, resp);
			if(node->rightChild != 0) //visit right child
				localFirstTraverse(node->rightChild, req, resp);
		}
	}
}

/// Handle a response returned by someone you requested a calculation from
void Smooth_TreePiece::receiveResponse(Response* resp) {
	if(resp->serialNumber < 0)
		cerr << thisIndex << ": What?  Negative serial number: " << resp->serialNumber << endl;
	//look up the pending response info
	map<int, PendingResponse>::iterator pending = pendingResponses.find(resp->serialNumber);
	if(pending == pendingResponses.end()) {
		cout << thisIndex << ": Smooth_TreePiece: Not Good: Hey, I should at least have my own entry!" << endl;
		cout << "Type of response: " << typeid(resp).name() << " Serial Number: " << resp->serialNumber << endl;
		return;
	}
	
	//merge in the new response, if it's not the starting response
	if(resp != pending->second.resp) {
		pending->second.resp->merge(resp);
		delete resp;
	}
		
	//if this is the last response expected, report back to calling chare
	if(pending->second.numResponsesPending == 1) {
		if(pending->second.requestingTreePiece == thisIndex) {
			//if the requesting chare is you, you have your answer, do something with it
			
			//modify the particle that made the initial request (starting serial number is particle index)
			pending->second.resp->updateInitialParticle(myParticles[pending->second.requestingSerialNumber]);
			
			//New: finishingCondition?
			if(++numComplete == totalNumPending) //the radii have been found for all my particles
				contribute(0, 0, CkReduction::concat, callback);
		} else {
			//prepare, then send
			pending->second.resp->prepareForSending();
			pending->second.resp->serialNumber = pending->second.requestingSerialNumber;
			smoothTreePieces[pending->second.requestingTreePiece].receiveResponse(pending->second.resp);
		}
		
		//delete the entry for this pending response
		delete pending->second.resp;
		pendingResponses.erase(pending);
	} else //decrement number of pending responses for this one
		pending->second.numResponsesPending--;

}

void Smooth_TreePiece::minmaxDensity(const CkCallback& cb) {
	callback = cb;
	double minmaxDensity[2];
	double logden;
	minmaxDensity[0] = minmaxDensity[1] = log10(myParticles[1].density);
	for(int i = 2; i < numParticles + 1; ++i) {
		logden = log10(myParticles[i].density);
		if(logden < minmaxDensity[0])
			minmaxDensity[0] = logden;
		if(logden > minmaxDensity[1])
			minmaxDensity[1] = logden;
	}
	contribute(2 * sizeof(double), minmaxDensity, minmaxReduction, cb);
}

void Smooth_TreePiece::makeDensityHistogram(const int numDensityBins, const double minDensity, const double maxDensity, const CkCallback& cb) {
	vector<int> densityHistogram(numDensityBins, 0);
	int bin;
	for(int i = 2; i < numParticles + 1; ++i) {
		bin = static_cast<int>(floor(numDensityBins * (log10(myParticles[i].density) - minDensity) / (maxDensity - minDensity)));
		if(bin >= numDensityBins)
			bin = numDensityBins - 1;
		else if(bin < 0)
			bin = 0;
		densityHistogram[bin]++;
	}
	contribute(densityHistogram.size() * sizeof(int), densityHistogram.begin(), CkReduction::sum_int, cb);
}

inline bool diskOrderCompare(const FullParticle& p1, const FullParticle& p2) {
	return p1.diskOrder < p2.diskOrder;	
}

void Smooth_TreePiece::saveInformation(const string& prefix, const CkCallback& cb) {
	mySerialNumber = 0;
	if(pendingResponses.size() != 0)
		cerr << thisIndex << " Ack! I still have responses!" << endl;
	
	//sort particles based on disk order
	sort(myParticles.begin() + 1, myParticles.end() - 1, diskOrderCompare);
	
	//open the output file
	fstream outgroup((prefix + ".radius").c_str(), ios::in | ios::out);
	
	//write all the radii in their respective locations
	int location;
	double value;
	for(int i = 1; i < numParticles + 1; ++i) {
		location = sizeof(int) + myParticles[i].diskOrder * sizeof(double);
		outgroup.seekp(location);
		outgroup.write(&(myParticles[i].smoothingRadius), sizeof(double));
	}
	outgroup.close();
		
	outgroup.open((prefix + ".density").c_str(), ios::in | ios::out);
	for(int i = 1; i < numParticles + 1; ++i) {
		location = sizeof(int) + myParticles[i].diskOrder * sizeof(double);
		outgroup.seekp(location);
		outgroup.write(&(myParticles[i].density), sizeof(double));
	}
	outgroup.close();

	outgroup.open((prefix + ".divv").c_str(), ios::in | ios::out);
	for(int i = 1; i < numParticles + 1; ++i) {
		location = sizeof(int) + myParticles[i].diskOrder * sizeof(double);
		outgroup.seekp(location);
		outgroup.write(&(myParticles[i].divv), sizeof(double));
	}
	outgroup.close();

	outgroup.open((prefix + ".curlv").c_str(), ios::in | ios::out);
	for(int i = 1; i < numParticles + 1; ++i) {
		location = sizeof(int) + myParticles[i].diskOrder * sizeof(Vector3D<double>);
		outgroup.seekp(location);
		outgroup.write(&(myParticles[i].curlv), sizeof(Vector3D<double>));
	}
	outgroup.close();

	outgroup.open((prefix + ".meanv").c_str(), ios::in | ios::out);
	for(int i = 1; i < numParticles + 1; ++i) {
		location = sizeof(int) + myParticles[i].diskOrder * sizeof(Vector3D<double>);
		outgroup.seekp(location);
		outgroup.write(&(myParticles[i].meanVelocity), sizeof(Vector3D<double>));
	}
	outgroup.close();

	outgroup.open((prefix + ".veldisp").c_str(), ios::in | ios::out);
	for(int i = 1; i < numParticles + 1; ++i) {
		location = sizeof(int) + myParticles[i].diskOrder * sizeof(double);
		outgroup.seekp(location);
		value = sqrt(myParticles[i].velDispSq);
		outgroup.write(&value, sizeof(double));
	}
	outgroup.close();
		
	outgroup.open((prefix + ".phase").c_str(), ios::in | ios::out);
	for(int i = 1; i < numParticles + 1; ++i) {
		location = sizeof(int) + myParticles[i].diskOrder * sizeof(double);
		outgroup.seekp(location);
		value = myParticles[i].density / pow(myParticles[i].velDispSq, 1.5);
		outgroup.write(&value, sizeof(double));
	}
	outgroup.close();
	
	outgroup.open((prefix + ".mach").c_str(), ios::in | ios::out);
	for(int i = 1; i < numParticles + 1; ++i) {
		location = sizeof(int) + myParticles[i].diskOrder * sizeof(double);
		outgroup.seekp(location);
		value = myParticles[i].meanVelocity.length() / sqrt(myParticles[i].velDispSq);
		outgroup.write(&value, sizeof(double));
	}
	outgroup.close();
	
	outgroup.open((prefix + ".speed").c_str(), ios::in | ios::out);
	for(int i = 1; i < numParticles + 1; ++i) {
		location = sizeof(int) + myParticles[i].diskOrder * sizeof(double);
		outgroup.seekp(location);
		value = myParticles[i].meanVelocity.length();
		outgroup.write(&value, sizeof(double));
	}
	outgroup.close();

	//re-sort particles back into key order
	sort(myParticles.begin() + 1, myParticles.end() - 1);
	
	contribute(0, 0, CkReduction::concat, cb);
}
