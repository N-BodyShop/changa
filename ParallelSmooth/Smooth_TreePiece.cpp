/** \file Smooth_TreePiece.cpp
 \author Graeme Lufkin (gwl@u.washington.edu)
*/

#include "Tree.h"

#include "Smooth_TreePiece.h"
#include "SPH_Requests.h"

using std::string;
using std::priority_queue;

void Smooth_TreePiece::findSmoothingRadius(const int n, const CkCallback& cb) {
	numNeighbors = n;
	numComplete = 0;
	callback = cb;
	
	if(numParticles <= numNeighbors) {
		cerr << thisIndex << ": Smooth_TreePiece: Fatal: I don't have enough particles to do initial radius estimate" << endl;
		cb.send(0);
		return;
	}
	
	TreeRequestRequiringResponse* req;
	
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
		if(pbc.xPeriod && maxr > pbc.xPeriod / 2)
			cerr << thisIndex << ": Smooth_TreePiece: Not good, initial radius guess is larger than the x-axis period: " << maxr << endl;
		if(pbc.yPeriod && maxr > pbc.yPeriod / 2)
			cerr << thisIndex << ": Smooth_TreePiece: Not good, initial radius guess is larger than the y-axis period: " << maxr << endl;
		if(pbc.zPeriod && maxr > pbc.zPeriod / 2)
			cerr << thisIndex << ": Smooth_TreePiece: Not good, initial radius guess is larger than the z-axis period: " << maxr << endl;
		
		//the original serial number is the particle number
		req = new NeighborSearchRequest(myParticles[i].position, maxr, thisIndex, i, numNeighbors);
		req->startingNode = root->lookupKey();
		handleTreeRequest(req);
	}
}

void Smooth_TreePiece::performSmoothOperation(const int opID, const CkCallback& cb) {
	numComplete = 0;
	callback = cb;
	
	TreeRequest* req;
	for(int i = 1; i < numParticles + 1; ++i) {
		switch(opID) {
			case DensityOperation:
				req = new DensityRequest(myParticles[i].position, 2 * myParticles[i].smoothingRadius, thisIndex, i, myParticles[i].mass);
				break;
			case VelDivCurlOperation:
				req = new VelDivCurlRequest(myParticles[i].position, 2 * myParticles[i].smoothingRadius, thisIndex, i, myParticles[i].mass, myParticles[i].density, myParticles[i].velocity);
				break;
			case BallCountOperation:
				req = new BallCountRequest(myParticles[i].position, 2 * myParticles[i].smoothingRadius, thisIndex, i);
				break;
			case VelDispOperation:
				req = new VelDispRequest(myParticles[i].position, 2 * myParticles[i].smoothingRadius, thisIndex, i, myParticles[i].mass, myParticles[i].density, myParticles[i].velocity, myParticles[i].meanVelocity, myParticles[i].divv);
				break;
			case UserOperation:
				//req = dm->createUserTreeRequest(myParticles[i], thisIndex, i);
				break;
			default:
				cerr << "Aah!" << endl;
		}
		req->startingNode = root->lookupKey();
		handleTreeRequest(req);
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
		} else if(nodeIter->second->chareID == 0) { //I do have the node's parent!
			//send empty response, letting the requestor know that I don't have anything to contribute
			if(req->isResponseRequired()) {
				TreeRequestRequiringResponse* myReq = dynamic_cast<TreeRequestRequiringResponse *>(req);
				smoothTreePieces[myReq->requestingTreePiece].receiveResponse(myReq->createResponse(myReq->requestingSerialNumber));
			}
		} else {
			int realChare = (nodeIter->second->chareID < 0 ? -nodeIter->second->chareID - 1 : nodeIter->second->chareID - 1);
			//pass the buck
			if(verbosity) {
				cout << "Type of request: " << typeid(req).name() << endl;
				if(req->isResponseRequired())
					cout << thisIndex << ": Node " << req->startingNode << " from " << dynamic_cast<TreeRequestRequiringResponse *>(req)->requestingTreePiece << " is not really me, passing the buck to " << realChare << endl;
				else
					cout << thisIndex << ": Node " << req->startingNode << " is not really me, passing the buck to " << realChare << endl;
			}
			smoothTreePieces[realChare].handleTreeRequest(req);
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
	
	//check if the sphere really intersects my node
	if(!nodeIter->second->box.intersects(pbc, req->s)) {
		cerr << thisIndex << ": Smooth_TreePiece: Not Good: Why did I get this request, then?" << endl;
		delete req;
		return;
	}

	//traverse tree, handling criterion, building self-response and sending off requests
	//preOrderTraverse(nodeIter->second, req, resp);
	localFirstTraverse(nodeIter->second, req, resp);

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
					resp->receiveContribution(kernel, *req, *p, offset);
			}
		} else {
			//visit left child
			if(node->chareID < 0) { //left child is non-local
				//does left child really intersect?
				TreeNode leftChild(node, node->key);
				if(leftChild.box.intersects(pbc, req->s)) {
					if(req->isResponseRequired())
						pendingResponses[resp->serialNumber].numResponsesPending++;
					req->startingNode = node->leftChildLookupKey();
					smoothTreePieces[-node->chareID - 1].handleTreeRequest(req);
				}
			} else if(node->leftChild != 0)
				preOrderTraverse(node->leftChild, req, resp);
			
			//visit right child
			if(node->chareID > 0) { //right child is non-local
				//does right child really intersect?
				TreeNode rightChild(node, node->rightChildKey());
				if(rightChild.box.intersects(pbc, req->s)) {
					if(req->isResponseRequired())
						pendingResponses[resp->serialNumber].numResponsesPending++;
					req->startingNode = node->rightChildLookupKey();
					smoothTreePieces[node->chareID - 1].handleTreeRequest(req);
				}
			} else if(node->rightChild != 0)
				preOrderTraverse(node->rightChild, req, resp);
		}
	}
}

void Smooth_TreePiece::localFirstTraverse(TreeNode* node, TreeRequest* req, Response* resp) {
	if(node->box.intersects(pbc, req->s)) {
		if(node->isBucket()) {
			for(FullParticle* p = node->beginBucket; p != node->endBucket; ++p) {
				//calculate offset, pass to response->receiveContribution if p is in the sphere
				Vector3D<double> offset = pbc.offset(p->position, req->s.origin);
				if(offset.length() <= req->s.radius)
					resp->receiveContribution(kernel, *req, *p, offset);
			}
		} else {
			//visit left child
			if(node->chareID < 0) { //left child is non-local
				//visit right child first
				if(node->rightChild != 0)
					localFirstTraverse(node->rightChild, req, resp);
				//then visit non-local left child
				//does left child really intersect?
				TreeNode leftChild(node, node->key);
				if(leftChild.box.intersects(pbc, req->s)) {
					if(req->isResponseRequired())
						pendingResponses[resp->serialNumber].numResponsesPending++;
					req->startingNode = node->leftChildLookupKey();
					smoothTreePieces[-node->chareID - 1].handleTreeRequest(req);
				}
			} else {
				if(node->leftChild != 0)
					localFirstTraverse(node->leftChild, req, resp);
				//visit right child
				if(node->chareID > 0) { //right child is non-local
					//does right child really intersect?
					TreeNode rightChild(node, node->rightChildKey());
					if(rightChild.box.intersects(pbc, req->s)) {
						if(req->isResponseRequired())
							pendingResponses[resp->serialNumber].numResponsesPending++;
						req->startingNode = node->rightChildLookupKey();
						smoothTreePieces[node->chareID - 1].handleTreeRequest(req);
					}
				} else if(node->rightChild != 0)
					localFirstTraverse(node->rightChild, req, resp);
			}
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
			if(++numComplete == numParticles) //the radii have been found for all my particles
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
