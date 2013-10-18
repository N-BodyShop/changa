#ifndef __FASTMULTIPOLE_H__
#define __FASTMULTIPOLE_H__
class FMMCompute : public Compute {
  public:
  FMMCompute() : Compute(FMM) {}

  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi);

  int openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID, State *state);

  // book keeping on notifications
  void nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp);
  void nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int bucket);
  void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp, Tree::NodeKey &remoteBucket);

  void initState(State *state);
  void stateReady(State *, TreePiece *, int chunk, int start, int end);
//  void printUndlist(DoubleWalkState *state, int level, TreePiece *tp);
//  void printClist(DoubleWalkState *state, int level, TreePiece *tp);
  void reassoc(void *cE, int activeRung, Opt *o);
  State *getNewState(int d1, int d2);
  State *getNewState(int d1);
  State *getNewState();
  void freeState(State *state);
  void freeDoubleWalkState(DoubleWalkState *state);

  private:

  void addChildrenToCheckList(GenericTreeNode *node, int reqID, int chunk, int awi, State *s, CheckList &chklist, TreePiece *tp);
  void addNodeToInt(GenericTreeNode *node, int offsetID, DoubleWalkState *s);
  void addLocalParticlesToInt(GravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s);
  void addRemoteParticlesToInt(ExternalGravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s);
template <typename ParticleT>
  void addParticlesToChkList(ParticleT *part,
                               int nPart, int reqID,
                               CheckList &chklist, TreePiece *tp);
  DoubleWalkState *allocDoubleWalkState();
};

#endif
