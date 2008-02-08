
class State{
};

class DoubleWalkState : public State{
  public: 
  CkQ<Tree::NodeKey> chklist;
  CkQ<Tree::NodeKey> clist;
  CkQ<Tree::NodeKey> plist;
};

class NullState : public State{
};
