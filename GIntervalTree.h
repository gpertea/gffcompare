#ifndef E_INTERVAL_TREE
#define E_INTERVAL_TREE

//#include "TStack.hh"
#include "GBase.h"
#include "GVec.hh"
//#include <math.h>
//#include <limits.h>
//#include <iostream>

// This is an interval tree implementation based on red-black-trees
// as described in the book _Introduction_To_Algorithms_ by Cormen, Leisserson, and Rivest.

#ifndef MAX_INT
#define MAX_INT INT_MAX // some architectures define INT_MAX not MAX_INT
#endif

class GIntervalTreeNode {
  friend class GIntervalTree;
protected:
  GSeg* storedInterval;
  int key;
  int high;
  int maxHigh;
  int red; /* if red=0 then the node is black */
  GIntervalTreeNode* left;
  GIntervalTreeNode* right;
  GIntervalTreeNode* parent;
public:
  void Print(GIntervalTreeNode*,
	     GIntervalTreeNode*) const;
  GIntervalTreeNode():storedInterval(NULL), key(0), high(0),maxHigh(0),red(0),
		  left(NULL), right(NULL), parent(NULL) {}
  GIntervalTreeNode(GSeg * newInterval): storedInterval (newInterval),
      key(newInterval->start), high(newInterval->end) ,
      maxHigh(high), red(0), left(NULL), right(NULL), parent(NULL) {  }
  ~GIntervalTreeNode() {}
};

struct G_ITRecursionNode {
public:
  // this structure stores the information needed when we take the
  // right branch in searching for intervals but possibly come back
  // and check the left branch as well.
  GIntervalTreeNode * start_node;
  unsigned int parentIndex;
  int tryRightBranch;
} ;

class GIntervalTree {
private:
  unsigned int recursionNodeStackSize;
  G_ITRecursionNode * recursionNodeStack;
  unsigned int currentParent;
  unsigned int recursionNodeStackTop;
protected:
  // A sentinel is used for root and for nil.  root->left should always
  // point to the node which is the root of the tree.  nil points to a
  // node which should always be black but has arbitrary children and
  // parent and no key or info; These sentinels are used so
  // that the root and nil nodes do not require special treatment in the code
  GIntervalTreeNode* root;
  GIntervalTreeNode* nil;
  void LeftRotate(GIntervalTreeNode*);
  void RightRotate(GIntervalTreeNode*);
  void TreeInsertHelp(GIntervalTreeNode*);
  void TreePrintHelper(GIntervalTreeNode*) const;
  void FixUpMaxHigh(GIntervalTreeNode*);
  void DeleteFixUp(GIntervalTreeNode*);
  void CheckMaxHighFields(GIntervalTreeNode*) const;
  int CheckMaxHighFieldsHelper(GIntervalTreeNode* y,
			const int currentHigh,
			int match) const;
public:
  GIntervalTree();
  ~GIntervalTree();
  void Print() const;
  GSeg* DeleteNode(GIntervalTreeNode*);
  GIntervalTreeNode * Insert(GSeg*);
  GIntervalTreeNode * GetPredecessorOf(GIntervalTreeNode*) const;
  GIntervalTreeNode * GetSuccessorOf(GIntervalTreeNode*) const;
  GVec<GSeg*> * Enumerate(int low, int high) ;
};


#endif
