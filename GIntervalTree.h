#ifndef E_INTERVAL_TREE
#define E_INTERVAL_TREE

#include "TStack.hh"
#include "GBase.h"
#include <math.h>
#include <limits.h>
#include <iostream>

//  The interval_tree.h and interval_tree.cc files contain code for 
//  interval trees implemented using red-black-trees as described in 
//  the book _Introduction_To_Algorithms_ by Cormen, Leisserson, 
//  and Rivest.  

//  CONVENTIONS:  
//                Function names: Each word in a function name begins with 
//                a capital letter.  An example function name is
//                CreateRedTree(a,b,c). Furthermore, each function name 
//                should begin with a capital letter to easily distinguish 
//                them from variables. 
//
//                Variable names: Each word in a variable name begins with 
//                a capital letter EXCEPT the first letter of the variable 
//                name.  For example, int newLongInt.  Global variables have 
//                names beginning with "g".  An example of a global 
//                variable name is gNewtonsConstant. 


#ifndef MAX_INT
#define MAX_INT INT_MAX // some architectures define INT_MAX not MAX_INT
#endif

// The Interval class is an Abstract Base Class.  This means that no
// instance of the Interval class can exist.  Only classes which
// inherit from the Interval class can exist.  Furthermore any class
// which inherits from the Interval class must define the member
// functions GetLowPoint and GetHighPoint.
//
// The GetLowPoint should return the lowest point of the interval and
// the GetHighPoint should return the highest point of the interval.  
/*
class GSeg {
public:
  GSeg(int vs=0, int ve=0):start(vs), end(vs) {}
  virtual ~GSeg() { }
  int start;
  int end;
  //virtual int GetLowPoint() const = 0;
  //virtual int GetHighPoint() const = 0;
  //virtual void Print() const;
};
*/

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

struct it_recursion_node {
public:
  /*  this structure stores the information needed when we take the */
  /*  right branch in searching for intervals but possibly come back */
  /*  and check the left branch as well. */

  GIntervalTreeNode * start_node;
  unsigned int parentIndex;
  int tryRightBranch;
} ;


class GIntervalTree {
private:
  unsigned int recursionNodeStackSize;
  it_recursion_node * recursionNodeStack;
  unsigned int currentParent;
  unsigned int recursionNodeStackTop;
protected:
  /*  A sentinel is used for root and for nil.  These sentinels are */
  /*  created when ITTreeCreate is called.  root->left should always */
  /*  point to the node which is the root of the tree.  nil points to a */
  /*  node which should always be black but has arbitrary children and */
  /*  parent and no key or info.  The point of using these sentinels is so */
  /*  that the root and nil nodes do not require special cases in the code */
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
  TemplateStack<GSeg *> * Enumerate(int low, int high) ;
  void CheckAssumptions() const;
};


#endif
