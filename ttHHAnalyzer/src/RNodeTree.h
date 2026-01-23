/*
 * RNodeTree.h
 *
 *  Created on: Dec 13, 2018
 *      Author: suyong
 */

#ifndef RNODETREE_H_
#define RNODETREE_H_

#include <vector>
#include <string>
#include "ROOT/RDataFrame.hxx"

using namespace std;
using namespace ROOT::RDF;


class RNodeTree
{
public:
	RNodeTree(RNode *rn);
	~RNodeTree();
	RNodeTree *addDaughter(RNode *rn, string asidx);
	RNodeTree *getParent(string asidx);

	// to get RNode of cut step d, in a linear chain
	RNodeTree *getRNodeTree(int d);
	// to get RNode of specific cut, if you have selections
	// that branch off at different stages, you must use this
	RNodeTree *getRNodeTree(string idx);
	RNode *getRNode();
	RNode *getRNode(int d);
	RNode *getRNode(string idx);
	void setRNode(RNode *rn);
	string getIndex();
	bool isLeaf();
	int getDaughtersize();
	void getRNodeLeafs(vector<RNodeTree *> &rntv);
	void Print();

private:
	void setParent(RNodeTree *rnt);
	RNode *thisnode;
	RNodeTree *parenttree;
	string _idx;
	vector<RNodeTree *> daughters;
};


#endif /* RNODETREE_H_ */
