/*
 * RNodeTree.cpp
 *
 *  Created on: Dec 13, 2018
 *      Author: suyong
 */
#include "RNodeTree.h"

RNodeTree::RNodeTree(RNode *rn)
:thisnode(rn), parenttree(NULL), _idx(""), daughters(0)
{

}

RNodeTree::~RNodeTree()
{

}

RNodeTree *RNodeTree::addDaughter(RNode *rn, string asidx)
{
	/*
	RNodeTree *rootnode=this;
	while(rootnode->parenttree!=NULL)
	{
		rootnode = parenttree;
	}
	RNodeTree *parentnode = rootnode->getRNodeTree(asidx.substr(0, asidx.length()-1));
	*/

	RNodeTree *parentnode = this->getParent(asidx);

	RNodeTree *rnt = NULL;
	if (parentnode != NULL)
	{
		if ( int(parentnode->daughters.size()) == std::stoi(asidx.substr(asidx.length()-1, 1)))
		rnt = new RNodeTree(rn);
		rnt->setParent(parentnode);
		rnt->_idx = asidx;
		parentnode->daughters.push_back(rnt);
	}
	else
	{
		std::cout << "Problem adding node "+asidx << std::endl;
	}
	return rnt;
}

RNodeTree *RNodeTree::getParent(string asidx)
{
	RNodeTree *rootnode=this;
	while(rootnode->parenttree!=NULL)
	{
		rootnode = parenttree;
	}
	RNodeTree *parentnode = rootnode->getRNodeTree(asidx.substr(0, asidx.length()-1));

	return parentnode;
}

void RNodeTree::setParent(RNodeTree *rnt)
{
	parenttree = rnt;
}

RNodeTree *RNodeTree::getRNodeTree(int d)
{
	RNodeTree *result = NULL;

	RNodeTree *rootnode=this;
	while(rootnode->parenttree!=NULL)
	{
		rootnode = parenttree;
	}

	if (d==0) result = rootnode;
	else {
		string idx = "";
		for (auto i=0; i<d; i++)
		{
			idx += "0";
		}
		result = rootnode->getRNodeTree(idx);
	}

	return result;
}

RNodeTree *RNodeTree::getRNodeTree(string idx)
{
	RNodeTree *result = NULL;

	// recursively call
	if (idx.size()==0) result = this;
	else
	{
		int nodenumber=std::stoi(idx.substr(0, 1));
		RNodeTree *daughter = daughters[nodenumber];
		result = daughter->getRNodeTree(idx.substr(1));
	}

	return result;
}

RNode *RNodeTree::getRNode()
{
	return this->thisnode;
}

RNode *RNodeTree::getRNode(int d)
{
	RNode *result = NULL;

	RNodeTree *rnt = getRNodeTree(d);
	if (rnt != NULL) result = rnt->thisnode;
	return result;
}


RNode *RNodeTree::getRNode(string idx)
{
	RNode *result = NULL;

	RNodeTree *rnt = getRNodeTree(idx);
	if (rnt != NULL) result = rnt->thisnode;
	return result;
}

string RNodeTree::getIndex(){
	return _idx;
}

void RNodeTree::setRNode(RNode *rn)
{
	this->thisnode = rn;
}

bool RNodeTree::isLeaf()
{
	if (this->getDaughtersize()==0) return true;
	else return false;
}

int RNodeTree::getDaughtersize()
{
	return daughters.size();
}

void RNodeTree::getRNodeLeafs(vector<RNodeTree *> &rntv)
{
	// dive into the daughters until leaf is met
	if (this->isLeaf())
	{
		rntv.push_back(this);
	}
	else
	{
		for (auto adaughter: daughters)
		{
			adaughter->getRNodeLeafs(rntv);
		}
	}

}

void RNodeTree::Print()
{
	cout << this->_idx << endl;
	for (auto ad : daughters)
	{
		ad->Print();
	}
}

