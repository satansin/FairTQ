/* ----------------------------------------------------------------------------
This header file includes class Advanced RegionTree declaration.
It is data structure to maintain the cells in query plane with hyperplanes
---------------------------------------------------------------------------- */

#ifndef _CELL_TREE_H_
#define _CELL_TREE_H_

#include "header.h"
#include "lp_lib.h"
#include "rtree.h"
#include "rentry.h"
#include "rnode.h"
#include "filemem.h"
#include "global.h"
#include <cmath>

#include "skyline.h"

extern int objCnt;

typedef struct rBound
{
	int minimum;
	int maximum;

	rBound()
	{
		minimum = 0;
		maximum = objCnt - 1;
	}
};

typedef struct cell
{
	unordered_map<long int, bool> IntersectHP;
	unordered_set<long int> AboveHP;
	unordered_set<long int> BelowHP;
	

	///////////////////////////////////////////
	unordered_set<long int> TopKSet;
	vector<vector<float> > extremePoint;
	unordered_map<long int, unordered_set<long int>> newAncestors;
	unordered_map<long int, unordered_set<long int>> newdescendants;
	vector<vector<long int > > pairMap;

	unordered_map<long int, unordered_set<long int>> combAncestors;
	unordered_map<long int, unordered_set<long int>> combDescendants;

	float Distance;

	///////////////////////////////////////////
	rBound lu;
	int rank;
	bool isPruned;

	unordered_set<long int> focaled;
	unordered_set<long int> ignoreset;

	cell* left;
	cell* right;

	cell()
	{
		left = NULL;
		right = NULL;
		isPruned = false;
		rank = 0;
		Distance= 1000;
	}

	cell(cell* copy)
	{
		if (copy->isPruned == false)
		{
			left = NULL;
			right = NULL;
			rank = copy->rank;
			isPruned = copy->isPruned;

			for (auto iter = copy->IntersectHP.begin(); iter != copy->IntersectHP.end(); iter++)
			{
				IntersectHP[iter->first] = iter->second;
			}
			for( int i = 0; i < copy->pairMap.size() ; i++)
			{
				pairMap.push_back(copy->pairMap[i]);
			}
			for (auto iter = copy->focaled.begin(); iter != copy->focaled.end(); iter++)
			{
				focaled.insert(*iter);
			}
			for (auto iter = copy->TopKSet.begin(); iter != copy->TopKSet.end(); iter++)
			{
				TopKSet.insert(*iter);
			}
			for (auto iter = copy->ignoreset.begin(); iter != copy->ignoreset.end(); iter++)
			{
				ignoreset.insert(*iter);
			}

			for(auto iter = copy->BelowHP.begin() ; iter != copy->BelowHP.end(); iter++)
			{
				BelowHP.insert(*iter);
			}
			for(auto iter = copy->AboveHP.begin() ; iter != copy->AboveHP.end(); iter++)
			{
				AboveHP.insert(*iter);
			}
			Distance = copy->Distance;

			for (auto iter = copy->newdescendants.begin(); iter != copy->newdescendants.end(); iter++)
			{
				newdescendants[iter->first] = iter->second;
			}

			for (auto iter = copy->newAncestors.begin(); iter != copy->newAncestors.end(); iter++)
			{
				newAncestors[iter->first] = iter->second;
			}

			for (auto iter = copy->combAncestors.begin(); iter != copy->combAncestors.end(); iter++)
			{
				combAncestors[iter->first] = iter->second;
			}

			for (auto iter = copy->combDescendants.begin(); iter != copy->combDescendants.end(); iter++)
			{
				combDescendants[iter->first] = iter->second;
			}
		} 
	}

	void appendto(cell* node)
	{
		for (auto iter = IntersectHP.begin(); iter != IntersectHP.end(); iter++)
		{
			node->IntersectHP[iter->first] = iter->second;
		}

		for (auto iter = AboveHP.begin(); iter != AboveHP.end(); iter++)
		{
			node->AboveHP.insert(*iter);
		}

		for (auto iter = BelowHP.begin(); iter != BelowHP.end(); iter++)
		{
			node->BelowHP.insert(*iter);
		}
		
		///////////////////////////////////////////
		for (auto iter = TopKSet.begin(); iter != TopKSet.end(); iter++)
        {
            node->TopKSet.insert(*iter);
        }

		///////////////////////////////////////////
		node->rank = rank;
		node->isPruned = isPruned;
	}

	void copyleaf(cell* node)
	{
		for (auto iter = IntersectHP.begin(); iter != IntersectHP.end(); iter++)
		{
			node->IntersectHP[iter->first] = iter->second;
		}
		for (auto iter = BelowHP.begin(); iter != BelowHP.end(); iter++)
		{
			node->BelowHP.insert(*iter);
		}
		for (auto iter = focaled.begin(); iter != focaled.end(); iter++)
		{
			node->focaled.insert(*iter);
		}
		for (auto iter = ignoreset.begin(); iter != ignoreset.end(); iter++)
		{
			node->ignoreset.insert(*iter);
		}
		///////////////////////////////////////////
		for (auto iter = TopKSet.begin(); iter != TopKSet.end(); iter++)
        {
            node->TopKSet.insert(*iter);
        }	
        //if(left == NULL && right == NULL)
        
        	for( int i = 0; i < extremePoint.size() ; i++)
        	{
        		node->extremePoint.push_back(extremePoint[i]);
        	}

        	for(int i = 0 ; i < pairMap.size(); i++ )
        	{
        		node->pairMap.push_back(pairMap[i]);
        	}
        

        for(auto iter = newAncestors.begin(); iter != newAncestors.end(); iter++ )
        {
        	node->newAncestors[iter->first] = iter->second;
        }
        for(auto iter = newdescendants.begin(); iter != newdescendants.end(); iter++)
        {
        	node->newdescendants[iter->first] = iter->second;
        }

        for(auto iter = combAncestors.begin(); iter != combAncestors.end(); iter++ )
        {
        	node->combAncestors[iter->first] = iter->second;
        }
        for(auto iter = combDescendants.begin(); iter != combDescendants.end(); iter++)
        {
        	node->combDescendants[iter->first] = iter->second;
        }
		///////////////////////////////////////////
        node->Distance = Distance;
		node->rank = rank;
		node->isPruned = isPruned;
	}

	void release()
	{
		IntersectHP.clear();
		unordered_map<long int, bool>().swap(IntersectHP);
		AboveHP.clear();
		unordered_set<long int>().swap(AboveHP);
		BelowHP.clear();
		unordered_set<long int>().swap(BelowHP);
		///////////////////////////////////////////
		TopKSet.clear();

		unordered_set<long int>().swap(TopKSet);
		extremePoint.clear();
		vector< vector<float> >().swap(extremePoint);
		pairMap.clear();
		vector< vector<long int> >().swap(pairMap);
	}
	
	double cellSize()
	{
		double ret = 0;
		ret += IntersectHP.size() * 5;
		ret += AboveHP.size() * 4;
		ret += BelowHP.size() * 4;
		///////////////////////////////////////////
		ret += TopKSet.size() * 4;
		ret += extremePoint.size()*4; //???? wrong
		if (left != NULL)
			ret += 4;
		else if (right != NULL)
			ret + 4;
		return ret*1.0 / MB;
	}

	void addNewAncestor(long int id, long int inc) {
		if (this->newAncestors.find(id) == this->newAncestors.end()) {
			this->newAncestors[id] = {inc};
		} else {
			this->newAncestors[id].insert(inc);
		}
	}
	void addNewDescendant(long int id, long int dec) {
		if (this->newdescendants.find(id) == this->newdescendants.end()) {
			this->newdescendants[id] = {dec};
		} else {
			this->newdescendants[id].insert(dec);
		}
	}
	void addCombAncestor(long int id, long int inc) {
		if (this->combAncestors.find(id) == this->combAncestors.end()) {
			this->combAncestors[id] = {inc};
		} else {
			this->combAncestors[id].insert(inc);
		}
	}
	void addCombDescendant(long int id, long int dec) {
		if (this->combDescendants.find(id) == this->combDescendants.end()) {
			this->combDescendants[id] = {dec};
		} else {
			this->combDescendants[id].insert(dec);
		}
	}
};

extern unordered_map<long int, cell*> cellID;
extern unordered_map<long int, RtreeNode*> ramTree;

typedef struct gNode
{
	long int rID;
	unordered_set<long int> rDominator;
	gNode(long int recordID) :rID(recordID){}; //!!modified!!there is an error in original code
};


struct cellCompare
{
	bool operator()(const cell* a, const cell* b) const
	{
		return a->Distance < b->Distance;
	}
};

struct cellpairCompare
{
	bool operator()(const pair<cell*, long int> &a, const pair<cell*, long int> &b) const
	{
		return a.first->rank < b.first->rank;
	}
};

typedef struct Score
{
	float min;
	float max;
	Score() 
	{
		min = 0;
		max = 1;
	}
};

class cellTree
{
public:
	cellTree();
	cellTree(cell* node);
	cellTree(int size);
	~cellTree();

	void releaseCell(cell* node);

	// void insert(vector<long long int> &hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult, unordered_set<long int> ancestors, vector<long int> descendants, vector<cell*>& subregions, long int focal, int k , float maxSP, vector<int> protectedAttribute, int numOfKind);
	void insert(vector<long long int> &hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult, const unordered_set<long int>& ancestors, vector<cell*>& subregions, long int focal, int k , float maxSP, const vector<int>& protectedAttribute, int numOfKind);
	bool isFeasible(unordered_map<long int, bool>& touchhs, long int hpid, bool sideindicator);
	void lpModel(lprec* lp, int& dimen, unordered_map<long int, bool>& touchhs);
	void addHP(lprec* model, long int hpid, bool sideindicator);

	void updateRank(cell* node, const int mink, long int hpid);
	void inserthp(long int & hpid, const int mink, cell* node, unordered_map<long int, bool>& touchhs, vector<cell*>& subregions);
	
	// added by Zheng
	void inserthp2(long long int & hpid, const int mink, cell* node, unordered_map<long int, bool>& touchhs, vector<cell*>& subregions, int k , float maxSP, const vector<int>& protectedAttribute, long int focal, int numOfKind);
	int pendingIntersection(vector<vector<float> >& extremePoint, long long  int hpid);
	void updateRelation(long long int & hpid, long int & focal, cell* node,  bool sideindicator);
	// vector<vector<double> > updateExtreme2d( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long int hpid, bool sideindicator, vector<vector<long int > >& pairMap);
	// vector<vector<double> > updateExtreme( vector<vector<float> >& extremePoint, unordered_map<long int, bool>& touchhs,long long int hpid, bool sideindicator, vector<vector<long int > >& pairMap);
	// vector<vector<double> > updateExtreme4d( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long int hpid, bool sideindicator, vector<vector<long int > >& pairMap);
	// vector<vector<double> > updateExtreme5d( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long int hpid, bool sideindicator, vector<vector<long int > >& pairMap);
	void updateExtreme2d( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long int hpid, bool sideindicator, vector<vector<long int > >& pairMap, vector<vector<float> > &ret);
	void updateExtreme( vector<vector<float> >& extremePoint, unordered_map<long int, bool>& touchhs,long long int hpid, bool sideindicator, vector<vector<long int > >& pairMap, vector<vector<float> > &ret);
	void updateExtreme4d( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long int hpid, bool sideindicator, vector<vector<long int > >& pairMap, vector<vector<float> > &ret);
	void updateExtreme5d( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long int hpid, bool sideindicator, vector<vector<long int > >& pairMap, vector<vector<float> > &ret);
	// added by Zheng


	void collectLeaf(vector<cell*>& leaves, const int& mink);
	void dfsTraversal(cell* node, cell* leaf, vector<cell*>& leaves);


	void opt_insert(vector<long int> &hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult);
	void opt_inserthp(long int & hpid, const int mink, cell* node, cell* all);
	int isDominatorInserted(unordered_set<long int>& rdominators, cell* all);

	void maintainDAG(unordered_set<long int>& skylines, unordered_set<long int>& removeSL, float* PG[], const int dimen);
	void markSingular(vector<cell*>& leaves, unordered_set<long int>& a_skylines, const int& mink, vector<cell*>& finalResult, unordered_set<long int>& singular);

	void updateCellTree(cell* node);

	int Lemma2(vector<cell*>& cells);

	void lpModelwithSIDELEN(lprec* lp, int& dimen, unordered_map<long int, bool>& touchhs);
	void addHPwithSIDELEN(lprec* model, long int hpid, bool sideindicator);
	void halfspace2polytope(cell* ret, string outfile, string hsfile);

	// look-ahead techniques

	void collectLeaf(vector<pair<cell*, long int>>& leaves, const int& mink);
	void dfsTraversal(cell* node, cell* leaf, vector<pair<cell*, long int>>& leaves);
	void markSingular(vector<pair<cell*, long int>>& leaves, unordered_set<long int>& a_skylines, const int& mink, vector<cell*>& finalResult, unordered_set<long int>& singular);
	void updateSkyline(unordered_set<long int>& skylines, vector<pair<cell*, long int>>& leaves, unordered_set<long int>& singular);

	void rankBound(vector<pair<cell*, long int>>& leaves, Rtree& rt, Point& a_pt, const int& mink, vector<cell*>& finalResult);
	void cellBound(cell* item, Rtree& rt, Point& a_pt, rBound& itemR, const int& mink);
	void recordScore(lprec* lp, RtreeNodeEntry* e, Score& rScore, bool& isLeafNode);
	void focalScore(lprec* lp, Point& pt, Score& pScore);

	void findCellMBR(lprec* lp, vector<float>& cl, vector<float>& cu);
	void approxiRecordScore(RtreeNodeEntry* e, Score& entryR, vector<float>& cl, vector<float>& cu, bool& isLeafNode);
	void rtreeRAM(Rtree& rt, unordered_map<long int, RtreeNode*>& ramTree);


	// optimized look-ahead techniques
	void rankBoundOpt(vector<pair<cell*, long int>>& leaves, Rtree& rt, Point& a_pt, const int& mink, vector<cell*>& finalResult);
	void cellBoundOpt(cell* item, Rtree& rt, Point& a_pt, rBound& itemR, const int& mink);

	double lpOptObj(lprec* lp, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, const bool& isMax);
	void fastOptObj(vector<float>& cu, vector<float>& cl, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, Score& objScore);

	void rankBoundOptFULL(vector<pair<cell*, long int>>& leaves, Rtree& rt, Point& a_pt, const int& mink, vector<cell*>& finalResult);
	void cellBoundOptFULL(cell* item, Rtree& rt, Point& a_pt, rBound& itemR, const int& mink);
	double lpOptObjFULL(lprec* lp, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, const bool& isMax);
	void fastOptObjFULL(vector<float>& cu, vector<float>& cl, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, Score& objScore);

	cell* root;
	double treeSize;

private:
	unordered_map<long int, gNode*> dagNode;
};
#endif
