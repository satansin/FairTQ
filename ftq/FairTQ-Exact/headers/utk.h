#ifndef _UTOK_H_
#define _UTOK_H_


#include "header.h"
#include "cellTree.h"
#include "skyline.h"
//#include "Array.h"
extern unordered_map<long long int, long long int> RecordIDtoHalfPlaneID;
extern vector<vector<float>> HalfSpaces; // halfspace 
extern unordered_map<long int, RtreeNode*> ramTree; // load Rtree to main-memory

struct daGraphCompare
{
	bool operator() (const pair<long int, unordered_set<long int>> a, const  pair<long int, unordered_set<long int>> b) const
	{
		return a.second.size() > b.second.size();
	}
};

class UTK
{

public:
	UTK();
	~UTK();

	// pivot of region
	void pivotRegion(vector<float>& R, vector<float>& pivot);
	vector<double> find_feasible(vector<vector<float> > hyperplane);
	vector<double> newHope( unordered_map<long int, bool>& utility_hyperplane);


	// r-dominance
	bool isRdominated(const int  dim, vector<float>& R, float focal[], vector<float>& entry, bool& fDe);
	vector<long int> TopKQuery(vector<float>& pivot, int dim, float* PG[], int k);
	vector<long int> TopKQuery(vector<double>& pivot, int dim, float* PG[], int k);

	// count region dominator
	bool countRegionDominator(int dimen, float pt[], vector<long int>& rskyband, float* PG[], vector<float>& R, const int k, unordered_set<long int>& dominators);

	// order score for r skyband
	float orderScore(vector<float>& pivot, float entry[], int dimen);

	// maintain dominace graph
	void maintainGraph(const int recordID, float* PG[], vector<float>& region, int dim);

	// r-skyband and r-dominace graph
	void rskyband(vector<float>& region, const int dimen, Rtree& a_rtree, vector<long int>& rskyband, float* PG[], const int k);
	void onionToRskyband(vector<float>& region, const int dimen, Rtree& a_rtree, vector<long int>& rskyband, float* PG[], const int k, vector<long int>klayers );
	// domination layers
	void rskybandLayer();

	void computeHS(Point& pt, long long int b, float* PG[], int dim);

	// rsa solution
	void rsa(vector<float>& region, unordered_map<int, cell*>& utkRet, const int k, const int dim, float* PG[], Rtree& a_rtree);

	// dirll optimization
	bool drill(Point& focal, vector<cell*>& leaves, vector<cell*>& finalResult, const int k, float* PG[]);

	// functions for jaa
	void findAnchor(vector<float>& pivot, unordered_set<long int>& ignoreset, long int& focal, unordered_set<long int>& focaled, const int dim, float* PG[], const int k);
	
	void jointArrangement(long int& epoch, cellTree* sol, const long int focal, vector<cell*>& exactutk, unordered_set<long int>& focaled,
		const int k, const int dim, unordered_set<long int>& ignoreset, float* PG[], Rtree& rtree, vector<cell*>& subregions,
		vector<unordered_map<long int, bool>>& intersection, vector<int>& protectedAttribute, vector<vector<int>>& generalization,
		vector<vector<vector<float> > >& partitionPoint, float & maxSP, float & minLoss ,float & minDistance,
		float alpha, float beta, int numOfKind, float SPinTheory, vector<float> originalF, vector<int> kindCount, int num_levels, float* bestF, int* bestTopK, int algorithm, bool* relaxed_found);

	void findPivot(cell* subreg, vector<float>& pivot, vector<float> originalF);

	void initHalfspace(cell* subreg);

	void insertingRecords(cellTree* candidate, vector<long int>& processRecords, const long int focal, float* PG[], Point& pt, long int& epoch, const int dim);

	void InsertRecords(cellTree* subreg, vector<long long int>& processRecords, const long int focal, float* PG[], Point& pt, long int& epoch, const int dim,  vector<cell*>& subregions);

	// jaa solution
	void jaa(vector<float>& region, vector<cell*>& exactutk, int& k, const int dim, float* PG[], Rtree& a_rtree, vector<unordered_map<long int, bool> >& intersection,
		vector<int>& protectedAttribute, vector<vector<int>>& generalization, vector<vector<vector<float> > >& partitionPoint, float SPinTheory,
		float alpha, float beta, int numOfKind, vector<float> originalF, float & maxSP, float & minDistance, float & minLoss,
		vector<int> kindCount, int num_levels, float* bestF, int* bestTopK, int algorithm, bool* relaxed_found);

	double treeSize;

private:
	cell* root;
	unordered_map<long int, vector<long int>> layers;
	unordered_map<long int, unordered_set<long int>> daGraph;
	unordered_map<long int, unordered_set<long int>> dominationship;
};
#endif

