#include "utk.h"
#include "skyline.h"

extern "C"{
#include "/usr/local/GNU/glpk/include/glpk.h"
}

UTK::UTK()
{
	treeSize = 0.0;
}

UTK::~UTK()
{
	layers.clear();
	daGraph.clear();
	dominationship.clear();
}

void UTK::pivotRegion(vector<float>& R, vector<float>& pivot)
{
	int dim = R.size() / 2;
	pivot.clear(); 
	for (int i = 0; i < dim; i++)
		pivot.push_back((R[2 * i] + R[2 * i + 1]) / 2.0);
}
       
bool UTK::isRdominated(const int  dim, vector<float>& R, float focal[], vector<float>& entry, bool& fDe)
{
	// generate HP;
	vector<float> tmpHS;
	float entry_d = entry[dim - 1];
	float focal_d = focal[dim - 1];
	for (int d = 0; d < dim - 1; d++)
	{
		tmpHS.push_back((focal[d] - focal_d) - (entry[d] - entry_d));
	}
	tmpHS.push_back(entry_d - focal_d);


	// verify each vertic
	int posCount = 0;
	int negCount = 0;
	int totVertics = pow(2, dim - 1);
	for (int i = 0; i < totVertics; i++)
	{
		stringstream ss;
		string tmpS = bitset<MAXDIMEN>(i).to_string();
		tmpS = tmpS.substr(tmpS.size() - (dim - 1), dim - 1);

		float sum = 0;
		for (int si = 0; si < tmpS.size(); si++)
		{
			if (tmpS[si] == '0')
				sum += R[si * 2] * tmpHS[si];
			else if (tmpS[si] == '1')
				sum += R[si * 2 + 1] * tmpHS[si];
			else
				cout << "bug here!!!" << endl;
		}

		if (sum > tmpHS[dim - 1])
			negCount++;
		else
			posCount++;
	}

	// dominationship
	if (negCount == totVertics)
	{
		fDe = true;
		return true;
	}
	else if (posCount == totVertics)
	{
		fDe = false;
		return true;
	}
	else
		return false;
}

bool UTK::countRegionDominator(int dimen, float pt[], vector<long int>& rskyband, float* PG[], vector<float>& R, const int k, unordered_set<long int>& dominators)
{
	vector<float> record(dimen,0);
	dominators.clear();
	bool fDe;
	int count = 0; 

	for (int i = 0; i < rskyband.size(); i++)
	{
		for (int di = 0; di < dimen; di++)
			record[di] = (PG[rskyband[i]][di] + PG[rskyband[i]][di + dimen]) / 2.0;
		int flagSame = 1;
		for(int di = 0; di < dimen; di++)
		{
			if(abs(record[di] - pt[di]) >= 0.0001)
			{
				flagSame = 0;
				break;
			}
		}
		if(flagSame) continue;
		fDe = false;

		for (int dii = 0; dii < dimen; dii++)
		{
			//if (record[dii] + SIDELEN < pt[dii])
			if (record[dii]  < pt[dii])
			{
				fDe = true;
				break;
			}
		}
		if (fDe == false)
		{
			count++;
			//cout<<rskyband[i]<<endl;
			dominators.insert(rskyband[i]);
		}

	}
	if (count < k)
	{
		return true;
	}
	else
	{
		return false;
	}
}

// new for FairTQ
bool judge(const pair<float, long int> a, const pair<float, long int> b) {
    return a.first>b.first;
}
// new for FairTQ
bool judgeDouble(const pair<double, long int> a, const pair<double, long int> b) {
    return a.first>b.first;
}

// new for FairTQ
vector<long int> UTK::TopKQuery(vector<float>& pivot, int dim, float* PG[], int k)
{
	vector<pair<float, long int>> anchors;

	float U_RANGE = 1000.0;
	for (auto iter = daGraph.begin() ; iter != daGraph.end() ; iter++)
	{
		long int index = iter->first;
		float ret = 0;
		float weight = 0;

		for (int di = 0; di < dim ; di++)
		{
			weight += pivot[di];
			ret += pivot[di] * (PG[index][di] + PG[index][di + dim]) / 2.0;
		}
		ret += (U_RANGE - weight)*(PG[index][dim - 1] + PG[index][dim - 1 + dim]) / 2.0;

		anchors.push_back(make_pair(ret, index));
	}

	sort(anchors.begin(), anchors.end(), judge);
	vector<long int> TopKSet;

	for(int i = 0; i< k; i++)
	{
		TopKSet.push_back(anchors[i].second);
	}

	return TopKSet;
}

// new for FairTQ
vector<long int> UTK::TopKQuery(vector<double>& pivot, int dim, float* PG[], int k)
{
	vector<pair<double, long int>> anchors;

	float U_RANGE = 1000.0;
	for (auto iter = daGraph.begin() ; iter != daGraph.end() ; iter++)
	{
		long int index = iter->first;
		double ret = 0;
		double weight = 0;

		for (int di = 0; di < dim ; di++)
		{
			weight += pivot[di];
			ret += pivot[di] * (PG[index][di] + PG[index][di + dim]) / 2.0;
		}
		ret += (U_RANGE - weight)*(PG[index][dim - 1] + PG[index][dim - 1 + dim]) / 2.0;

		anchors.push_back(make_pair(ret, index));
	}

	sort(anchors.begin(), anchors.end(), judgeDouble);
	vector<long int> TopKSet;

	for(int i = 0; i< k; i++)
	{
		TopKSet.push_back(anchors[i].second);
	}

	return TopKSet;
}

float UTK::orderScore(vector<float>& pivot, float entry[], int dim)
{
	float ret=0;
	float weight = 0;
	float U_RANGE = 1000.0;
	for (int i = 0; i < dim - 1; i++)
	{
		weight += pivot[i];
		ret += pivot[i] * entry[i];
	}
	ret += (U_RANGE - weight)*entry[dim-1];
	return ret;
}

void UTK::maintainGraph(const int recordID, float* PG[], vector<float>& region, int dim)
{
	long int tmpID = recordID;
	vector<float> record(dim, 0);
	float* focal = new float[dim];
	bool fDe = true;

	for (int i = 0; i < dim; i++)
		focal[i] = (PG[recordID][i] + PG[recordID][i + dim]) / 2.0;

	for (auto iter = daGraph.begin(); iter != daGraph.end(); iter++)
	{
		for (int i = 0; i < dim; i++)
			record[i] = (PG[iter->first][i] + PG[iter->first][i + dim]) / 2.0;
		if (isRdominated(dim, region, focal, record, fDe) && fDe==false)
		{
			daGraph[recordID].insert(iter->first);
		}
	}
	unordered_set<long int> newSet;
	daGraph[tmpID] = newSet;
}

void UTK::rskyband(vector<float>& region, const int dimen, Rtree& a_rtree, vector<long int>& rskyband, float* PG[], const int k)
{
	vector<float> pivot;
	pivotRegion(region, pivot);
	for (int j =0 ; j <dimen-1; j++)
		cout<<"piovt: "<<pivot[j]<<endl;
	

	float ORIGNIN[MAXDIMEN];
	float mindist;
	for (int i = 0; i < dimen; i++)
		ORIGNIN[i] = 1;


	RtreeNode* node;
	priority_queue<pair<float, int>> heap;
	int NegPageid;

	float pt[MAXDIMEN];
	float maxscore;
	int pageID;
	float tmpScore;
	unordered_set<long int> dominators;

	heap.push(make_pair(INFINITY, a_rtree.m_memory.m_rootPageID));
	int count = 0;
	while (!heap.empty())
	{
		tmpScore = heap.top().first;
		pageID = heap.top().second;
		heap.pop();

		if (pageID > MAXPAGEID)      ///>= mine
		{
			//cout<<"pageID:"<<pageID<<endl;	
			for (int j = 0; j < dimen; j++)
			{
				
				pt[j] = (PG[pageID - MAXPAGEID][j] + PG[pageID - MAXPAGEID][j + dimen]) / 2.0;
				//cout<<pt[j];	
			}
		
			if (countRegionDominator(dimen, pt, rskyband, PG, region, k, dominators)) 
			{
				rskyband.push_back(pageID - MAXPAGEID);
				daGraph[pageID - MAXPAGEID] = dominators;
			}
		}
		else
		{
			node = a_rtree.m_memory.loadPage(pageID);
			//cout<<"iD"<<pageID<<endl;
			if (node->isLeaf())
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int j = 0; j < dimen; j++)
					{
						pt[j] = node->m_entry[i]->m_hc.getLower()[j] + SIDELEN;
						//cout<< pt[j]<<" ";
					}					//cout<<endl;				
					if (countRegionDominator(dimen, pt, rskyband, PG, region, k, dominators))
					{
						maxscore = orderScore(pivot, pt, dimen);
						//mindist = minDist(pt, ORIGNIN, dimen);
						//cout<< "maxscore:"<<maxscore<<endl;
						heap.push(make_pair(maxscore, node->m_entry[i]->m_id + MAXPAGEID));
						//cout<< node->m_entry[i]->m_id<<endl;
					}
				}
			}
			else
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int j = 0; j < dimen; j++)
					{
						pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
					}
					if (countRegionDominator(dimen, pt, rskyband, PG, region, k, dominators))
					{
						maxscore = orderScore(pivot, pt, dimen);
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.push(make_pair(maxscore, node->m_entry[i]->m_id));
						//cout<<"idhere:"<< node->m_entry[i]->m_id<<endl;
					}
				}
			}
		}
	}
}

// new for FairTQ
void UTK::onionToRskyband(vector<float>& region, const int dimen, Rtree& a_rtree, vector<long int>& rskyband, float* PG[], const int k, vector<long int>klayers )
{
	vector<float> pivot;
	pivotRegion(region, pivot);

	daGraph.clear();
	rskyband.assign(klayers.begin(), klayers.end());

	unordered_set<long int> dominators;
	float pt[MAXDIMEN];

	for(int i =0 ; i < rskyband.size(); i++)
	{
		for (int j = 0; j < dimen; j++)
		{
			pt[j] = (PG[rskyband[i]][j] + PG[rskyband[i]][j + dimen]) / 2.0;	
		}

		if (countRegionDominator(dimen, pt, rskyband, PG, region, k, dominators)) 
		{
			daGraph[rskyband[i]] = dominators;
		}

	}
	rskyband.clear();
	int count = 0;
	for(auto iter = daGraph.begin(); iter!= daGraph.end(); iter++)
	{
		rskyband.push_back(iter->first);
		count++;
	}

	// cout<<"count in rskyband: "<<count<<endl;
}

void UTK::rskybandLayer()
{
	unordered_set<long int> initSet;
	//cout<<"layers: "<<layers.size()<<endl;
	for (auto iter = layers.begin(); iter != layers.end(); iter++)
	{
		iter->second.clear();
	}
	//cout<<"layers: "<<layers.size()<<endl;
	for (auto iter = daGraph.begin(); iter != daGraph.end(); iter++)
	{
		layers[iter->second.size()].push_back(iter->first);
		if (dominationship.find(iter->first) == dominationship.end())
		{
			dominationship[iter->first] = initSet;
		}
		for (auto second = iter->second.begin(); second != iter->second.end(); second++)
		{
			dominationship[*second].insert(iter->first);
		}
	}
}

bool UTK::drill(Point& focal, vector<cell*>& leaves, vector<cell*>& finalResult, const int k, float* PG[])
{
	vector<float> pivot;
	int dim = focal.m_dimen;
	multimap<float, long int> anchors;
	float ret = 0;
	float weight = 0;

	for (int i = 0; i < leaves.size(); i++)
	{
		//findPivot(leaves[i], pivot); not use now
		
		for (int i = 0; i < layers.size(); i++)
		{
			for (auto ai = layers[i].begin(); ai != layers[i].end(); ai++)
			{
				ret = 0;
				weight = 0;
				for (int di = 0; di < dim - 1; di++)
				{
					weight += pivot[di];
					ret += pivot[di] * (PG[*ai][di] + PG[*ai][di + dim]) / 2;
				}
				ret += (1 - weight)*(PG[*ai][dim - 1] + PG[*ai][dim - 1 + dim]) / 2;
				
				if (anchors.size() < k)
				{
					anchors.insert(make_pair(ret, *ai));
				}
				else
				{
					if (ret>anchors.begin()->first)
					{
						anchors.erase(anchors.begin());
						anchors.insert(make_pair(ret, *ai));
					}
				}
			}
		}
		
		ret = 0;
		weight = 0;
		for (int di = 0; di < dim - 1; di++)
		{
			weight += pivot[di];
			ret += pivot[di] * focal.m_coor[di];
		}
		ret += (1 - weight)*focal.m_coor[dim-1];

		if (ret>anchors.begin()->first)
		{
			finalResult.push_back(leaves[i]);
			return true;
		}
	}
	return false;
}

void UTK::computeHS(Point& pt, long long int b, float* PG[], int dim)
{
	if (RecordIDtoHalfPlaneID.find(b) == RecordIDtoHalfPlaneID.end())
	{
		vector<float> entry(dim, 0);
		// cout<<"b\%objCnt: "<<b%objCnt<<endl;
        long long int id=b%objCnt;
        if(id==0){
            id=objCnt;
        }
		float U_RANGE = 1000.0;
		for (int i = 0; i < dim; i++)
		{
			entry[i] = (PG[id][i] + PG[id][i + dim]) / 2.0;
		}

		vector<float> tmpHS;
		float entry_d = entry[dim - 1];
		float focal_d = pt.m_coor[dim - 1];
		for (int d = 0; d < dim - 1; d++)
		{
			tmpHS.push_back((entry[d] - entry_d) - (pt.m_coor[d] - focal_d));
		}
		tmpHS.push_back( U_RANGE * (focal_d - entry_d));
		tmpHS.push_back(b);         
		long int size1=  HalfSpaces.size(); 
		HalfSpaces.push_back(tmpHS); 
		long int size2=  HalfSpaces.size();
		if(size2 - size1 != 1) cout<<"out of size: "<<HalfSpaces.size()<<endl;
		RecordIDtoHalfPlaneID.insert(make_pair(b, HalfSpaces.size()));
		//cout<<"test: "<< RecordIDtoHalfPlaneID[10000000];
	}
}

void UTK::rsa(vector<float>& region, unordered_map<int, cell*>& utkRet, const int k, const int dim, float* PG[], Rtree& rtree)
{
	vector<pair<long int, unordered_set<long int>>> orderedgNode;
	for (auto iter = daGraph.begin(); iter != daGraph.end(); iter++)
	{
		orderedgNode.push_back(make_pair(iter->first, iter->second));
	}
	sort(orderedgNode.begin(), orderedgNode.end(), daGraphCompare());
	rskybandLayer();

	for (int nodei = 0; nodei < orderedgNode.size() && !daGraph.empty(); nodei++)
	{

		cout << "============================================" << endl;
		cout << "No. " << nodei + 1 << ", Record ID: " << orderedgNode[nodei].first << endl;

		if (daGraph.find(orderedgNode[nodei].first)==daGraph.end())
		{
			continue;
		}

		vector<cell*> finalResult;
		vector<cell*> leaves;
		unordered_set<long int> removeSL;
		unordered_set<long int> singular;
		vector<long int> processRecords;

		HalfSpaces.clear();
		RecordIDtoHalfPlaneID.clear();

		// init focal record
		pair<long int, unordered_set<long int>> iter = orderedgNode[nodei];
		Point pt;
		pt.m_dimen = dim;
		for (int j = 0; j < dim; j++)
		{
			pt.m_coor[j] = (PG[iter.first][j] + PG[iter.first][j + dim]) / 2.0;
		}
		int updatek = k - iter.second.size() - 1;

		initHS(dim, region);
		cellTree* sol = new cellTree(2 * (dim - 1));

		for (int ki = 0; ki < k; ki++)
		{	
			if (ki == 0)
			{
				for (int i = 0; i < layers[ki].size(); i++)
				{
					if (iter.first != layers[ki][i] && iter.second.find(layers[ki][i]) == iter.second.end())
					{
						processRecords.push_back(layers[0][i]);
						computeHS(pt, layers[0][i], PG, dim);
					}
				}
				cout << "(1) Skylines(" << k << "): " << processRecords.size() << endl;
				//sol->insert(processRecords, updatek, rtree, pt, finalResult); no use now -- zz
			}
			else
			{
				sol->opt_insert(processRecords, updatek, rtree, pt, finalResult);
			}
			sol->collectLeaf(leaves, updatek);
			if (treeSize < sol->treeSize)
			{
				treeSize = sol->treeSize;
			}

			if (leaves.empty())
			{
				cout << "(3) Singular, #Ret: " << finalResult.size() << endl;
				break;
			}

			sol->markSingular(leaves, removeSL, updatek, finalResult, singular);
			if (finalResult.size() != 0)
			{
				cout << "(3) Singular, #Ret: " << finalResult.size() << endl;
				long int utkID = iter.first;
				utkRet[utkID] = finalResult[0]; 
				for (auto idIter = daGraph[utkID].begin(); idIter != daGraph[utkID].end(); idIter++)
				{
					if (daGraph.find(*idIter) != daGraph.end())
					{
						cout << "remove: " << *idIter << endl;
						utkRet[*idIter] = finalResult[0]; // need update
						daGraph.erase(*idIter);
					}
				}
				daGraph.erase(utkID);
				break;
			}
			else if (ki == k)
			{
				cout << "please check why this happens" << endl;
				break;
			}
			else // obtain next batch
			{
					processRecords.clear();
					for (int i = 0; i < layers[ki + 1].size(); i++)
					{
						bool isNext = false;
						for (auto newIter = removeSL.begin(); newIter != removeSL.end(); newIter++)
						{
							if (daGraph[layers[ki + 1][i]].find(*newIter) != daGraph[layers[ki + 1][i]].end())
							{
								isNext = true;
								break;
							}
						}
						if (isNext)
						{
							processRecords.push_back(layers[ki + 1][i]);
							computeHS(pt, layers[ki + 1][i], PG, dim);
						}
					}
				//}
			}
		}
		sol->~cellTree();
	}
}

void UTK::findAnchor(vector<float>& pivot, unordered_set<long int>& ignoreset, long int& focal, unordered_set<long int>& focaled, const int dim, float* PG[], const int k)
{
	multimap<float, long int> anchors;
	float U_RANGE = 1000.0;
	for (auto iter = daGraph.begin(); iter != daGraph.end(); iter++)
	{
		if (ignoreset.find(iter->first) == ignoreset.end() && focaled.find(iter->first) == focaled.end()) // ignore these dominators
		{
			float ret = 0;
			float weight = 0;

			for (int di = 0; di < dim - 1; di++)
			{
				weight += pivot[di];
				ret += pivot[di] * (PG[iter->first][di] + PG[iter->first][di + dim]) / 2.0;
			}
			ret += (U_RANGE - weight)*(PG[iter->first][dim - 1] + PG[iter->first][dim - 1 + dim]) / 2.0;

			if (anchors.size() < k)
			{
				anchors.insert(make_pair(ret, iter->first));
			}
			else
			{
				if (ret>anchors.begin()->first)
				{
					anchors.erase(anchors.begin());
					anchors.insert(make_pair(ret, iter->first));
				}
			}
		}
	}
	
	focal = anchors.begin()->second;
}

void UTK::findPivot(cell* subreg, vector<float>& pivot, vector<float> originalF)
{
    if (subreg->isPruned == false) 
    {
        for (auto iter = subreg->IntersectHP.begin(); iter != subreg->IntersectHP.end(); iter++)
        {
        	long long int iter_first_id = iter->first%objCnt;
        	if (iter_first_id == 0) {
        		iter_first_id = objCnt;
        	} // may not needed since it comes from the original UTK code?
        	// cout << "In findPivot, iter->first = " << iter->first << ", iter_first_id = " << iter_first_id << endl;
            if (iter->second == true &&  subreg->TopKSet.find(iter_first_id)!=subreg->TopKSet.end())
            {
                if (subreg->ignoreset.find(iter_first_id) == subreg->ignoreset.end())
                {
                    subreg->ignoreset.insert(iter_first_id);
                }
            }
        }
        for (auto iter = subreg->BelowHP.begin(); iter != subreg->BelowHP.end(); iter++)
        {
        	long long int iter_id = *iter%objCnt;
        	if (iter_id == 0) {
        		iter_id = objCnt;
        	} // may not needed since it comes from the original UTK code?
        	// cout << "In findPivot, *iter = " << (*iter) << ", iter_id = " << iter_id << endl;
            if (subreg->ignoreset.find(iter_id) == subreg->ignoreset.end())
            {
                subreg->ignoreset.insert(iter_id);
            }
        }
    }

    int dimen = HalfSpaces[0].size() - 2;
    //newHope(subreg->IntersectHP);
    pivot = QuadProgFULL(subreg->IntersectHP, originalF, dimen+1 );

}

// seems not used
void UTK::insertingRecords(cellTree* subreg, vector<long int>& processRecords, const long int focal, float* PG[], Point& pt, long int& epoch, const int dim)
{
	long long int recordID;
	processRecords.clear();
	unordered_set<long int> invalidRecords;
	invalidRecords.insert(focal);
	
	for (auto iter = subreg->root->focaled.begin(); iter != subreg->root->focaled.end(); iter++) // ignore focaled records
	{
		invalidRecords.insert(*iter);
	}
	for (auto iiter = dominationship[focal].begin(); iiter != dominationship[focal].end(); iiter++) // Lemma 1: ignore records which are r-dominated by focal record
	{
		invalidRecords.insert(*iiter);
	}
	for (auto iter = subreg->root->ignoreset.begin(); iter != subreg->root->ignoreset.end(); iter++)
	{
		invalidRecords.insert(*iter);
		for (auto iiter = daGraph[*iter].begin(); iiter != daGraph[*iter].end(); iiter++)
		{
			invalidRecords.insert(*iiter);
		}
	}

	for (auto iter = daGraph.begin(); iter != daGraph.end(); iter++)
	{
		if (invalidRecords.find(iter->first) == invalidRecords.end())
		{
			recordID = epoch*objCnt + iter->first;
			processRecords.push_back(recordID);
			computeHS(pt, recordID, PG, dim);
		}
	}
} 

void UTK::InsertRecords(cellTree* subreg, vector<long long int>& processRecords, const long int focal, float* PG[], Point& pt, long int& epoch, const int dim, vector<cell*>& subregions)
{
	long long int recordID;
	processRecords.clear();
	unordered_set<long int> invalidRecords;
	invalidRecords.insert(focal);
	if (subreg->root->isPruned == true)
	{
		subreg->root->rank = subreg->root->ignoreset.size();
		if(subreg->root->ignoreset.size() != subreg->root->TopKSet.size())
			cout<<"size: "<< subreg->root->ignoreset.size()<<" "<<subreg->root->TopKSet.size()<<endl;
		subreg->root->isPruned = false; 
		for (auto iter = subreg->root->focaled.begin(); iter != subreg->root->focaled.end(); iter++) // ignore focaled records
		{
			invalidRecords.insert(*iter);
			for (auto iiter = dominationship[*iter].begin(); iiter != dominationship[*iter].end(); iiter++) // Lemma 1: ignore records which are r-dominated by focal record
			{
				invalidRecords.insert(*iiter);
			}
		}
		for (auto iter = subreg->root->ignoreset.begin(); iter != subreg->root->ignoreset.end(); iter++) // skip ignored records
		{
			invalidRecords.insert(*iter);
		}
		for (auto iiter = dominationship[focal].begin(); iiter != dominationship[focal].end(); iiter++) // Lemma 1: ignore records which are r-dominated by focal record
		{
			invalidRecords.insert(*iiter);
		}

	}
	else
	{
		for (auto iter = subreg->root->focaled.begin(); iter != subreg->root->focaled.end(); iter++) // ignore focaled records
		{
			invalidRecords.insert(*iter);
			for (auto iiter = dominationship[*iter].begin(); iiter != dominationship[*iter].end(); iiter++) // Lemma 1: ignore records which are r-dominated by focal record
			{
				invalidRecords.insert(*iiter);
			}

		}
		for (auto iiter = dominationship[focal].begin(); iiter != dominationship[focal].end(); iiter++) // Lemma 1: ignore records which are r-dominated by focal record
		{
			invalidRecords.insert(*iiter);
		}
		
		for (auto iter = subreg->root->ignoreset.begin(); iter != subreg->root->ignoreset.end(); iter++)
		{
			invalidRecords.insert(*iter);
		}
	}
	for (auto iter = daGraph.begin(); iter != daGraph.end(); iter++)
	{
		if (invalidRecords.find(iter->first) == invalidRecords.end())
		{
			recordID = epoch*objCnt + iter->first;
			processRecords.push_back(recordID);
			computeHS(pt, recordID, PG, dim);
		}
	}
}

// new for FairTQ
vector<double> UTK::find_feasible(vector<vector<float> > hyperplane)
{
	int M = hyperplane.size();
	int D = hyperplane[0].size() - 1;

	// D + 2variables: D for dim, 2 for additional var for feasible
	int* ia = new int[1 + (D + 2) * M];  //TODO: delete
	int* ja = new int[1 + (D + 2) * M];  //TODO: delete
	double* ar = new double[1 + (D + 2) * M];   //TODO: delete
	int i, j;
	double epsilon = 0.0000000000001;


	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_prob_name(lp, "find_feasible");
	glp_set_obj_dir(lp, GLP_MAX);


	glp_add_rows(lp, M);  // add D rows: q_1...q_D
							  // Add rows q_1 ... q_D
	for (i = 1; i <= M; i++) {
		char buf[10];
		sprintf(buf, "q%d", i);
		glp_set_row_name(lp, i, buf);
		glp_set_row_bnds(lp, i, GLP_UP, 0, 0); // qi = 0
	}
	

	glp_add_cols(lp, D + 2);    // add D columns: v[1] ... v[D]
								// Add col v[1] ... v[D]
	for (i = 1; i <= D + 2; i++) {
		char buf[10];
		sprintf(buf, "v%d", i);

		glp_set_col_name(lp, i, buf);

		if(i <= D)
			glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0); // -infty <= v[i] < infty
		else if (i == D + 1)
			glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0); // 0 <= v[i] < infty
		else
			glp_set_col_bnds(lp, i, GLP_UP, 0.0, D+1);

		if(i == D + 2)
			glp_set_obj_coef(lp, i, 1);  // objective: 0
		else
			glp_set_obj_coef(lp, i, 0.0);  // objective: 0
	}


	int counter = 1;
	// set value on row q1 ... qD
	for (i = 1; i <= M; i++) {
		for (j = 1; j <= D + 2; j++) {

			ia[counter] = i; ja[counter] = j;
			
			if(j <= D)
			{
				ar[counter++] = hyperplane[i-1][j-1];
				//printf("%lf ", hyperplane[i-1]->normal->coord[j-1]);
			}
			else if (j == D+1)
			{
				ar[counter++] = hyperplane[i-1][j-1];
				//printf("%lf ", hyperplane[i-1]->offset);
			}
			else if (j == D+2)
			{
				ar[counter++] = 1;
				//printf("1.00000\n");
			}
		}
	}

	// loading data  
	glp_load_matrix(lp, counter - 1, ia, ja, ar);

										  // running simplex
	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_OFF; // turn off all message by glp_simplex 

	glp_simplex(lp, &parm);


	vector<double> feasible_pt(D);
	double w1, w2;
	w1 = glp_get_col_prim(lp, D+1);
	w2 = glp_get_col_prim(lp, D+1);

	if(w1 < 0 || w2 < 0 || w1 == 0.0 || w2 == 0.0)
	{
		printf("LP feasible error.\n");
		//return NULL;
	}
	for (i = 0; i < D; i++)
	{
		double v = glp_get_col_prim(lp, i + 1);
		//printf("w%d = %lf\n", i + 1, v);
		feasible_pt[i] = v / w1;
	}

	//printf("solution status: %d\n",glp_get_status(lp));
	//printf("solution is unbounded: %d\n", GLP_UNBND);
	//printf("solution is optimal: %d\n", GLP_OPT);
	//printf("no feasible solution: %d\n", GLP_NOFEAS);
	//printf("return: %lf\n", glp_get_obj_val(lp));
	//for (i = 0; i < D + 2; i++)
	//{
	//	double v = glp_get_col_prim(lp, i + 1);
	//	printf("w%d = %lf\n", i + 1, v);
	//}

	glp_delete_prob(lp); // clean up
	delete[]ia;
	delete[]ja;
	delete[]ar;

	return feasible_pt;
}

// new for FairTQ
vector<double> UTK::newHope( unordered_map<long int, bool>& utility_hyperplane)
{
	vector<vector<float> > hyperplanes;
	int dim = HalfSpaces[0].size() - 2;
	for (auto iter = utility_hyperplane.begin(); iter != utility_hyperplane.end(); iter++)
	{
		long int hsID = RecordIDtoHalfPlaneID[iter->first] - 1;

		vector<float> hyperplane;
		for (int j = 0; j < dim; j++)
		{
			if(iter->second == false)
				hyperplane.push_back(HalfSpaces[hsID][j]);
			else
				hyperplane.push_back(-HalfSpaces[hsID][j]);
		}
		if(iter->second == false)
			hyperplane.push_back(-HalfSpaces[hsID][dim]);
			//fprintf(wPtr, "%.6f ",-HalfSpaces[hsID][dim]);
		else
			hyperplane.push_back(HalfSpaces[hsID][dim]);


		hyperplanes.push_back(hyperplane);
	}
	return find_feasible(hyperplanes);
}

bool checkAlphaApp(unordered_set<long int>& TopKSet, int k, vector<int>& protectedAttribute, float alpha, int numOfKind)
{
	vector<int> countKind(numOfKind);
	long int objcnt = protectedAttribute.size();


	for(auto iter = TopKSet.begin();  iter!= TopKSet.end() ;++iter)
    {
        // countKind[protectedAttribute[*iter%objcnt]]++; // by Hao
        countKind[protectedAttribute[*iter - 1]]++;
    }

    int flag = 1;

    for(int i=0; i<numOfKind; i++)
    {
    	if(countKind[i]< alpha* float(k))
    		return 0;

    }
    return flag;
}

void UTK::jointArrangement(long int& epoch, cellTree* sol, const long int focal, vector<cell*>& exactutk, unordered_set<long int>& focaled,
	const int k, const int dim, unordered_set<long int>& ignoreset, float* PG[], Rtree& rtree, vector<cell*>& subregions, vector<unordered_map<long int, bool>>& intersection,
	vector<int>& protectedAttribute , vector<vector<int>>& generalization, vector<vector<vector<float> > >& partitionPoint, float & maxSP, float & minLoss ,float & minDistance,
	float alpha, float beta, int numOfKind, float SPinTheory, vector<float> originalF, vector<int> kindCount, int num_levels, float* bestF, int* bestTopK, int algorithm, bool* relaxed_found)
{
	// cout << "Start jointArrangement for selected focal (anchor): " << focal << endl;

	// cout << "dominationship of focal:";
	// for (auto & v: dominationship[focal]) {
	// 	cout << " " << v;
	// }
	// cout << endl;

	// cout << "Processing subreg:" << endl;
	// cout << "\tIntersectHP (" << sol->root->IntersectHP.size() << ")";
	// for (auto & v: sol->root->IntersectHP) {
	// 	cout << " (" << v.first << ", " << v.second << ")";
	// }
	// cout << endl;
	// cout << "\tTopKSet (" << sol->root->TopKSet.size() << ")";
	// for (auto & v: sol->root->TopKSet) {
	// 	cout << " " << v;
	// }
	// cout << endl;
	// cout << "\tignoreset (" << ignoreset.size() << ")";
	// for (auto & v: ignoreset) {
	// 	cout << " " << v;
	// }
	// cout << endl;
	
	float sp = 0; // unused
	float ratio = 0; // unused

	vector<cell*> finalResult;
	vector<cell*> leaves;
	unordered_set<long int> removeSL;
	unordered_set<long int> singular;
	vector<long long int> processRecords;

	// init focal record
	Point pt;
	pt.m_dimen = dim;
	for (int j = 0; j < dim; j++)
	{
		pt.m_coor[j] = (PG[focal][j] + PG[focal][j + dim]) / 2.0;
	}

	// added by Zheng
	if(sol->root->extremePoint.size() == 0) {
		// cout << "Handle extremePoint for dim = " << dim << endl;
		if(dim == 2)
		{
			for(int i=0; i < 2; i++ )
			{
				vector<float> vertices;
				int id1 = -(i+1);
				int hsID1 = RecordIDtoHalfPlaneID[id1] - 1;
				vertices ={abs(HalfSpaces[hsID1][dim-1])};

				//cout<< abs(HalfSpaces[hsID1][dim-1])<<endl;
				sol->root->extremePoint.push_back(vertices);
				vector<long int > newpair;
				newpair = {id1};
				sol->root->pairMap.push_back(newpair);
			}
		}
		if(dim == 3)
		{
			for(int i = 0 ; i< 2*dim - 2; i++)
			{
				vector<float> vertices;
				int id1 = -(i/2+1);
				int hsID1 = RecordIDtoHalfPlaneID[id1] - 1;
			
				int id2 = -(i%2+3);
				int hsID2 = RecordIDtoHalfPlaneID[id2] - 1;

				vertices = {abs(HalfSpaces[hsID1][dim-1]), abs(HalfSpaces[hsID2][dim-1])};
				// cout<<abs(HalfSpaces[hsID1][dim-1])<<" "<<abs(HalfSpaces[hsID2][dim-1])<<endl;
				sol->root->extremePoint.push_back(vertices);
				vector<long int > newpair;
				newpair = {id1, id2};
				sol->root->pairMap.push_back(newpair);
			}

		}

		if(dim == 4)
		{
		
			for(int i = 0; i < 8; i++)
			{
				vector<float> vertices;
				int id1 = -(i/4+1);
				int hsID1 = RecordIDtoHalfPlaneID[id1] - 1;
			
				int id2 = -((i%4)/2+3);
				int hsID2 = RecordIDtoHalfPlaneID[id2] - 1;

				int id3 = -(i%2+5);
				int hsID3 = RecordIDtoHalfPlaneID[id3] - 1;

				vertices = {abs(HalfSpaces[hsID1][dim-1]), abs(HalfSpaces[hsID2][dim-1]) , abs(HalfSpaces[hsID3][dim-1])};
				// cout<<abs(HalfSpaces[hsID1][dim-1])<<" "<<abs(HalfSpaces[hsID2][dim-1])<<" "<<abs(HalfSpaces[hsID3][dim-1])<<endl;
				sol->root->extremePoint.push_back(vertices);
				vector<long int > newpair;
				newpair = {id1, id2, id3};
				sol->root->pairMap.push_back(newpair);
			}
		
		}
		if(dim == 5)
		{
		
			for(int i = 0; i < 16; i++)
			{
				vector<float> vertices;
				int id1 = -(i/8+1);
				int hsID1 = RecordIDtoHalfPlaneID[id1] - 1;
			
				int id2 = -((i%8)/4+3);
				int hsID2 = RecordIDtoHalfPlaneID[id2] - 1;

				int id3 = -((i%4)/2+5);
				int hsID3 = RecordIDtoHalfPlaneID[id3] - 1;

				int id4 = -(i%2+7);
				int hsID4 = RecordIDtoHalfPlaneID[id4] - 1;

				vertices = {abs(HalfSpaces[hsID1][dim-1]), abs(HalfSpaces[hsID2][dim-1]) , abs(HalfSpaces[hsID3][dim-1]) , abs(HalfSpaces[hsID4][dim-1])};
				// cout<<abs(HalfSpaces[hsID1][dim-1])<<" "<<abs(HalfSpaces[hsID2][dim-1])<<" "<<abs(HalfSpaces[hsID3][dim-1])<<" "<<abs(HalfSpaces[hsID4][dim-1])<<endl;
				sol->root->extremePoint.push_back(vertices);
				vector<long int > newpair;
				newpair = {id1, id2, id3, id4};
				sol->root->pairMap.push_back(newpair);
			}
		
		}
		// cout << "Finish Handle extremePoint for dim = " << dim << endl;
	}

	unordered_set<long int> ancestors;

	for (auto iter = daGraph[focal].begin(); iter != daGraph[focal].end(); iter++) {
		ancestors.insert(*iter);
		if (sol->root->newAncestors.find(*iter) != sol->root->newAncestors.end()) {
			for (auto iiter = sol->root->newAncestors[*iter].begin(); iiter != sol->root->newAncestors[*iter].end(); iiter++) {
				ancestors.insert(*iiter);
			}
		}	
	}

	// cout << "ancestors before InsertRecords (" << ancestors.size() << ")";
	// for (auto & v: ancestors) {
	// 	cout << " " << v;
	// }
	// cout << endl;

	// vector<long int> descendants; // seems not used
	// added by Zheng

	InsertRecords(sol, processRecords, focal, PG, pt, epoch, dim, subregions);

	// cout << "processRecords after InsertRecords in epoch " << epoch << " (" << processRecords.size() << "):";
	// for (auto & v: processRecords) {
	// 	long long int v_id = v % objCnt;
	// 	if (v_id == 0) {
	// 		v_id = objCnt;
	// 	}
	// 	cout << " " << v << "(" << v_id << ")";
	// }
	// cout << endl;

	// cout << "In jointArrangement, print r-dominance table before insert" << endl;
	// for (auto & v: daGraph) {
	// 	cout << v.first << ":";
	// 	for (auto & vv: v.second) {
	// 		cout << " " << vv;
	// 	}
	// 	cout << endl;
	// } // this never changed
	
	if (processRecords.size() == 0) {
		// cout<< "Case 3: rank "<<sol->root->rank<<" Ignoreset: "<<sol->root->ignoreset.size()<<endl;
	} else {
		//cout << "(1) Insert Records: " << processRecords.size() << endl;
		
		// cout << "call sol->insert with mink = " << k << ", ancestors:";
		// // for (auto &v: ancestors) {
		// // 	cout << " " << v;
		// // }
		// cout << ", focal = " << focal << endl;

		// sol->insert(processRecords, k, rtree, pt, finalResult, ancestors,descendants, subregions, focal, k , maxSP, protectedAttribute, numOfKind);
		sol->insert(processRecords, k, rtree, pt, finalResult, ancestors, subregions, focal, k , maxSP, protectedAttribute, numOfKind);
		sol->collectLeaf(leaves, INFINITY);

		if (treeSize < sol->treeSize) {
			treeSize = sol->treeSize;
		}

		// cout << "Finally, after hyperplanes insertion (" << leaves.size() << "):" << endl;

		long int objCnt = protectedAttribute.size();
		//cout<<"leaves: "<<leaves.size()<<endl;

		for (int ci = 0; ci < leaves.size(); ci++) {
			long long int focalID = focal%objCnt;
			// cout << "In jointArrangement, focalID: " << focalID << endl;
			if (focalID == 0) {
				focalID = objCnt;
			}

			if (leaves[ci]->TopKSet.find(focalID) == leaves[ci]->TopKSet.end()) {
				leaves[ci]->TopKSet.insert(focalID);
				leaves[ci]->rank++;
			}

			// cout << "\tleaves[" << ci << "]->IntersectHP (" << leaves[ci]->IntersectHP.size() << ")";
			// for (auto & v: leaves[ci]->IntersectHP) {
			// 	cout << " (" << v.first << ", " << v.second << ")";
			// }
			// cout << endl;

			// cout << "\tleaves[" << ci << "]->TopKSet (" << leaves[ci]->TopKSet.size() << ")";
			// for (auto & v: leaves[ci]->TopKSet) {
			// 	cout << " " << v;
			// }
			// cout << endl;

			// cout << "\tleaves[" << ci << "]->newAncestors:" << endl;
			// for (auto & v: leaves[ci]->newAncestors) {
			// 	cout << "\t\t" << v.first << ":";
			// 	for (auto & vv: v.second) {
			// 		cout << " " << vv;
			// 	}
			// 	cout << endl;
			// }

			// cout << "\tleaves[" << ci << "]->newdescendants:" << endl;
			// for (auto & v: leaves[ci]->newdescendants) {
			// 	cout << "\t\t" << v.first << ":";
			// 	for (auto & vv: v.second) {
			// 		cout << " " << vv;
			// 	}
			// 	cout << endl;
			// }

			// // cout << "\tleaves[" << ci << "]->combAncestors:" << endl;
			// // for (auto & v: leaves[ci]->combAncestors) {
			// // 	cout << "\t\t" << v.first << ":";
			// // 	for (auto & vv: v.second) {
			// // 		cout << " " << vv;
			// // 	}
			// // 	cout << endl;
			// // }

			// // cout << "\tleaves[" << ci << "]->combDescendants:" << endl;
			// // for (auto & v: leaves[ci]->combDescendants) {
			// // 	cout << "\t\t" << v.first << ":";
			// // 	for (auto & vv: v.second) {
			// // 		cout << " " << vv;
			// // 	}
			// // 	cout << endl;
			// // }

			if (leaves[ci]->TopKSet.size() == k) {
				// cout << "\t-> Equal-to case" << endl;

				leaves[ci]->ignoreset.insert(focal);

				if (algorithm == 1) {
					// cout << "\t\tcheck approx1 scenario" << endl;

					if (checkAlphaApp(leaves[ci]->TopKSet, k, protectedAttribute, beta, numOfKind)) {

						vector<float> pivot_approx1 = QuadProgFULL(leaves[ci]->IntersectHP, originalF, dim);
						float distance_approx1 = Dist(pivot_approx1, originalF, dim);

						if (!(*relaxed_found)) {
							*relaxed_found = true;
							// cout << "\t\tFound a beta-relaxed-fair for the first time for beta = " << beta << ", we thus only focus on beta-relaxed-fair" << endl;

							minDistance = distance_approx1;
							for (int di = 0; di < dim; di++) {
								bestF[di] = pivot_approx1[di];
							}
							int ki = 0;
							for (auto & v: leaves[ci]->TopKSet) {
								bestTopK[ki] = v;
								ki++;
								if (ki >= k) {
									break;
								}
							}
							// cout << "\t\tDirectly update minDistance for the first time:" << endl;
							// cout<<"\t\t\tminDistance: "<<minDistance<<endl;
							// cout<<"\t\t\tpivot: " <<pivot_approx1[0]<<" "<<pivot_approx1[1]<<endl;
						} else {
							// cout << "\t\tFound a beta-relaxed-fair again for beta = " << beta << endl;

							if (minDistance > distance_approx1) {
								minDistance = distance_approx1;
								for (int di = 0; di < dim; di++) {
									bestF[di] = pivot_approx1[di];
								}
								int ki = 0;
								for (auto & v: leaves[ci]->TopKSet) {
									bestTopK[ki] = v;
									ki++;
									if (ki >= k) {
										break;
									}
								}

								// cout << "\t\tSmaller beta-relaxed-fair distance found:" << endl;
								// cout<<"\t\t\tminDistance: "<<minDistance<<endl;
								// cout<<"\t\t\tpivot: " <<pivot_approx1[0]<<" "<<pivot_approx1[1]<<endl;
							}
						}
					}

					if (*relaxed_found) {
						continue; // if we have switched to relaxed found mode, either from this check or previous check, we skipped the original exact steps
					}
				}

				// cout << "\t\tcomputeSP for leaves[" << ci << "]->TopKSet:";
				// for (auto & v: leaves[ci]->TopKSet) {
				// 	cout << " " << v << " (" << protectedAttribute[v - 1] << ")";
				// }
				// cout << endl;

				float SP = getAlphaFairness(leaves[ci]->TopKSet, PG, protectedAttribute, dim, generalization, partitionPoint, alpha, numOfKind, kindCount, num_levels);
				// cout << "SP: " << SP << endl;

				if (SP > maxSP || abs(SP - maxSP) <= 0.0001 * maxSP) { // added by Hao, since SP could be very small, cannot set this fixed value
					if(abs(SP - maxSP) <= 0.0001 * maxSP && minDistance <= 0.0001) { // SP = maxSP
						// cout << "\t\tFind equal SP but minDistance is already 0, thus skip" << endl;
						continue; // why do we skip distance = 0?
					}

					vector<float> pivot = QuadProgFULL(leaves[ci]->IntersectHP, originalF, dim);
					float distance = Dist(pivot, originalF, dim); // moving this inside this if, since we do not need to compute distance if SP is smaller

					// cout << "\t\tSP = " << SP << ", maxSP = " << maxSP << ", SP-theory = " << SPinTheory << ", distance = " << distance << endl;

					if (SP > maxSP) {
						maxSP = SP;
						minDistance = distance;

						for (int di = 0; di < dim; di++) {
							bestF[di] = pivot[di];
						}
						int ki = 0;
						for (auto & v: leaves[ci]->TopKSet) {
							bestTopK[ki] = v;
							ki++;
							if (ki >= k) {
								break;
							}
						}

						// cout << "\t\tFind larger SP, update:" << endl;
						// cout<<"\t\t\tminDistance: "<<minDistance<<endl;
						// cout<<"\t\t\tpivot: " <<pivot[0]<<" "<<pivot[1]<<endl;
						// cout<<"\t\t\tmaxSP: "<<maxSP<< " (SP-theory: "<<SPinTheory << ")"<< endl;
					} else {
						if (minDistance > distance) {
							minDistance = distance;

							for (int di = 0; di < dim; di++) {
								bestF[di] = pivot[di];
							}
							int ki = 0;
							for (auto & v: leaves[ci]->TopKSet) {
								bestTopK[ki] = v;
								ki++;
								if (ki >= k) {
									break;
								}
							}

							// cout << "\t\tFind equal SP but smaller distance, update:" << endl;
							// cout<<"\t\t\tminDistance: "<<minDistance<<endl;
							// cout<<"\t\t\tpivot: " <<pivot[0]<<" "<<pivot[1]<<endl;
							// cout<<"\t\t\tmaxSP: "<<maxSP<< " (SP-theory: "<<SPinTheory << ")"<< endl;
						}
					}
				}
			}
			else if (leaves[ci]->TopKSet.size() < k) {
				// cout << "\t-> Less-than case" << endl;

				for (auto iter = leaves[ci]->TopKSet.begin(); iter != leaves[ci]->TopKSet.end(); iter++)
					leaves[ci]->ignoreset.insert(*iter);

				//newHope(leaves[ci]->IntersectHP);
				
				vector<float> pivot = QuadProgFULL(leaves[ci]->IntersectHP, originalF, dim );

				float distance = Dist(pivot, originalF, dim);
				leaves[ci]->Distance = distance;

				if (distance > minDistance && abs(maxSP - SPinTheory) / SPinTheory < 0.0001 && maxSP > 0) { // added by Hao, since SP could be very small, cannot set this fixed value
					// cout << "\t\tDisgard this leaf since distance > minDistance && maxSP = " << maxSP << ", SPinTheory = " << SPinTheory << endl;
					sol->releaseCell(leaves[ci]);
				// } else if (distance > minDistance && abs(maxSP - SPinTheory) / SPinTheory < 0.01 && maxSP > 0) {

				} else {
					if (checkLemma3(leaves[ci]->TopKSet, k, maxSP, protectedAttribute, numOfKind)) {
						// cout << "\t\tcheckLemma3 complete, add this leaf into subregions" << endl;
						subregions.push_back(leaves[ci]);
					}
					else {
						// cout << "\t\tcheckLemma3 failed, prune this leaf" << endl;
					}
				}
				
				//subregions.push_back(leaves[ci]);
			} else {
				// cout << "\t-> Greater-than case" << endl;

				vector<float> pivot = QuadProgFULL(leaves[ci]->IntersectHP, originalF, dim );

				float distance = Dist(pivot, originalF, dim);
				leaves[ci]->Distance = distance;
				leaves[ci]->focaled.insert(focal);
			
    			
				if (distance > minDistance && abs(maxSP - SPinTheory) / SPinTheory < 0.0001 && maxSP > 0) { // added by Hao, since SP could be very small, cannot set this fixed value
					// cout << "\t\tDisgard this leaf since distance > minDistance && maxSP = " << maxSP << ", SPinTheory = " << SPinTheory << endl;
					sol->releaseCell(leaves[ci]);
				} else {
					leaves[ci]->BelowHP.clear();
					leaves[ci]->AboveHP.clear();
					leaves[ci]->TopKSet.clear();
					// cout << "\t\tAdd previous ignoreset to leaves[ci]->TopKSet:";
					for (auto iter = leaves[ci]->ignoreset.begin(); iter != leaves[ci]->ignoreset.end(); ++iter) {
						// cout << " " << (*iter);
						leaves[ci]->TopKSet.insert(*iter);
					}
					// cout << endl;

					// cout << "\t\tAdd this leaf into subregions" << endl;
					subregions.push_back(leaves[ci]);
				}
			}
		}
	}
	sol->~cellTree();

	// cout << "End jointArrangement" << endl;
}

void UTK::initHalfspace(cell* subreg)
{

}

void UTK::jaa(vector<float>& region, vector<cell*>& exactutk, int& k, const int dim, float* PG[], Rtree& a_rtree,
	vector<unordered_map<long int, bool>>& intersection, vector<int>&protectedAttribute,  vector<vector<int>>& generalization,
	vector<vector<vector<float> > >& partitionPoint, float SPinTheory, float alpha, float beta, int numOfKind, vector<float> originalF,
	float & maxSP, float & minDistance, float & minLoss, vector<int> kindCount, int num_levels, float* bestF, int* bestTopK, int algorithm, bool* relaxed_found)
{
	// cout << "Region input:";
	// for (auto & v: region) {
	// 	cout << " " << v;
	// }
	// cout << endl;

	rskybandLayer(); // update layer

	// obtain pivot
	vector<float> pivot;
	pivotRegion(region, pivot);
	// cout << "First return of pivotRegion:";
	// for (auto & v: pivot) {
	// 	cout << " " << v;
	// }
	// cout << endl;

	// obtain anchor
	unordered_set<long int> ignoreset;
	unordered_set<long int> focaled;
	long int focal = 0;
	findAnchor(pivot, ignoreset, focal, focaled, dim, PG, k);
	// cout << "First return of findAnchor: " << focal << endl;

	// cout << "ignoreset:";
	// for (auto & v: ignoreset) {
	// 	cout << " " << v;
	// }
	// cout << endl;
	// cout << "focaled:";
	// for (auto & v: focaled) {
	// 	cout << " " << v;
	// }
	// cout << endl;

	// process focal
	vector<cell*> subregions;

	HalfSpaces.clear();
	RecordIDtoHalfPlaneID.clear();
	initHS(dim, region);
	cellTree* candidate = new cellTree(2 * (dim - 1));
	long int epoch = 0;

	// cout << "Start of epoch: " << epoch << endl;
	// cout << "Initial HalfSpaces (" << HalfSpaces.size() << "):" << endl;
	// for (auto & hs: HalfSpaces) {
	// 	for (int hsi = 0; hsi < hs.size(); hsi++) {
	// 		cout << hs[hsi] << " ";
	// 		if (hsi == 0) {
	// 			cout << "x + ";
	// 		} else if (hsi == 1) {
	// 			cout << "y = ";
	// 		} else if (hsi == 2) {
	// 			cout << ", ";
	// 		}
	// 	}
	// 	cout << endl;
	// }

	// for (auto & i_daGraph: daGraph) {
	// 	auto & id_daGraph = i_daGraph.first;
	// 	auto & init_ancs = i_daGraph.second;
	// 	if (candidate->root->combAncestors.find(id_daGraph) == candidate->root->combAncestors.end()) {
	// 		candidate->root->combAncestors[id_daGraph] = init_ancs; // for initial add, this is always true
	// 	}
	// 	for (auto & init_anc: init_ancs) {
	// 		if (candidate->root->combDescendants.find(init_anc) == candidate->root->combDescendants.end()) {
	// 			candidate->root->combDescendants[init_anc] = {id_daGraph};
	// 		} else {
	// 			candidate->root->combDescendants[init_anc].insert(id_daGraph);
	// 		}
	// 	}
	// }

	// cout << "Initial combined anc-des for root:" << endl;
	// cout << "\tanc:" << endl;
	// for (auto & v: candidate->root->combAncestors) {
	// 	cout << "\t\t" << v.first << ":";
	// 	for (auto & vv: v.second) {
	// 		cout << " " << vv;
	// 	}
	// 	cout << endl;
	// }
	// cout << "\tdec:" << endl;
	// for (auto & v: candidate->root->combDescendants) {
	// 	cout << "\t\t" << v.first << ":";
	// 	for (auto & vv: v.second) {
	// 		cout << " " << vv;
	// 	}
	// 	cout << endl;
	// }

	jointArrangement(epoch, candidate, focal, exactutk, focaled, k, dim, ignoreset, PG, a_rtree,
		subregions, intersection,protectedAttribute, generalization, partitionPoint,
		maxSP, minLoss, minDistance, alpha, beta, numOfKind, SPinTheory, originalF, kindCount, num_levels, bestF, bestTopK, algorithm, relaxed_found);

	cell* subreg = NULL;
	vector<cellTree*> cellTreeSet;
	while (!subregions.empty()) {
		if (epoch % 500 == 1) {
			// cout << endl << "epoch = " << epoch << ", remain subregions: " << subregions.size() << endl;
			// cout << "HalfSpaces (" << HalfSpaces.size() << "):" << endl;
			
			// for (auto & hs: HalfSpaces) {
			// 	for (int hsi = 0; hsi < hs.size(); hsi++) {
			// 		cout << hs[hsi] << " ";
			// 		if (hsi == 0) {
			// 			cout << "x + ";
			// 		} else if (hsi == 1) {
			// 			cout << "y = ";
			// 		} else if (hsi == 2) {
			// 			cout << ", ";
			// 		}
			// 	}
			// 	cout << endl;
			// }
		}

		epoch++;

		sort(subregions.begin(), subregions.end(), cellCompare());

		vector<cell*> ::iterator it = subregions.begin();
		subreg = *it;
		subregions.erase(it);
		pivot.clear(); 
		long int objCnt = protectedAttribute.size(); // Whether this equals to global objCnt? Seems true
		
		if (subreg->Distance > minDistance && abs(maxSP - SPinTheory) / SPinTheory < 0.0001 && maxSP > 0) { // added by Hao, since SP could be very small, cannot set this fixed value
			// cout << "Disgard this region since distance > minDistance && maxSP = " << maxSP << ", SPinTheory = " << SPinTheory << endl;
			continue;
		}

		findPivot(subreg, pivot, originalF);
		// cout << "Return of findPivot:"; // pivot[0]<<" "<<pivot[1]<<endl;
		// for (auto & v: pivot) {
		// 	cout << " " << v;
		// }
		// cout  << endl;

		if (subreg->ignoreset.size() >=k) {
			//cout<<"RANK: "<<subreg->rank<<endl;
			// cout<< "Case 1: TopKSet: "<< subreg->TopKSet.size()<<" Ignoreset: "<<subreg->ignoreset.size()<<" "<<subreg->isPruned <<endl;

			exactutk.push_back(subreg);
		} else {
			//cout<<"subreg->ignoreset.size(): "<<subreg->ignoreset.size()<<endl;
			findAnchor(pivot, subreg->ignoreset, focal, subreg->focaled, dim, PG, k - subreg->ignoreset.size());
			// cout << "Return of findAnchor: " << focal << endl;
			
			candidate = new cellTree(subreg);
			jointArrangement(epoch, candidate, focal, exactutk, subreg->focaled, k, dim, subreg->ignoreset, PG, a_rtree,
				subregions, intersection, protectedAttribute, generalization, partitionPoint,
				maxSP, minLoss, minDistance, alpha, beta, numOfKind, SPinTheory, originalF, kindCount, num_levels, bestF, bestTopK, algorithm, relaxed_found);
		}
	}
}
