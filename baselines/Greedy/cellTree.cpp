#include "cellTree.h"
#include <iostream> 
#include <algorithm> 
#include <vector>
using namespace std;

extern vector<vector<float>> HalfSpaces;
extern unordered_map<long int, long int> RecordIDtoHalfPlaneID;
extern double totalSpaceCost;
extern unordered_map<long int, RtreeNode*> ramTree;

float min(float a, float b){
    if(a > b)
        return b;
    else
        return a;
}

float max(float a, float b){
    if(a > b)
        return a;
    else
        return b;
}

cellTree::cellTree()
{
	root = new cell();
	treeSize = 0;
}

cellTree::cellTree(cell* node)
{
	root = new cell();
	node->copyleaf(root);
	treeSize = 0;
}

cellTree::cellTree(int size)
{
	root = new cell();

	for (int i = 0; i < size; i++)
	{
		root->IntersectHP[-(i+1)] = false;
	}
	treeSize = 0;
}

cellTree::~cellTree()
{
	for (auto iter = dagNode.begin(); iter != dagNode.end(); iter++)
	{
		iter->second->rDominator.clear();
		unordered_set<long int>().swap(iter->second->rDominator);
		delete iter->second;
	}
	dagNode.clear();
	unordered_map<long int, gNode*>().swap(dagNode);
	releaseCell(root);
}

void cellTree::releaseCell(cell* node)
{
	cell* tmp;
	for (tmp = node; tmp != nullptr && tmp->left != nullptr; tmp = tmp->left);

	while (node != nullptr)
	{
		for (tmp = node; tmp != nullptr && tmp->left != nullptr; tmp = tmp->left);
		cell* old = node;
		node = node->left;
		old->IntersectHP.clear();
		unordered_map<long int, bool>().swap(old->IntersectHP);
		old->BelowHP.clear();
		unordered_set<long int>().swap(old->BelowHP);
		old->AboveHP.clear();
		unordered_set<long int>().swap(old->AboveHP);
		///////////////////////////////////////////
		old->TopKSet.clear();
        unordered_set<long int>().swap(old->TopKSet);
        old->extremePoint.clear();
        vector<vector<float>>().swap(old->extremePoint);

        old->focaled.clear();
		unordered_set<long int>().swap(old->focaled);
		old->pairMap.clear();
		vector<vector<long int > >().swap(old->pairMap);
		delete old;
	}
}

void cellTree::lpModel(lprec* lp, int& dimen, unordered_map<long int, bool>& touchhs)
{
	double row[MAXDIMEN];
	if (lp == NULL)
	{
		fprintf(stderr, "Unable to create new LP model\n");
		exit(0);
	}
	float U_RANGE = 1000.0;
	// constraints in each dimension
	for (int di = 0; di < dimen; di++)
	{
		row[0] = 0;
		for (int dii = 1; dii < dimen + 1; dii++)
		{
			if (dii - 1 == di)
			{
				row[dii] = 1;
			}
			else
			{
				row[dii] = 0;
			}
		}

		add_constraint(lp, row, GE, 0);
		add_constraint(lp, row, LE, U_RANGE);
	}

	// in reduced space, sum_{q_i} should less than 1
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = 1;
	}
	add_constraint(lp, row, GE, 0);
	add_constraint(lp, row, LE, U_RANGE);

	// constraints in intersected hyperplanes
	for (unordered_map<long int, bool>::iterator iter = touchhs.begin(); iter != touchhs.end(); iter++)
	{
		addHP(lp, iter->first, iter->second);
	}

	// set scale 
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
}

void cellTree::lpModelwithSIDELEN(lprec* lp, int& dimen, unordered_map<long int, bool>& touchhs)
{
	double row[MAXDIMEN];
	if (lp == NULL)
	{
		fprintf(stderr, "Unable to create new LP model\n");
		exit(0);
	}

	// constraints in each dimension
	for (int di = 0; di < dimen; di++)
	{
		row[0] = 0;
		for (int dii = 1; dii < dimen + 1; dii++)
		{
			if (dii - 1 == di)
			{
				row[dii] = 1;
			}
			else
			{
				row[dii] = 0;
			}
		}
		add_constraint(lp, row, GE, SIDELEN*0.1);
		add_constraint(lp, row, LE, 1-SIDELEN*0.1);
	}

	// in reduced space, sum_{q_i} should less than 1
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = 1;
	}
	add_constraint(lp, row, GE, SIDELEN*0.1);
	add_constraint(lp, row, LE, 1-SIDELEN*0.1);

	// constraints in intersected hyperplanes
	for (unordered_map<long int, bool>::iterator iter = touchhs.begin(); iter != touchhs.end(); iter++)
	{
		addHPwithSIDELEN(lp, iter->first, iter->second);
	}

	// set scale 
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
}

void cellTree::addHP(lprec* model, long int hpid, bool sideindicator)
{
	int hsID = RecordIDtoHalfPlaneID[hpid] - 1;
	int dimen = HalfSpaces[0].size() - 2;
	double row[MAXDIMEN];

	row[0] = 0;
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = HalfSpaces[hsID][dii - 1];
	}
	if (sideindicator == false)
	{
		add_constraint(model, row, LE, HalfSpaces[hsID][dimen]);
	}
	else if (sideindicator == true)
	{
		add_constraint(model, row, GE, HalfSpaces[hsID][dimen]);
	}
	else
	{
		std::cout << "Unable to detect half plane direction!!!" << endl;
	}
}

void cellTree::addHPwithSIDELEN(lprec* model, long int hpid, bool sideindicator)
{
	int hsID = RecordIDtoHalfPlaneID[hpid] - 1;
	int dimen = HalfSpaces[0].size() - 2;


	double row[MAXDIMEN];

	row[0] = 0;
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = HalfSpaces[hsID][dii - 1];
	}
	if (sideindicator == false)
	{
		add_constraint(model, row, LE, HalfSpaces[hsID][dimen]-SIDELEN*0.1);
	}
	else if (sideindicator == true)
	{
		add_constraint(model, row, GE, HalfSpaces[hsID][dimen]+SIDELEN*0.1);
	}
	else
	{
		std::cout << "Unable to detect half plane direction!!!" << endl;
	}
}

bool cellTree::isFeasible(unordered_map<long int, bool>& touchhs, long int hpid, bool sideindicator)
{
	double row[MAXDIMEN];
	int dimen = HalfSpaces[0].size() - 2;

	lprec *lp = make_lp(0, dimen);
	
	lpModel(lp, dimen, touchhs);
	addHP(lp, hpid, sideindicator);

	row[0] = 0;
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = -1;
	}
	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);

	int ret = solve(lp);
	get_variables(lp, row);

	// for reduced space
	delete_lp(lp);
	if (ret == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}





// vector<vector<double> > cellTree::updateExtreme4d( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long int hpid, bool sideindicator, vector<vector<long int > >& pairMap)
void cellTree::updateExtreme4d( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long int hpid, bool sideindicator, vector<vector<long int > >& pairMap, vector<vector<float> > &candidate)
{


	int hsID = RecordIDtoHalfPlaneID[hpid] - 1;
	int dimen = HalfSpaces[0].size() - 2;

	// vector<vector<double> > candidate;
	candidate.clear();
	vector<vector<long int > > newMap;
	vector<long int> stillIndex;

	for (int i = 0; i < extremePoint.size(); ++i)
	{

		double dotResult = 0;
		for(int j = 0; j< dimen; j++)
			dotResult += extremePoint[i][j] * HalfSpaces[hsID][j];

		// vector<double> stillExtreme = {extremePoint[i][0], extremePoint[i][1], extremePoint[i][2]};
		vector<float> stillExtreme = {extremePoint[i][0], extremePoint[i][1], extremePoint[i][2]};

		if(sideindicator == false)
		{
			if(dotResult < HalfSpaces[hsID][dimen] )
			{
				candidate.push_back(stillExtreme);
				for(int j = 0 ; j < pairMap[i].size(); j++)
					stillIndex.push_back(pairMap[i][j]);
				newMap.push_back(pairMap[i]);
			}
		}
		else
		{
			if(dotResult > HalfSpaces[hsID][dimen]  )
			{
				candidate.push_back(stillExtreme);
				
				for(int j = 0 ; j < pairMap[i].size(); j++)
					stillIndex.push_back(pairMap[i][j]);
				newMap.push_back(pairMap[i]);
			}
		}	
				
		
	}
	
	int newPointCount = 0;

	vector<long int > IDset;
	for(auto iter = touchhs.begin(); iter != touchhs.end(); ++iter)
	{
		IDset.push_back(iter->first);
	}
	for (int i = 0; i < IDset.size()-1; ++i)
	{
		for(int j = i+1; j < IDset.size(); ++j)
		{
			long int ID = IDset[i];
			long int planeID = RecordIDtoHalfPlaneID[ID] - 1;

			long int ID2 = IDset[j];
			long int planeID2 = RecordIDtoHalfPlaneID[ID2] - 1;

			double a1 = HalfSpaces[planeID][0];
			double b1 = HalfSpaces[planeID][1];
			double c1 = HalfSpaces[planeID][2];
			double d1 = -HalfSpaces[planeID][3];

			double a2 = HalfSpaces[hsID][0];
			double b2 = HalfSpaces[hsID][1];
			double c2 = HalfSpaces[hsID][2];
			double d2 = -HalfSpaces[hsID][3];

			double a3 = HalfSpaces[planeID2][0];
			double b3 = HalfSpaces[planeID2][1];
			double c3 = HalfSpaces[planeID2][2];
			double d3 = -HalfSpaces[planeID2][3];
			double m1,m2,n1,n2,k1,k2,x,y,z;


			m1= -a1*b2 + a2*b1;
			m2= -a2*b3 + a3*b2;
			n1= -a1*c2  + a2*c1;
			n2=-a2*c3+a3*c2;
			k1= -a1*d2+a2*d1;
			k2= - a2*d3+a3*d2;
			
			z=(m2*k1-m1*k2)/(m1*n2-m2*n1 );
			y=(-k1 - n1*z)/(m1 );
			x=(-d1-c1*z-b1*y)/(a1 );

			//cout<<"before x, y: "<<x <<" "<<y<<endl;
			int flag = 1;
			//cout<<"before: "<<x<<" "<<y<<" "<<z<<endl;
			if (x> 0 && y>0  & z>0)
			{

				for(auto iiter = touchhs.begin(); iiter != touchhs.end(); ++iiter)
				{
					if(iiter->first == IDset[i]|| iiter->first == IDset[j] || iiter->first == hpid  )
						continue;
					double 	dotResult2 = 0;
					long int pendID = RecordIDtoHalfPlaneID[iiter->first] - 1;
					dotResult2 = x * HalfSpaces[pendID][0] + y * HalfSpaces[pendID][1] + z * HalfSpaces[pendID][2];

					if(dotResult2 == HalfSpaces[pendID][3])
						{
							flag = 1;
							break;
						}

					if(iiter->second == false)
					{
						
						if(dotResult2 < HalfSpaces[pendID][3] && flag == 1 )
							flag = 1;
						else
						{
							flag = 0;
							break;
						}
					}
					else
					{
						if(dotResult2 > HalfSpaces[pendID][3]  && flag == 1)
							flag = 1;
						else
						{
							flag = 0;
							break;
						}
					}

				}
				if(flag)
				{
				
					int flagSame = 1;
					for(int i = 0; i< candidate.size();i++)
					{
						if(abs(candidate[i][0]-x)<0.0001 && abs(candidate[i][1]-y)<0.0001 && abs(candidate[i][2]-z)<0.0001)
						{
							flagSame = 0;
							break;
						}
					}
					if(flagSame)
					{
						newPointCount++;
						// vector<double> newPoint={x,y, z};
						vector<float> newPoint={x,y, z};
						candidate.push_back(newPoint);

						stillIndex.push_back(hpid);
						stillIndex.push_back(IDset[i]);
						stillIndex.push_back(IDset[j]);
						vector<long int > newPair = {hpid, IDset[i], IDset[j]};
						newMap.push_back(newPair);
					

						// vector<double>().swap(newPoint);
						vector<float>().swap(newPoint);
					}


				}
			}
			else
			{
				continue;
			}
		}
	}
	pairMap.clear();
	for(int i = 0; i < newMap.size(); i++)
	{
		pairMap.push_back(newMap[i]);
	}

	unordered_map<long int, bool> temp;
	for(auto iter= touchhs.begin(); iter != touchhs.end(); iter++)
	{
		int flag = 0;
		for(int i = 0 ; i < stillIndex.size() ; i++)
		{
			if(iter->first == stillIndex[i])
			{
				flag = 1;
				break;
			}
		}
		if(flag )
		{
			temp[iter->first] = iter->second;
		}
	}
	touchhs.clear();
	for(auto iter = temp.begin(); iter != temp.end(); iter++)
	{
		touchhs[iter->first] = iter->second;
	}
	// return candidate;
} 


void cellTree::updateRank(cell* node, const int mink, long int hpid)
{
	long long int hpid_used = hpid%objCnt;
	if (hpid_used == 0) {
		hpid_used = objCnt;
	}
	if(node->TopKSet.find(hpid_used) == node->TopKSet.end()) 
	{
		node->TopKSet.insert(hpid_used);
		node->rank++;
		// // Hao added
		// // cout << "node->TopKSet.insert: " << hpid_used << ", rank from " << (node->rank - 1) << " to " << node->rank << endl;
		// for (auto & v: node->combAncestors[hpid_used]) {
		// 	if (node->TopKSet.find(v) == node->TopKSet.end()) {
		// 		node->TopKSet.insert(v);
		// 		node->rank++;
		// 		// // Hao added
		// 		// cout << "node->TopKSet.insert: " << v << ", rank from " << (node->rank - 1) << " to " << node->rank << endl;
		// 	}
		// }
	}

	if(node->BelowHP.find(hpid)== node->BelowHP.end())
		node->BelowHP.insert(hpid);	

	if (node->TopKSet.size() >= mink)
	{
		node->isPruned = true;
		
		if (node->left != NULL)
		{
			releaseCell(node->left);
			node->left = NULL;
		}
		if (node->right != NULL)
		{
			releaseCell(node->right);
			node->right = NULL;
		}
	}
	else
	{
		if (node->left != NULL)
			updateRank(node->left, mink,hpid);
		if (node->right != NULL)
			updateRank(node->right, mink,hpid);
	}
}

void cellTree::inserthp(long int & hpid, const int mink, cell* node, unordered_map<long int, bool>& touchhs, vector<cell*>& subregions)
{
	cout<<"not used in here"<<endl;
}

void cellTree::updateRelation(long long int & hpid, long int & focal, cell* node,  bool sideindicator)
{

	// unordered_set<long int> initSet;

	if(sideindicator == false)
	{
		long int anID = hpid%objCnt;
		if (anID == 0) {
			anID = objCnt;
		}

		// if (node->newAncestors.find(focal) == node->newAncestors.end())
		// {
		// 	node->newAncestors[focal] = initSet;
		// }

		// node->newAncestors[focal].insert(anID);
					
		// if(node->newdescendants.find(anID) == node->newdescendants.end())
		// 	node->newdescendants[anID] = initSet;
		// node->newdescendants[anID].insert(focal);

		// use new functions instead
		node->addNewAncestor(focal, anID);
		node->addNewDescendant(anID, focal);
		// cout << "addNewAncestor " << focal << ": " << anID << endl;
		// cout << "addNewDescendant " << anID << ": " << focal << endl;

		// // combine all the new relations into comb*
		// node->addCombAncestor(focal, anID);
		// // cout << "addCombAncestor " << focal << ": " << anID << endl;
		// for (auto & anc_of_anID: node->combAncestors[anID]) {
		// 	if (anc_of_anID != focal) {
		// 		node->addCombAncestor(focal, anc_of_anID);
		// 		// cout << "addCombAncestor " << focal << ": " << anc_of_anID << endl;
		// 	}
		// }
		// node->addCombDescendant(anID, focal);
		// // cout << "addCombDescendant " << anID << ": " << focal << endl;
		// for (auto & dec_of_focal: node->combDescendants[focal]) {
		// 	if (dec_of_focal != anID) {
		// 		node->addCombDescendant(anID, dec_of_focal);
		// 		// cout << "addCombDescendant " << anID << ": " << dec_of_focal << endl;
		// 	}
		// }
	}
	else
	{
		long int deID = hpid%objCnt;
		if (deID == 0) {
			deID = objCnt;
		}

		// if (node->newAncestors.find(deID) == node->newAncestors.end())
		// {
		// 	node->newAncestors[deID] = initSet;
		// }

		// node->newAncestors[deID].insert(focal);
		
		// if(node->newdescendants.find(focal) == node->newdescendants.end())
		// 	node->newdescendants[focal] = initSet;
		// node->newdescendants[focal].insert(deID);

		// use new functions instead
		node->addNewAncestor(deID, focal);
		node->addNewDescendant(focal, deID);
		// cout << "addNewAncestor " << deID << ": " << focal << endl;
		// cout << "addNewDescendant " << focal << ": " << deID << endl;

		// // combine all the new relations into comb*
		// node->addCombAncestor(deID, focal);
		// // cout << "addCombAncestor " << deID << ": " << focal << endl;
		// for (auto & anc_of_deID: node->combAncestors[deID]) {
		// 	if (anc_of_deID != focal) {
		// 		node->addCombAncestor(deID, anc_of_deID);
		// 		// cout << "addCombAncestor " << deID << ": " << anc_of_deID << endl;
		// 	}
		// }
		// node->addCombDescendant(focal, deID);
		// // cout << "addCombDescendant " << focal << ": " << deID << endl;
		// for (auto & dec_of_deID: node->combDescendants[deID]) {
		// 	if (dec_of_deID != focal) {
		// 		node->addCombDescendant(focal, dec_of_deID);
		// 		// cout << "addCombDescendant " << focal << ": " << dec_of_deID << endl;
		// 	}
		// }
	}
}


// vector<vector<double> > cellTree::updateExtreme( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long long int hpid, bool sideindicator, vector<vector<long int > >& pairMap)
void cellTree::updateExtreme( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long long int hpid, bool sideindicator, vector<vector<long int > >& pairMap, vector<vector<float> > &candidate)
{


    long long int hsID = RecordIDtoHalfPlaneID[hpid] - 1;

    int dimen = HalfSpaces[0].size() - 2;
    // vector<vector<double> > candidate;
    candidate.clear();
    vector<vector<long int > > newMap;
    vector<long int> stillIndex;


    for (int i = 0; i < extremePoint.size(); ++i)
    {

        double dotResult = 0;
        for(int j = 0; j< dimen; j++)
            dotResult += extremePoint[i][j] * HalfSpaces[hsID][j];

        // vector<double> stillExtreme = {extremePoint[i][0], extremePoint[i][1]};
        vector<float> stillExtreme = {extremePoint[i][0], extremePoint[i][1]};
        if(sideindicator == false)
        {
            if(dotResult < HalfSpaces[hsID][dimen]   )
            {
                candidate.push_back(stillExtreme);
                for(int j = 0 ; j < pairMap[i].size(); j++)
                    stillIndex.push_back(pairMap[i][j]);
                newMap.push_back(pairMap[i]);
            }
        }
        else
        {
            if(dotResult > HalfSpaces[hsID][dimen]  )
            {

                candidate.push_back(stillExtreme);

                for(int j = 0 ; j < pairMap[i].size(); j++)
                    stillIndex.push_back(pairMap[i][j]);
                newMap.push_back(pairMap[i]);
            }
        }
    }
    int index = candidate.size();


    //if(index == 0)
    //	return extremePoint;


    int newPointCount = 0;
    vector<long int > IDset;
    for(auto iter = touchhs.begin(); iter != touchhs.end(); ++iter)
    {
        IDset.push_back(iter->first);
    }
    for(  int i = 0; i < IDset.size(); i++)
    {
        long int ID = IDset[i];
        long int planeID = RecordIDtoHalfPlaneID[ID] - 1;
        //cout<<"planeID: "<<planeID<<endl;
        double a1 = HalfSpaces[planeID][0];
        double b1 = HalfSpaces[planeID][1];
        double c1 = -HalfSpaces[planeID][2];

        double a2 = HalfSpaces[hsID][0];
        double b2 = HalfSpaces[hsID][1];
        double c2 = -HalfSpaces[hsID][2];

        double y  = (c1 * a2 - c2 * a1) / (a1 * b2 - a2 * b1);
        double x  = (c2 * b1 - c1 * b2) / (a1 * b2 - a2 * b1);
        int flag = 1;

        if (x> 0 && y>0 )
        {

            for(auto iiter = touchhs.begin(); iiter != touchhs.end(); ++iiter)
            {
                if(iiter->first == IDset[i])
                    continue;
                double 	dotResult2 = 0;
                long int pendID = RecordIDtoHalfPlaneID[iiter->first] - 1;
                dotResult2 = x * HalfSpaces[pendID][0] + y * HalfSpaces[pendID][1];

                if(iiter->second == false)
                {
                    if( dotResult2 <= HalfSpaces[pendID][2]   )
                        flag = 1;
                    else
                    {
                        flag = 0;
                        break;
                    }
                }
                else
                {
                    if(dotResult2 >= HalfSpaces[pendID][2])
                        flag = 1;
                    else
                    {
                        flag = 0;
                        break;
                    }
                }

            }
            if(flag)
            {
                newPointCount++;
                // vector<double> newPoint={x,y};
                vector<float> newPoint={x,y};
                int flagSame = 0;
                for(int i=0; i < candidate.size(); i++)
                {

                    if( abs(candidate[i][0] - x) <= 0.00001 && abs(candidate[i][1] - y) <= 0.00001)
                    {
                        flagSame = 1;
                        break;
                    }
                }
                if(flagSame)
                {
                    newPointCount--;
                    continue;
                }

                candidate.push_back(newPoint);


                stillIndex.push_back(hpid);
                stillIndex.push_back(IDset[i]);
                vector<long int > newPair = {hpid, IDset[i] };
                newMap.push_back(newPair);

                if(newPointCount == dimen)
                    break;
            }
        }
        else
        {
            continue;
        }


    }

    pairMap.clear();
    for(int i = 0; i < newMap.size(); i++)
    {
        pairMap.push_back(newMap[i]);
    }

    unordered_map<long int, bool> temp;
    for(auto iter= touchhs.begin(); iter != touchhs.end(); iter++)
    {
        int flag = 0;
        for(int i = 0 ; i < stillIndex.size() ; i++)
        {
            if(iter->first == stillIndex[i])
            {
                flag = 1;
                break;
            }
        }
        if(flag )
        {
            temp[iter->first] = iter->second;
        }
    }
    touchhs.clear();
    for(auto iter = temp.begin(); iter != temp.end(); iter++)
    {
        touchhs[iter->first] = iter->second;
    }
    // return candidate;
}

// vector<vector<double> > cellTree::updateExtreme5d( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long int hpid, bool sideindicator, vector<vector<long int > >& pairMap)
void cellTree::updateExtreme5d( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long int hpid, bool sideindicator, vector<vector<long int > >& pairMap, vector<vector<float> > &candidate)
{


	int hsID = RecordIDtoHalfPlaneID[hpid] - 1;
	int dimen = HalfSpaces[0].size() - 2;
	// vector<vector<double> > candidate;
	candidate.clear();
	vector<vector<long int > > newMap;
	vector<long int> stillIndex;
	
	for (int i = 0; i < extremePoint.size(); ++i)
	{
	
		double dotResult = 0;
		for(int j = 0; j< dimen; j++)
			dotResult += extremePoint[i][j] * HalfSpaces[hsID][j];

		// vector<double> stillExtreme = {extremePoint[i][0], extremePoint[i][1], extremePoint[i][2], extremePoint[i][3]};
		vector<float> stillExtreme = {extremePoint[i][0], extremePoint[i][1], extremePoint[i][2], extremePoint[i][3]};

		if(sideindicator == false)
		{
			if(dotResult < HalfSpaces[hsID][dimen] )
			{
				candidate.push_back(stillExtreme);
				for(int j = 0 ; j < pairMap[i].size(); j++)
					stillIndex.push_back(pairMap[i][j]);
				newMap.push_back(pairMap[i]);
			}
		}
		else
		{
			if(dotResult > HalfSpaces[hsID][dimen]  )
			{
				candidate.push_back(stillExtreme);
				
				for(int j = 0 ; j < pairMap[i].size(); j++)
					stillIndex.push_back(pairMap[i][j]);
				newMap.push_back(pairMap[i]);
			}
		}	
				
		
	}
	
	int newPointCount = 0;

	vector<long int > IDset;
	for(auto iter = touchhs.begin(); iter != touchhs.end(); ++iter)
	{
		IDset.push_back(iter->first);
	}
	for (int i = 0; i < IDset.size()-2; ++i)
	{
		for(int j = i+1; j < IDset.size()-1; ++j)
		{
			for(int k = j+1; k < IDset.size(); ++k)
			{
				long int ID1 = IDset[i];
				long int planeID1 = RecordIDtoHalfPlaneID[ID1] - 1;

				long int ID2 = IDset[j];
				long int planeID2 = RecordIDtoHalfPlaneID[ID2] - 1;

				long int ID3 = IDset[k];
				long int planeID3 = RecordIDtoHalfPlaneID[ID3] - 1;
		
				double a10 = HalfSpaces[planeID1][0];
				double b10 = HalfSpaces[planeID1][1];
				double c10 = HalfSpaces[planeID1][2];
				double d10 = HalfSpaces[planeID1][3];
				double e10 = -HalfSpaces[planeID1][4];

				double a20 = HalfSpaces[hsID][0];
				double b20 = HalfSpaces[hsID][1];
				double c20 = HalfSpaces[hsID][2];
				double d20 = -HalfSpaces[hsID][3];
				double e20 = -HalfSpaces[hsID][4];

				double a30 = HalfSpaces[planeID2][0];
				double b30 = HalfSpaces[planeID2][1];
				double c30 = HalfSpaces[planeID2][2];
				double d30 = -HalfSpaces[planeID2][3];
				double e30 = -HalfSpaces[planeID2][4];

				double a40 = HalfSpaces[planeID3][0];
				double b40 = HalfSpaces[planeID3][1];
				double c40 = HalfSpaces[planeID3][2];
				double d40 = HalfSpaces[planeID3][3];
				double e40 = -HalfSpaces[planeID3][4];
				

				double a1, a2, a3, b1, b2, b3, c1,c2, c3, d1, d2,d3;


				a1 = a10*d20 - a20*d10;
				a2 = a10*d30 - a30*d10;
				a3 = a10*d40 - a40*d10;
				b1 = b10*d20 - b20*d10;
				b2 = b10*d30 - b30*d10;
				b3 = b10*d40 - b40*d10;
				c1 = c10*d20 - c20*d10;
				c2 = c10*d30 - c30*d10;
				c3 = c10*d40 - c40*d10;
				d1 = e10*d20 - e20*d10;
				d2 = e10*d30 - e30*d10;
				d3 = e10*d40 - e40*d10;
				double m1,m2,n1,n2,k1,k2,x,y,z, w;


				m1= -a1*b2 + a2*b1;
				m2= -a2*b3 + a3*b2;
				n1= -a1*c2  + a2*c1;
				n2=-a2*c3+a3*c2;
				k1= -a1*d2+a2*d1;
				k2= - a2*d3+a3*d2;
				
				z=(m2*k1-m1*k2)/(m1*n2-m2*n1);
				y=(-k1 - n1*z)/(m1);
				x=(-d1-c1*z-b1*y)/(a1);
				w = (-e10 -c10*z - b10*y - a10*x )/(d10);
				

				int flag = 1;
				if (x> 0 && y>0  & z>0 & w>0 )
				{

					for(auto iiter = touchhs.begin(); iiter != touchhs.end(); ++iiter)
					{
						if(iiter->first == IDset[i] || iiter->first == IDset[j] || iiter->first == IDset[k] || iiter->first == hpid)
							continue;
						double 	dotResult2 = 0;
						long int pendID = RecordIDtoHalfPlaneID[iiter->first] - 1;
						dotResult2 = x * HalfSpaces[pendID][0] + y * HalfSpaces[pendID][1] + z * HalfSpaces[pendID][2] + w * HalfSpaces[pendID][3];

						if(iiter->second == false)
						{
							if(dotResult2 == HalfSpaces[pendID][4])
							{
								flag = 1;
								break;
							}
							if(dotResult2 < HalfSpaces[pendID][4] && flag == 1 )
								flag = 1;
							else
							{
								flag = 0;
								break;
							}
						}
						else
						{
							if(dotResult2 == HalfSpaces[pendID][4])
							{
								flag = 1;
								break;
							}
							if(dotResult2 > HalfSpaces[pendID][4]  && flag == 1)
								flag = 1;
							else
							{
								flag = 0;
								break;
							}
						}

					}
					if(flag)
					{
					
						int flagSame = 1;
						for(int c = 0; c< candidate.size();c++)
						{
							if(abs(candidate[c][0]-x) < 0.001 && abs(candidate[c][1]- y)< 0.001 && abs(candidate[c][2]- z)< 0.001 && abs(candidate[c][3]- w )< 0.001)
							{
								flagSame = 0;
								break;
							}
						}
						if(flagSame)
						{
							newPointCount++;
							// vector<double> newPoint={x,y, z, w};
							vector<float> newPoint={x,y, z, w};
							candidate.push_back(newPoint);

							stillIndex.push_back(hpid);
							stillIndex.push_back(ID1);
							stillIndex.push_back(ID2);
							stillIndex.push_back(ID3);
							vector<long int > newPair = {hpid, ID1, ID2, ID3};
							newMap.push_back(newPair);
						

							// vector<double>().swap(newPoint);
							vector<float>().swap(newPoint);
						}


					}
				}
				else
				{
					continue;
				}
			}
		}
	}

	pairMap.clear();
	for(int i = 0; i < newMap.size(); i++)
	{
		pairMap.push_back(newMap[i]);
	}

	unordered_map<long int, bool> temp;
	for(auto iter= touchhs.begin(); iter != touchhs.end(); iter++)
	{
		int flag = 0;
		for(int i = 0 ; i < stillIndex.size() ; i++)
		{
			if(iter->first == stillIndex[i])
			{
				flag = 1;
				break;
			}
		}
		if(flag )
		{
			temp[iter->first] = iter->second;
		}
	}
	touchhs.clear();
	for(auto iter = temp.begin(); iter != temp.end(); iter++)
	{
		touchhs[iter->first] = iter->second;
	}
	// return candidate;
} 

int cellTree::pendingIntersection(vector<vector<float> >& extremePoint, long long int hpid)
{
	int flag = -1;
	long long int hsID = RecordIDtoHalfPlaneID[hpid] - 1;
	if(hsID<0) return -1;
	int dimen = HalfSpaces[0].size() - 2;

	for (int i = 0; i < extremePoint.size(); ++i)
	{

		double dotResult = 0;
		for(int j = 0; j< dimen; j++)
			dotResult += extremePoint[i][j] * HalfSpaces[hsID][j];
		if (abs(dotResult -HalfSpaces[hsID][dimen]) <0.01 )
			continue;



		if(dotResult < HalfSpaces[hsID][dimen] && flag != 0)
		{
			flag = 1;
		}
		
		else if(dotResult > HalfSpaces[hsID][dimen] && flag != 1)
		{
			flag = 0;
		}
		else
		{
			return 2;
		}
		
	}
	return flag;
}

// vector<vector<double> > cellTree::updateExtreme2d( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long int hpid, bool sideindicator, vector<vector<long int > >& pairMap)
void cellTree::updateExtreme2d( vector<vector<float> > &extremePoint, unordered_map<long int, bool>& touchhs,long int hpid, bool sideindicator, vector<vector<long int > >& pairMap, vector<vector<float> > &candidate)
{


	int hsID = RecordIDtoHalfPlaneID[hpid] - 1;
	int dimen = HalfSpaces[0].size() - 2;
	// vector<vector<double> > candidate;
	candidate.clear();
	vector<vector<long int > > newMap;
	vector<long int> stillIndex;

	for (int i = 0; i < extremePoint.size(); ++i)
	{
	
		double dotResult = 0;
		for(int j = 0; j< dimen; j++)
			dotResult += extremePoint[i][j] * HalfSpaces[hsID][j];

		// vector<double> stillExtreme = {extremePoint[i][0]};
		vector<float> stillExtreme = {extremePoint[i][0]};

		if(sideindicator == false)
		{
			if(dotResult < HalfSpaces[hsID][dimen] )
			{

				int flagSame = 1;
				for(int i=0; i < candidate.size(); i++)
				{
					if(abs(stillExtreme[0] - candidate[i][0]) < 0.01)
					{
						flagSame = 0;
						break;
					}
				}
				if(flagSame)
					candidate.push_back(stillExtreme);

			}
		}
		else
		{
			if(dotResult > HalfSpaces[hsID][dimen]  )
			{
				int flagSame = 1;
				for(int i=0; i < candidate.size(); i++)
				{
					if(abs(stillExtreme[0] - candidate[i][0]) < 0.01)
					{
						flagSame = 0;
						break;
					}
				}
				if(flagSame)
				candidate.push_back(stillExtreme);

			}
		}	
				
		
	}

	float a2 = HalfSpaces[hsID][0];
	float b2 = HalfSpaces[hsID][1];
	double x= b2/a2;
	// vector<double> newPoint={abs(x)};
	vector<float> newPoint={abs(x)};

	int flagSame = 1;
	for(int i=0; i < candidate.size(); i++)
	{
		if(abs(newPoint[0] - candidate[i][0]) < 0.01)
		{
			flagSame = 0;
			break;
		}
	}
	if(flagSame)
		candidate.push_back(newPoint);
	// return candidate;
} 


void cellTree::inserthp2(long int long & hpid, const int mink, cell* node, unordered_map<long int, bool>& touchhs, vector<cell*>& subregions, int k , float maxSP, const vector<int>& protectedAttribute, long int focal, int numOfKind)
{
	if (node->isPruned == false)
	{
		int flag = pendingIntersection(node->extremePoint, hpid);
		
		if (flag == 0)
		{
			updateRank(node, mink, hpid);
			updateRelation(hpid, focal, node, false);
		} 
		else if (flag == 1)
		{
			node->AboveHP.insert(hpid);
			updateRelation(hpid, focal, node, true);
		}
		else if (flag == 2)
		{
			if (node->left == NULL&&node->right == NULL)
			{
				if (node->TopKSet.size() >= mink)   
				{
					node->isPruned = true;
				}
				else
				{
					int dim = HalfSpaces[0].size() - 2+1;
				
					node->left = new cell(node);
					node->right = new cell(node);
					
					// // commented by Hao 20221217
					// vector<vector<double> > newExtremeLeft;
					// vector<vector<double> > newExtremeRight;

					
					// if(dim == 2)
					// 	newExtremeLeft= updateExtreme2d( node->extremePoint,  node->left->IntersectHP,hpid, false, node->left->pairMap);
					// else if(dim == 3)
					// 	newExtremeLeft= updateExtreme( node->extremePoint,  node->left->IntersectHP,hpid, false, node->left->pairMap);
					// else if(dim == 4)
					// 	newExtremeLeft= updateExtreme4d( node->extremePoint,  node->left->IntersectHP,hpid, false, node->left->pairMap);
					// else if(dim == 5)
					// 	newExtremeLeft= updateExtreme5d( node->extremePoint,  node->left->IntersectHP,hpid, false, node->left->pairMap);
					
					
					// if(dim == 2)
					// 	newExtremeRight=updateExtreme2d( node->extremePoint,  node->right->IntersectHP,hpid, true, node->right->pairMap);
					// else if(dim == 3)
					// 	newExtremeRight=updateExtreme( node->extremePoint,  node->right->IntersectHP,hpid, true, node->right->pairMap);
					// else if(dim == 4)
					// 	newExtremeRight=updateExtreme4d( node->extremePoint,  node->right->IntersectHP,hpid, true, node->right->pairMap);
					// else if(dim == 5)
					// 	newExtremeRight=updateExtreme5d( node->extremePoint,  node->right->IntersectHP,hpid, true, node->right->pairMap);
					

			
					// node->left->extremePoint = newExtremeLeft;
					// node->left->IntersectHP[hpid] = false;
					// vector<vector<double> >().swap(newExtremeLeft);

					// node->right->IntersectHP[hpid] = true;
					// node->right->extremePoint = newExtremeRight;
					// vector<vector<double> >().swap(newExtremeRight);
					// // commented by Hao 20221217

					// added by Hao 20221217
					if(dim == 2)
						updateExtreme2d( node->extremePoint,  node->left->IntersectHP,hpid, false, node->left->pairMap, node->left->extremePoint);
					else if(dim == 3)
						updateExtreme( node->extremePoint,  node->left->IntersectHP,hpid, false, node->left->pairMap, node->left->extremePoint);
					else if(dim == 4)
						updateExtreme4d( node->extremePoint,  node->left->IntersectHP,hpid, false, node->left->pairMap, node->left->extremePoint);
					else if(dim == 5)
						updateExtreme5d( node->extremePoint,  node->left->IntersectHP,hpid, false, node->left->pairMap, node->left->extremePoint);
					
					
					if(dim == 2)
						updateExtreme2d( node->extremePoint,  node->right->IntersectHP,hpid, true, node->right->pairMap, node->right->extremePoint);
					else if(dim == 3)
						updateExtreme( node->extremePoint,  node->right->IntersectHP,hpid, true, node->right->pairMap, node->right->extremePoint);
					else if(dim == 4)
						updateExtreme4d( node->extremePoint,  node->right->IntersectHP,hpid, true, node->right->pairMap, node->right->extremePoint);
					else if(dim == 5)
						updateExtreme5d( node->extremePoint,  node->right->IntersectHP,hpid, true, node->right->pairMap, node->right->extremePoint);
					
					node->left->IntersectHP[hpid] = false;

					node->right->IntersectHP[hpid] = true;
					// added by Hao 20221217
					
					updateRelation(hpid, focal, node->left, true);
					updateRelation(hpid, focal, node->right, false); // they are original under if (node->right->TopKSet.size() >= mink) {...}

					long long int hpid_used = hpid%objCnt;
					if (hpid_used == 0) {
						hpid_used = objCnt;
					}
					node->right->TopKSet.insert(hpid_used);
					// // Hao added
					// cout << "node->right->TopKSet.insert: " << hpid_used << endl;

					// for (auto & v: node->right->combAncestors[hpid_used]) {
					// 	if (node->right->TopKSet.find(v) == node->right->TopKSet.end()) {
					// 		node->right->TopKSet.insert(v);
					// 		// // Hao added
					// 		// cout << "node->right->TopKSet.insert: " << v << endl;
					// 	}
					// }

					if (node->right->TopKSet.size() >= mink)
					{
						node->right->isPruned = true;
					}
				}
			}
			else if (node->left != NULL&&node->right != NULL)
			{
				if (node->left->isPruned == false)
				{
					unordered_map<long int, bool> leftths = touchhs;
					inserthp2(hpid, mink, node->left, leftths, subregions, k, maxSP, protectedAttribute, focal, numOfKind);
				}
				if (node->right->isPruned == false)
				{
					unordered_map<long int, bool> rightths = touchhs;
					inserthp2(hpid, mink, node->right, rightths, subregions, k, maxSP, protectedAttribute, focal, numOfKind);
				}
				if (node->left->isPruned == true && node->right->isPruned == true)
				{
					node->isPruned = true;
				}
			}		
		}
	}
}

void print_tree_leaf(const cell* tree_leaf) {
	cout << "\t\tIntersectHP (" << tree_leaf->IntersectHP.size() << ")";
	for (auto & v: tree_leaf->IntersectHP) {
		cout << " (" << v.first << ", " << v.second << ")";
	}
	cout << endl;

	cout << "\t\tTopKSet (" << tree_leaf->TopKSet.size() << ")";
	for (auto & v: tree_leaf->TopKSet) {
		cout << " " << v;
	}
	cout << endl;

	cout << "\t\tnewAncestors:" << endl;
	for (auto & v: tree_leaf->newAncestors) {
		cout << "\t\t\t" << v.first << ":";
		for (auto & vv: v.second) {
			cout << " " << vv;
		}
		cout << endl;
	}

	// cout << "\t\tcombAncestors:" << endl;
	// for (auto & v: tree_leaf->combAncestors) {
	// 	cout << "\t\t\t" << v.first << ":";
	// 	for (auto & vv: v.second) {
	// 		cout << " " << vv;
	// 	}
	// 	cout << endl;
	// }

	cout << "\t\tnewdescendants:" << endl;
	for (auto & v: tree_leaf->newdescendants) {
		cout << "\t\t\t" << v.first << ":";
		for (auto & vv: v.second) {
			cout << " " << vv;
		}
		cout << endl;
	}

	// cout << "combDescendants:" << endl;
	// for (auto & v: tree_leaf->combDescendants) {
	// 	cout << "\t" << v.first << ":";
	// 	for (auto & vv: v.second) {
	// 		cout << " " << vv;
	// 	}
	// 	cout << endl;
	// }
}

// void cellTree::insert(vector<long long int> &hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult, unordered_set<long int> ancestors,vector<long int> descendants,vector<cell*>& subregions, long int focal, int k , float maxSP, vector<int>& protectedAttribute, int numOfKind)
void cellTree::insert(vector<long long int> &hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult, const unordered_set<long int>& ancestors, vector<cell*>& subregions, long int focal, int k , float maxSP, const vector<int>& protectedAttribute, int numOfKind)
{
	unordered_map<long int, bool> hs;
	for(auto iter = ancestors.begin(); iter != ancestors.end(); iter++ )
	{
		if(root->isPruned == false)
		{
			updateRank(root, mink, *iter);
		}
		else
			break;
	}
	for (int pos = 0; pos < hps.size(); pos++) 
	{

		if(root->isPruned == false)
		{
			if(ancestors.size() > 0)
			{
				long long int hps_pos_id = hps[pos]%objCnt;
				if (hps_pos_id == 0) {
					hps_pos_id = objCnt;
				}

				if(find(ancestors.begin(), ancestors.end(), hps_pos_id)  == ancestors.end())
				{
					inserthp2(hps[pos], mink, root, hs, subregions, k , maxSP, protectedAttribute, focal, numOfKind);

					// // Hao added
					// cout << "print tree leaves after inserthp2 " << hps[pos] << endl;
					// vector<cell*> leaves;
					// this->collectLeaf(leaves, INFINITY);
					// for (int li = 0; li < leaves.size(); li++) {
					// 	cout << "\tleaves[" << li << "]" << endl;
					// 	print_tree_leaf(leaves[li]);
					// }
				}
				else
				{
					continue;
				}
			}
			else
			{
				inserthp2(hps[pos], mink, root, hs, subregions, k , maxSP, protectedAttribute, focal, numOfKind);

				// // Hao added
				// cout << "print tree leaves after inserthp2 " << hps[pos] << endl;
				// vector<cell*> leaves;
				// this->collectLeaf(leaves, INFINITY);
				// for (int li = 0; li < leaves.size(); li++) {
				// 	cout << "\tleaves[" << li << "]" << endl;
				// 	print_tree_leaf(leaves[li]);
				// }
			}
		}
		else
		{
			updateCellTree(root);
			break;
		}
	}

	updateCellTree(root);
}

void cellTree::collectLeaf(vector<cell*>& leaves, const int& mink)
{
	leaves.clear();
	//if (root->isPruned == false)
	//{
		treeSize += root->cellSize();
		if (root->left != NULL)
		{
			cell* leaf = new cell();
			root->copyleaf(leaf);
			dfsTraversal(root->left, leaf, leaves);
		}
		if (root->right != NULL)
		{
			cell* leaf = new cell();
			root->copyleaf(leaf);
			dfsTraversal(root->right, leaf, leaves);
		}
		if (root->left == NULL&& root->right == NULL)
		{

			cell* leaf = new cell();
			root->copyleaf(leaf);
			leaves.push_back(leaf); 
		}
	//}
	sort(leaves.begin(), leaves.end(), cellCompare());
}

void cellTree::collectLeaf(vector<pair<cell*, long int>>& leaves, const int& mink)
{
	leaves.clear();
	if (root->isPruned == false)
	{
		if (root->left != NULL && root->left->isPruned == false)
		{
			cell* leaf = new cell();
			dfsTraversal(root->left, leaf, leaves);
		}
		if (root->right != NULL && root->right->isPruned == false)
		{
			cell* leaf = new cell();
			dfsTraversal(root->right, leaf, leaves);
		}
	}
	sort(leaves.begin(), leaves.end(), cellpairCompare());
}

void cellTree::dfsTraversal(cell* node, cell* leaf, vector<cell*>& leaves)
{
	//node->copyleaf(leaf);
	treeSize += node->cellSize();
	if (node->left == NULL && node->right == NULL)
	{
		node->copyleaf(leaf);
		leaves.push_back(leaf);
	}
	else
	{
		if (node->left != NULL)
		{
			cell* left = new cell();
			//leaf->copyleaf(left);
			dfsTraversal(node->left, left, leaves);
		}
		if (node->right != NULL)
		{
			cell* right = new cell();
			dfsTraversal(node->right, right, leaves);
		}
	}
}

void cellTree::dfsTraversal(cell* node, cell* leaf, vector<pair<cell*, long int>>& leaves)
{
	node->copyleaf(leaf);

	if (node->left == NULL && node->right == NULL)
	{
		leaves.push_back(make_pair(leaf,leaves.size()));
		cellID[leaves.size() - 1] = node;
	}
	else
	{
		if (node->left->isPruned == false)
		{
			cell* left = new cell();
			leaf->copyleaf(left);
			dfsTraversal(node->left, left, leaves);
		}
		if (node->right->isPruned == false)
		{
			dfsTraversal(node->right, leaf, leaves);
		}
		else
		{
			leaf->release();
			delete leaf;
		}
	}
}

void cellTree::opt_insert(vector<long int> &hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult)
{
	for (int pos = 0; pos < hps.size(); pos++)
	{
		if (root->isPruned == false)
		{
			unordered_map<long int, bool> leftths;
			//inserthp(hps[pos], mink, root->left, leftths); not used now --zz
			unordered_map<long int, bool> rightths;
			//inserthp(hps[pos], mink, root->right, rightths); not used now --zz
		}
	}
	updateCellTree(root);
}

void cellTree::opt_inserthp(long int & hpid, const int mink, cell* node, cell* all)
{
	if (node->isPruned == false)
	{
		node->appendto(all);
		int indicator = isDominatorInserted(dagNode[hpid]->rDominator, all);
		if (indicator == -1) // found negative halfspace from hpid's dominators
		{
			node->AboveHP.insert(hpid);
			totalSpaceCost += 4;
			releaseCell(all);
			delete all;
		}
		else if (indicator == 1) // found positive halfspace from hpid's dominators
		{
			//inserthp(hpid, mink, node, all->IntersectHP);
			cout<<"not used now --zz"<<endl;
		}
		else if (indicator == 0) // go to its children
		{
			if (node->left->isPruned == false)
			{
				cell* left = new cell();
				all->appendto(left);
				opt_inserthp(hpid, mink, node->left, left);
			}
			if (node->right->isPruned == false)
			{
				opt_inserthp(hpid, mink, node->left, all);
			}
			else 
			{
				node->isPruned = true;
				all->release();
				delete all;
			}
		}
		else
		{
			cout << "there is a bug in opt_inserthp" << endl;
		}
	}
}

void cellTree::maintainDAG(unordered_set<long int>& skylines, unordered_set<long int>& removeSL, float* PG[], const int dimen)
{
	if (removeSL.size() == 0)
	{
		for (auto iter = skylines.begin(); iter != skylines.end(); iter++)
		{
			gNode* tmp = new gNode((long int) (*iter));
			dagNode[*iter] = tmp;
		}
	}
	else
	{
		float oldR[MAXDIMEN], newR[MAXDIMEN];
		for (auto iter = skylines.begin(); iter != skylines.end(); iter++)
		{
			if (dagNode.find(*iter) == dagNode.end())
			{
				gNode* tmp = new gNode((long int)*iter);
				dagNode[*iter] = tmp;
				totalSpaceCost += 5;
				for (auto riter = removeSL.begin(); riter != removeSL.end(); riter++)
				{
					for (int i = 0; i < dimen; i++)
					{
						newR[i] = (PG[*iter][i] + PG[*iter][i + dimen]) / 2;
						oldR[i] = (PG[*riter][i] + PG[*riter][i + dimen]) / 2;
						if (dominateRecords(oldR, newR, dimen))
						{
							tmp->rDominator.insert(*riter);
						}
					}
				}
			}
			else
			{
				// this node has been inserted in dag
			}
		}
	}
}

int cellTree::isDominatorInserted(unordered_set<long int>& rdominators, cell* all)
{
	for (auto iter = rdominators.begin(); iter != rdominators.end(); iter++)
	{
		if (all->AboveHP.find(*iter) != all->AboveHP.end())
		{
			return -1;
		}
		if (all->BelowHP.find(*iter) != all->BelowHP.end())
		{
			return 1;
		}
		if (all->IntersectHP.find(*iter) != all->IntersectHP.end())
		{
			return all->IntersectHP[*iter] == true ? 1 : -1;
		}
	}
	return 0;
}

void cellTree::markSingular(vector<cell*>& leaves, unordered_set<long int>& a_skylines, const int& mink, vector<cell*>& finalResult, unordered_set<long int>& singular)
{
	long int recordID;
	a_skylines.clear();
	bool flag = false;

	for (int i = 0; i < leaves.size(); i++)
	{
		if (leaves[i]->isPruned == false)
		{
			if (leaves[i]->rank <= mink)
			{
				for (auto iter = leaves[i]->BelowHP.begin(); iter != leaves[i]->BelowHP.end(); iter++)
				{
					recordID = *iter;
					if (singular.find(recordID) == singular.end())
					{
						flag = true;
						if (a_skylines.find(recordID) == a_skylines.end())
						{
							a_skylines.insert(recordID);
						}
					}
				}

				for (auto iter = leaves[i]->IntersectHP.begin(); iter != leaves[i]->IntersectHP.end(); iter++)
				{
					if (iter->second == true)
					{
						recordID = iter->first;
						if (singular.find(recordID) == singular.end())
						{
							flag = true;
							if (a_skylines.find(recordID) == a_skylines.end())
							{
								a_skylines.insert(recordID);
							}
						}
					}
				}

				if (flag == false)
				{
					leaves[i]->isPruned = true;
					cell* tmp = new cell();
					leaves[i]->copyleaf(tmp);
					finalResult.push_back(tmp);
					return; // return as we are sure that record r will be utk;
				}
				leaves[i]->release();
				delete leaves[i];
			}
			else
			{
				leaves[i]->release();
				delete leaves[i];
			}
		}
	}
	for (auto iter = a_skylines.begin(); iter != a_skylines.end(); iter++)
	{
		singular.insert(*iter);
	}
}

void cellTree::markSingular(vector<pair<cell*, long int>>& leaves, unordered_set<long int>& a_skylines, const int& mink, vector<cell*>& finalResult, unordered_set<long int>& singular)
{
	long int recordID;
	bool flag = false;

	for (int i = 0; i < leaves.size(); i++)
	{
		if (leaves[i].first->isPruned == false)
		{
			if (leaves[i].first->rank < mink)
			{
				for (auto iter = leaves[i].first->BelowHP.begin(); iter != leaves[i].first->BelowHP.end(); iter++)
				{
					recordID = *iter;
					if (singular.find(recordID) == singular.end())
					{
						flag = true;
						if (a_skylines.find(recordID) == a_skylines.end())
						{
							a_skylines.insert(recordID);
						}
						break;
					}
				}

				for (auto iter = leaves[i].first->IntersectHP.begin(); iter != leaves[i].first->IntersectHP.end(); iter++)
				{
					if (iter->second == true)
					{
						recordID = iter->first;
						if (singular.find(recordID) == singular.end())
						{
							flag = true;
							if (a_skylines.find(recordID) == a_skylines.end())
							{
								a_skylines.insert(recordID);
							}
							break;
						}
					}
				}

				if (flag == false)
				{
					leaves[i].first->isPruned = true;
					cellID[leaves[i].second]->isPruned = true;
					cell* tmp = new cell();
					leaves[i].first->copyleaf(tmp);
					finalResult.push_back(tmp);
				}
			}
			else
			{
				leaves[i].first->release();
				
			}
		}
	}

}

void cellTree::updateSkyline(unordered_set<long int>& skylines, vector<pair<cell*, long int>>& leaves, unordered_set<long int>& singular)
{
	for (int i = 0; i < leaves.size(); i++)
	{
		if (leaves[i].first->isPruned == false)
		{
			for (auto iter = leaves[i].first->BelowHP.begin(); iter != leaves[i].first->BelowHP.end(); iter++)
			{
				if (singular.find(*iter) == singular.end())
				{
					singular.insert(*iter);
					if (skylines.find(*iter) != skylines.end())
					{
						skylines.erase(*iter);
					}
					else
					{
						cout << "please check what happened!" << endl;
					}
				}
			}

			for (auto iter = leaves[i].first->IntersectHP.begin(); iter != leaves[i].first->IntersectHP.end(); iter++)
			{
				if (iter->second == true && singular.find(iter->first) == singular.end())
				{
					singular.insert(iter->first);
					if (skylines.find(iter->first) != skylines.end())
					{
						skylines.erase(iter->first);
					}
					else
					{
						cout << "please check what happened!" << endl;
					}
				}
			}
		}
	}
}

void cellTree::updateCellTree(cell* node)
{
	if (node->isPruned == true || node == NULL)
		return ;
	else if (node->left == NULL&&node->right == NULL)
		return;
	else if (node->left->isPruned == true && node->right->isPruned == true)
	{
		node->isPruned = true;
		releaseCell(node->left);
		releaseCell(node->right);
		node->left = NULL;
		node->right = NULL;
		return;
	}
	else
	{
		updateCellTree(node->left);
		updateCellTree(node->right);
		if (node->left->isPruned == true && node->right->isPruned == true)
		{
			node->isPruned = true;
			releaseCell(node->left);
			releaseCell(node->right);
			node->left = NULL;
			node->right = NULL;
			return;
		}
	}
}

int cellTree::Lemma2(vector<cell*>& cells)
{
	int touching = 0;
	int count = 0;
	for (int i = 0; i < cells.size() && count<=100; i++)
	{
		double row[MAXDIMEN];
		int dimen = HalfSpaces[0].size() - 2;
		lprec *lp = make_lp(0, dimen);
		lpModel(lp, dimen, cells[i]->IntersectHP);
		touching += cells[i]->IntersectHP.size();

		row[0] = 0;
		for (int dii = 1; dii < dimen + 1; dii++)
		{
			row[dii] = -1;
		}
		set_obj_fn(lp, row);
		set_verbose(lp, IMPORTANT);
		set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);

		int ret = solve(lp);
		get_variables(lp, row);
		delete_lp(lp);
		count++;
	}
	return touching;
}

void cellTree::halfspace2polytope(cell* ret, string outfile, string hsfile)
{
	double row[MAXDIMEN];

	int dimen = HalfSpaces[0].size() - 2;
	lprec *lp = make_lp(0, dimen);

	lpModelwithSIDELEN(lp, dimen, ret->IntersectHP);

	row[0] = 0;
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = 1;
	}
	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);

	int tmp = solve(lp);
	get_variables(lp, row);

	if (tmp == 0)
	{
		// qhull
		FILE *fout1 = fopen(hsfile.c_str(), "w+");
		if (fout1 == NULL)
		{
			cout << "Error in opening file innerpoint.txt" << endl;
			getchar();
			exit(0);
		}

		fprintf(fout1, "%d 1\n", dimen);   //the dimensionality
		for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", row[i]); //the feasible point found by lp above   *****************
		fprintf(fout1, "\n");
		fprintf(fout1, "%d\n", dimen + 1);   //dimensionality + 1
		fprintf(fout1, "%d\n", ret->IntersectHP.size() + 2 * (dimen + 1));  //the total number of hyperplanes for intersection

		// constraints in each dimension
		for (int di = 0; di < dimen; di++)
		{
			row[0] = 0;
			for (int dii = 1; dii < dimen + 1; dii++)
			{
				if (dii - 1 == di)
				{
					row[dii] = 1;
				}
				else
				{
					row[dii] = 0;
				}
			}

			//add_constraint(lp, row, GE, 0);
			for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", -row[i + 1]); // dimension constraint  *****************
			fprintf(fout1, "0 \n ");
			for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", row[i + 1]); // dimension constraint *****************
			fprintf(fout1, "-1 \n ");
		}

		// in reduced space, sum_{q_i} should less than 1
		for (int dii = 1; dii < dimen + 1; dii++)
		{
			row[dii] = 1;
		}
		for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", -1 * row[i + 1]); // dimension constraint  *****************
		fprintf(fout1, "0 \n ");
		for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", row[i + 1]); // dimension constraint *****************
		fprintf(fout1, "-1 \n ");

		// constraints in intersected hyperplanes
		int hpid, hsID;
		bool sideindicator;
		for (unordered_map<long int, bool>::iterator iter = ret->IntersectHP.begin(); iter != ret->IntersectHP.end(); iter++)
		{
			hpid = iter->first;
			sideindicator = iter->second;
			hsID = RecordIDtoHalfPlaneID[hpid] - 1;

			row[0] = 0;
			for (int dii = 1; dii < dimen + 1; dii++)
			{
				row[dii] = HalfSpaces[hsID][dii - 1];
			}
			if (sideindicator == false)
			{
				for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", row[i + 1]); // dimension constraint *****************
				fprintf(fout1, "%f \n ", -1 * HalfSpaces[hsID][dimen]);
			}
			else if (sideindicator == true)
			{
				for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", -row[i + 1]); // dimension constraint  *****************
				fprintf(fout1, "%f \n ", HalfSpaces[hsID][dimen]);
			}
			else
			{
				std::cout << "Unable to detect half plane direction!!!" << endl;
			}
		}
		fclose(fout1);

		#ifdef MAC:
		char sys_string[4096] = "qhalf Fp < ";
		strcat(sys_string, hsfile.c_str());
		strcat(sys_string, " | qconvex FS > ");
		strcat(sys_string, outfile.c_str());
		#else
		char sys_string[4096] = "D:\\MonoTopK\\src\\qhull\\qhalf.exe Fp < ";
		strcat(sys_string, hsfile.c_str());
		strcat(sys_string, " | D:\\MonoTopK\\src\\qhull\\qconvex.exe FS > ");
		strcat(sys_string, outfile.c_str());
		#endif
		system(sys_string);
	}
	return;
}

void cellTree::rankBound(vector<pair<cell*, long int>>& leaves, Rtree& rt, Point& a_pt, const int& mink, vector<cell*>& finalResult)
{
	rBound itemR;

	for (int i = 0; i < leaves.size(); i++)
	{
		if (leaves[i].first->isPruned == false)
		{
			if (i % 100 == 0)
				cout << i << endl;
			cellBound(leaves[i].first, rt, a_pt, itemR, mink);
			if (itemR.minimum > mink)
			{
				leaves[i].first->isPruned = true;
			}
			else if (itemR.maximum <= mink)
			{
				leaves[i].first->isPruned = true;
				cellID[leaves[i].second]->isPruned = true;

				cell* tmp = new cell();
				leaves[i].first->copyleaf(tmp);
				tmp->rank = itemR.maximum;
				finalResult.push_back(tmp);
				
			}
		}
	}
}

void cellTree::rtreeRAM(Rtree& rt, unordered_map<long int, RtreeNode*>& ramTree)
{
	ramTree.clear();
	queue<long int> H;
	RtreeNode* node;
	H.push(rt.m_memory.m_rootPageID);
	long int pageID;
	while (!H.empty())
	{
		pageID = H.front();
		H.pop();
		node = rt.m_memory.loadPage(pageID);
		ramTree[pageID] = node;
		totalSpaceCost += 4104;
		if (node->isLeaf() == false)
		{
			for (int i = 0; i < node->m_usedspace; i++)
			{
				H.push(node->m_entry[i]->m_id);
			}
		}
	}
}

void cellTree::cellBound(cell* item, Rtree& rt, Point& a_pt, rBound& itemR, const int& mink)
{
	Score pScore, eScore;
	itemR.minimum = 0;
	itemR.maximum = objCnt - 1;
	int dimen = HalfSpaces[0].size() - 2;
	lprec* lp = make_lp(0, dimen);
	int isOptimal;
	lpModel(lp, dimen, item->IntersectHP);
	focalScore(lp, a_pt, pScore);

	vector<float> cl, cu;
	float clsum = 0, cusum = 0;
	findCellMBR(lp, cl, cu);
	for (int di = 0; di < dimen; di++)
	{
		clsum += cl[di];
		cusum += cu[di];
	}

    cl.push_back(max(1 - cusum, 0));
    cu.push_back(min(1 - clsum, 1));

	queue<long int> H;
	RtreeNode* node;
	H.push(rt.m_memory.m_rootPageID);
	long int pageID;
	bool isLeafNode;
	while (!H.empty())
	{
		pageID = H.front();
		H.pop();
		node = ramTree[pageID];
		if (node->isLeaf() == false)
		{
			isLeafNode = false;
			for (int i = 0; i < node->m_usedspace; i++)
			{
				approxiRecordScore(node->m_entry[i], eScore, cl, cu, isLeafNode);
				if (eScore.min >= pScore.max)
				{
					itemR.minimum += node->m_entry[i]->num_records;
					if (itemR.minimum > mink)
						return;
				}
				else if (eScore.max < pScore.min)
				{
					// the records in this rtree node will affect the rank of given cell "item"
					itemR.maximum -= node->m_entry[i]->num_records;
					if (itemR.maximum < mink)
						return;
				}
				else
				{
					H.push(node->m_entry[i]->m_id);
				}
			}
		}
		else
		{
			isLeafNode = true;
			for (int i = 0; i < node->m_usedspace; i++)
			{
				approxiRecordScore(node->m_entry[i], eScore, cl, cu, isLeafNode);
				if (eScore.min >= pScore.max)
				{
					itemR.minimum += 1;
					if (itemR.minimum > mink)
						return;
				}
				else if (eScore.max < pScore.min)
				{
					itemR.maximum -= 1;
					if (itemR.maximum < mink)
						return;
				}

			}
		}
	}
	delete_lp(lp);
}

void cellTree::focalScore(lprec* lp, Point& a_pt, Score& pScore)
{
	double row[MAXDIMEN];
	int dimen = HalfSpaces[0].size() - 2;
	int isOptimal;

	// obtain min score for focal record
	row[0] = 0;
	for (int dii = 0; dii < dimen; dii++)
	{
		row[dii + 1] = a_pt[dii] - a_pt[dimen];
	}
	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	isOptimal = solve(lp);
	get_variables(lp, row);

	if (isOptimal == 2)
	{
		cout << "focalScore LP_solver: please check what happened!!!" << endl;
	}

	pScore.min = a_pt[dimen];
	for (int dii = 0; dii < dimen; dii++)
	{
		pScore.min += (a_pt[dii] - a_pt[dimen]) * row[dii];
	}


	row[0] = 0;
	for (int dii = 0; dii < dimen; dii++)
	{
		row[dii + 1] = -(a_pt[dii] - a_pt[dimen]);
	}
	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	isOptimal = solve(lp);
	get_variables(lp, row);

	if (isOptimal == 2)
	{
		cout << "focalScore LP_solver: please check what happened!!!" << endl;
	}

	pScore.max = a_pt[dimen];
	for (int dii = 0; dii < dimen; dii++)
	{
		pScore.max += (a_pt[dii] - a_pt[dimen]) * row[dii];
	}
}

void cellTree::recordScore(lprec* lp, RtreeNodeEntry* e, Score& entryR, bool& isLeafNode)
{
	int dimen = HalfSpaces[0].size() - 2;
	double r_lower[MAXDIMEN];
	double r_upper[MAXDIMEN];
	int isOptimal;

	r_lower[0] = 0.0;
	r_upper[0] = 0.0;

	if (isLeafNode == false)
	{
		for (int di = 0; di < dimen; di++)
		{
			r_upper[di + 1] = -(e->m_hc.getUpper()[di] - e->m_hc.getUpper()[dimen]);
			r_lower[di + 1] = e->m_hc.getLower()[di] - e->m_hc.getLower()[dimen];
		}
	}
	else
	{
		for (int di = 0; di < dimen; di++)
		{
			r_upper[di + 1] = -((e->m_hc.getLower()[di] + e->m_hc.getUpper()[di]) / 2 - (e->m_hc.getUpper()[dimen] + e->m_hc.getLower()[dimen]) / 2);
			r_lower[di + 1] = ((e->m_hc.getLower()[di] + e->m_hc.getUpper()[di]) / 2 - (e->m_hc.getUpper()[dimen] + e->m_hc.getLower()[dimen]) / 2);
		}
	}

	set_obj_fn(lp, r_upper);
	set_verbose(lp, IMPORTANT);
	isOptimal = solve(lp);
	get_variables(lp, r_upper);

	if (isOptimal == 2)
	{
		cout << "recordScore LP_solver: please check what happened!!!" << endl;
	}

	if (isLeafNode == false)
	{
		entryR.max = e->m_hc.getUpper()[dimen];
		for (int dii = 0; dii < dimen; dii++)
		{
			entryR.max += (e->m_hc.getUpper()[dii] - e->m_hc.getUpper()[dimen]) * r_upper[dii];
		}
	}
	else
	{
		entryR.max = (e->m_hc.getUpper()[dimen] + e->m_hc.getLower()[dimen]) / 2;
		for (int dii = 0; dii < dimen; dii++)
		{
			entryR.max += ((e->m_hc.getUpper()[dii] + e->m_hc.getLower()[dii]) / 2 - (e->m_hc.getUpper()[dimen] + e->m_hc.getLower()[dimen]) / 2) * r_upper[dii];
		}
	}


	set_obj_fn(lp, r_lower);
	set_verbose(lp, IMPORTANT);
	isOptimal = solve(lp);
	get_variables(lp, r_lower);

	if (isOptimal == 2)
	{
		cout << "recordScore LP_solver: please check what happened!!!" << endl;
	}

	if (isLeafNode == false)
	{
		entryR.min = e->m_hc.getLower()[dimen];
		for (int dii = 0; dii < dimen; dii++)
		{
			entryR.min += (e->m_hc.getLower()[dii] - e->m_hc.getLower()[dimen]) * r_lower[dii];
		}
	}
	else
	{
		entryR.min = (e->m_hc.getUpper()[dimen] + e->m_hc.getLower()[dimen]) / 2;
		for (int dii = 0; dii < dimen; dii++)
		{
			entryR.min += ((e->m_hc.getUpper()[dii] + e->m_hc.getLower()[dii]) / 2 - (e->m_hc.getUpper()[dimen] + e->m_hc.getLower()[dimen]) / 2) * r_lower[dii];
		}
	}
}

void cellTree::findCellMBR(lprec* lp, vector<float>& cl, vector<float>& cu)
{
	double row[MAXDIMEN];
	int dimen = HalfSpaces[0].size() - 2;
	int isOptimal;
	cl.clear();
	cu.clear();

	// obtain min score for focal record
	for (int i = 0; i < dimen; i++)
	{
		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (dii == i)
				row[dii + 1] = 1;
			else
				row[dii + 1] = 0;
		}
		set_obj_fn(lp, row);
		set_verbose(lp, IMPORTANT);
		isOptimal = solve(lp);
		get_variables(lp, row);
		if (isOptimal == 2)
		{
			cout << "focalScore LP_solver: please check what happened!!!" << endl;
		}
		else
		{
			cl.push_back(row[i]);
		}


		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (dii == i)
				row[dii + 1] = -1;
			else
				row[dii + 1] = 0;
		}
		set_obj_fn(lp, row);
		set_verbose(lp, IMPORTANT);
		isOptimal = solve(lp);
		get_variables(lp, row);
		if (isOptimal == 2)
		{
			cout << "focalScore LP_solver: please check what happened!!!" << endl;
		}
		else
		{
			cu.push_back(row[i]);
		}
	}
}

void cellTree::approxiRecordScore(RtreeNodeEntry* e, Score& entryR, vector<float>& cl, vector<float>& cu, bool& isLeafNode)
{
	int dimen = HalfSpaces[0].size() - 2;
	entryR.max = 0;
	entryR.min = 0;

	if (isLeafNode == false)
	{
		for (int dii = 0; dii < dimen + 1; dii++)
		{
			entryR.max += e->m_hc.getUpper()[dii] * cu[dii];
			entryR.min += e->m_hc.getLower()[dii] * cl[dii];
		}
	}
	else
	{
		for (int dii = 0; dii < dimen + 1; dii++)
		{
			entryR.max += (e->m_hc.getUpper()[dii] + e->m_hc.getLower()[dii]) / 2 * cu[dii];
			entryR.min += (e->m_hc.getUpper()[dii] + e->m_hc.getLower()[dii]) / 2 * cl[dii];
		}
	}
}

void cellTree::rankBoundOpt(vector<pair<cell*, long int>>& leaves, Rtree& rt, Point& a_pt, const int& mink, vector<cell*>& finalResult)
{
	rBound itemR;

	for (int i = 0; i < leaves.size(); i++)
	{
		if (leaves[i].first->isPruned == false)
		{
			if (i % 100 == 0)
				cout << i << endl;
			cellBoundOpt(leaves[i].first, rt, a_pt, itemR, mink);
			if (itemR.minimum > mink)
			{
				leaves[i].first->isPruned = true;
			}
			else if (itemR.maximum <= mink)
			{
				leaves[i].first->isPruned = true;
				cellID[leaves[i].second]->isPruned = true;

				cell* tmp = new cell();
				leaves[i].first->copyleaf(tmp);
				tmp->rank = itemR.maximum;
				finalResult.push_back(tmp);
			}
		}
	}
}

void cellTree::cellBoundOpt(cell* item, Rtree& rt, Point& a_pt, rBound& itemR, const int& mink)
{
	itemR.minimum = 0;
	itemR.maximum = 0;
	int dimen = HalfSpaces[0].size() - 2;
	lprec* lp = make_lp(0, dimen);
	int isOptimal;
	lpModel(lp, dimen, item->IntersectHP);
	double obj = 0;

	vector<float> cl, cu;
	float clsum = 0, cusum = 0;
	findCellMBR(lp, cl, cu);
	
	for (int di = 0; di < dimen; di++)
	{
		clsum += cl[di];
		cusum += cu[di];
	}

	cl.push_back(max(1 - cusum, 0));
	cu.push_back(min(1 - clsum, 1));
	Score objScore;

	queue<long int> H;
	RtreeNode* node;
	H.push(rt.m_memory.m_rootPageID);
	long int pageID;
	bool isLeafNode;
	while (!H.empty())
	{
		pageID = H.front();
		H.pop();
		node = ramTree[pageID];
		if (node->isLeaf() == false)
		{
			isLeafNode = false;
			for (int i = 0; i < node->m_usedspace; i++)
			{
				fastOptObj(cu, cl, node->m_entry[i], a_pt, isLeafNode, objScore);
				if (objScore.min > 0)
				{
					itemR.minimum += node->m_entry[i]->num_records;
					itemR.maximum += node->m_entry[i]->num_records;
				}
				else if (objScore.max > 0)
				{
					itemR.maximum += node->m_entry[i]->num_records;
				}
				else
				{
					if ((obj = lpOptObj(lp, node->m_entry[i], a_pt, isLeafNode, false)) >= 0)
					{
						itemR.minimum += node->m_entry[i]->num_records;
						itemR.maximum += node->m_entry[i]->num_records;
					}
					else if ((obj = lpOptObj(lp, node->m_entry[i], a_pt, isLeafNode, true)) > 0)
					{
						itemR.maximum += node->m_entry[i]->num_records;
					}
					else
					{
						H.push(node->m_entry[i]->m_id);
					}
				}
			}
		}
		else
		{
			isLeafNode = true;
			for (int i = 0; i < node->m_usedspace; i++)
			{
				fastOptObj(cu, cl, node->m_entry[i], a_pt, isLeafNode, objScore);
				if (objScore.min >= 0)
				{
					itemR.minimum += 1;
					itemR.maximum += 1;
				}
				else if (objScore.max > 0)
				{
					itemR.maximum += 1;
				}
				else
				{
					if ((obj = lpOptObj(lp, node->m_entry[i], a_pt, isLeafNode, false)) >= 0)
					{
						itemR.minimum += 1;
						itemR.maximum += 1;
					}
					else if ((obj = lpOptObj(lp, node->m_entry[i], a_pt, isLeafNode, true)) > 0)
					{
						itemR.maximum += 1;
					}
				}
			}
		}
	}
	delete_lp(lp);
}

double cellTree::lpOptObj(lprec* lp, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, const bool& isMax)
{
	double ret;
	double row[MAXDIMEN];
	double optobj[MAXDIMEN];
	int dimen = HalfSpaces[0].size() - 2;

	if (isMax == false)
	{
		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (isLeafNode == false)
				row[dii + 1] = (e->m_hc.getLower()[dii] + -e->m_hc.getLower()[dimen]) - (a_pt[dii] - a_pt[dimen]);
			else
				row[dii + 1] = ((e->m_hc.getLower()[dii] + e->m_hc.getUpper()[dii]) / 2 - (e->m_hc.getLower()[dimen] + e->m_hc.getUpper()[dimen]) / 2) - (a_pt[dii] - a_pt[dimen]);
		}
		if (isLeafNode == false)
		{
			ret = e->m_hc.getLower()[dimen] - a_pt[dimen];
		}
		else
		{
			ret = (e->m_hc.getLower()[dimen] + e->m_hc.getUpper()[dimen]) / 2 - a_pt[dimen];
		}
		
	}
	else
	{
		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (isLeafNode == false)
				row[dii + 1] = -((e->m_hc.getUpper()[dii] - e->m_hc.getUpper()[dimen]) - (a_pt[dii] - a_pt[dimen]));
			else
				row[dii + 1] = -(((e->m_hc.getLower()[dii] + e->m_hc.getUpper()[dii]) / 2 - (e->m_hc.getLower()[dimen] + e->m_hc.getUpper()[dimen]) / 2) - (a_pt[dii] - a_pt[dimen]));
		}
		
		if (isLeafNode == false)
		{
			ret = (e->m_hc.getUpper()[dimen] - a_pt[dimen]);
		}
		else
		{
			ret = (e->m_hc.getLower()[dimen] + e->m_hc.getUpper()[dimen]) / 2 - a_pt[dimen];
		}
	}

	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);

	int tmp = solve(lp);
	get_variables(lp, optobj);
	if (tmp != 0)
	{
		cout << "optimized objective has errors, please check!" << endl;
	}
	else
	{
		double tmp = ret + get_objective(lp);
		if (isMax)
		{
			for (int dii = 0; dii < dimen; dii++)
			{
				ret += -1 * row[dii + 1] * optobj[dii];
			}
		}
		else
		{
			for (int dii = 0; dii < dimen; dii++)
			{
				ret += row[dii + 1] * optobj[dii];
			}
		}
	}
	return ret;
}

void cellTree::fastOptObj(vector<float>& cu, vector<float>& cl, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, Score& objScore)
{
	double weight = 0.0;
	double maxw = 0.0;
	double minw = 0.0;
	int dimen = HalfSpaces[0].size() - 1;
	
	for (int dii = 0; dii < dimen; dii++)
	{
		if (isLeafNode == false)
		{
			minw = (e->m_hc.getLower()[dii] - a_pt[dii]);
			maxw = (e->m_hc.getUpper()[dii] - a_pt[dii]);

			if (maxw < 0)
			{
				objScore.max += minw*cl[dii];
				objScore.min += maxw*cu[dii];
			}
			else if (minw <= 0)
			{
				objScore.max += maxw*cu[dii];
				objScore.min += minw*cu[dii];
			}
			else
			{
				objScore.max += maxw*cu[dii];
				objScore.min += minw*cl[dii];
			}
		}
		else
		{
			weight = ((e->m_hc.getLower()[dii] + e->m_hc.getUpper()[dii]) / 2 - a_pt[dii]);
			if (weight >= 0)
			{
				objScore.max += weight*cu[dii];
				objScore.min += weight*cl[dii];
			}
			else
			{
				objScore.max += weight*cl[dii];
				objScore.min += weight*cu[dii];
			}
		}
	}
	
}

void cellTree::cellBoundOptFULL(cell* item, Rtree& rt, Point& a_pt, rBound& itemR, const int& mink)
{
	itemR.minimum = 0;
	itemR.maximum = 0;
	int dimen = HalfSpaces[0].size() - 2;
	lprec* lp = make_lp(0, dimen);
	int isOptimal;
	lpModel(lp, dimen, item->IntersectHP);
	double obj = 0;

	vector<float> cl, cu;
	float clsum = 0, cusum = 0;
	findCellMBR(lp, cl, cu);
	Score objScore;

	queue<long int> H;
	RtreeNode* node;
	H.push(rt.m_memory.m_rootPageID);
	long int pageID;
	bool isLeafNode;
	while (!H.empty())
	{
		pageID = H.front();
		H.pop();
		node = ramTree[pageID];
		if (node->isLeaf() == false)
		{
			isLeafNode = false;
			for (int i = 0; i < node->m_usedspace; i++)
			{
				fastOptObjFULL(cu, cl, node->m_entry[i], a_pt, isLeafNode, objScore);
				if (objScore.min > 0)
				{
					itemR.minimum += node->m_entry[i]->num_records;
					itemR.maximum += node->m_entry[i]->num_records;
				}
				else if (objScore.max > 0)
				{
					itemR.maximum += node->m_entry[i]->num_records;
				}
				else
				{
					if ((obj = lpOptObjFULL(lp, node->m_entry[i], a_pt, isLeafNode, false)) >= 0)
					{
						itemR.minimum += node->m_entry[i]->num_records;
						itemR.maximum += node->m_entry[i]->num_records;
					}
					else if ((obj = lpOptObjFULL(lp, node->m_entry[i], a_pt, isLeafNode, true)) > 0)
					{
						itemR.maximum += node->m_entry[i]->num_records;
					}
					else
					{
						H.push(node->m_entry[i]->m_id);
					}
				}
			}
		}
		else
		{
			isLeafNode = true;
			for (int i = 0; i < node->m_usedspace; i++)
			{
				fastOptObjFULL(cu, cl, node->m_entry[i], a_pt, isLeafNode, objScore);
				if (objScore.min >= 0)
				{
					itemR.minimum += 1;
					itemR.maximum += 1;
				}
				else if (objScore.max > 0)
				{
					itemR.maximum += 1;
				}
				else
				{
					if ((obj = lpOptObjFULL(lp, node->m_entry[i], a_pt, isLeafNode, false)) >= 0)
					{
						itemR.minimum += 1;
						itemR.maximum += 1;
					}
					else if ((obj = lpOptObjFULL(lp, node->m_entry[i], a_pt, isLeafNode, true)) > 0)
					{
						itemR.maximum += 1;
					}
				}
			}
		}
	}
	delete_lp(lp);
}

double cellTree::lpOptObjFULL(lprec* lp, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, const bool& isMax)
{
	double ret = 0.0;
	double row[MAXDIMEN];
	double optobj[MAXDIMEN];
	int dimen = HalfSpaces[0].size() - 2;

	if (isMax == false)
	{
		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (isLeafNode == false)
				row[dii + 1] = (e->m_hc.getLower()[dii] - a_pt[dii]);
			else
				row[dii + 1] = ((e->m_hc.getLower()[dii] + e->m_hc.getUpper()[dii]) / 2 - a_pt[dii]);
		}
	}
	else
	{
		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (isLeafNode == false)
				row[dii + 1] = -(e->m_hc.getUpper()[dii] - a_pt[dii]);
			else
				row[dii + 1] = -((e->m_hc.getLower()[dii] + e->m_hc.getUpper()[dii]) / 2 - a_pt[dii]);
		}

	}

	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);

	int tmp = solve(lp);
	get_variables(lp, optobj);
	if (tmp != 0)
	{
		cout << "optimized objective has errors, please check!" << endl;
	}
	else
	{
		double tmp = ret + get_objective(lp);
		if (isMax)
		{
			for (int dii = 0; dii < dimen; dii++)
			{
				ret += -1 * row[dii + 1] * optobj[dii];
			}
		}
		else
		{
			for (int dii = 0; dii < dimen; dii++)
			{
				ret += row[dii + 1] * optobj[dii];
			}
		}
	}
	return ret;
}

void cellTree::rankBoundOptFULL(vector<pair<cell*, long int>>& leaves, Rtree& rt, Point& a_pt, const int& mink, vector<cell*>& finalResult)
{
	rBound itemR;

	for (int i = 0; i < leaves.size(); i++)
	{
		if (leaves[i].first->isPruned == false)
		{
			if (i % 100 == 0)
				cout << i << endl;
			cellBoundOptFULL(leaves[i].first, rt, a_pt, itemR, mink);
			if (itemR.minimum > mink)
			{
				leaves[i].first->isPruned = true;
			}
			else if (itemR.maximum <= mink)
			{
				leaves[i].first->isPruned = true;
				cellID[leaves[i].second]->isPruned = true;

				cell* tmp = new cell();
				leaves[i].first->copyleaf(tmp);
				tmp->rank = itemR.maximum;
				finalResult.push_back(tmp);
			}
		}
	}
}

void cellTree::fastOptObjFULL(vector<float>& cu, vector<float>& cl, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, Score& objScore)
{
	double weight = 0.0;
	double maxw = 0.0;
	double minw = 0.0;
	int dimen = HalfSpaces[0].size() - 2;

	for (int dii = 0; dii < dimen; dii++)
	{
		if (isLeafNode == false)
		{
			minw = (e->m_hc.getLower()[dii] - a_pt[dii]);
			maxw = (e->m_hc.getUpper()[dii] - a_pt[dii]);

			if (maxw < 0)
			{
				objScore.max += minw*cl[dii];
				objScore.min += maxw*cu[dii];
			}
			else if (minw <= 0)
			{
				objScore.max += maxw*cu[dii];
				objScore.min += minw*cu[dii];
			}
			else
			{
				objScore.max += maxw*cu[dii];
				objScore.min += minw*cl[dii];
			}
		}
		else
		{
			weight = ((e->m_hc.getLower()[dii] + e->m_hc.getUpper()[dii]) / 2 - a_pt[dii]);
			if (weight >= 0)
			{
				objScore.max += weight*cu[dii];
				objScore.min += weight*cl[dii];
			}
			else
			{
				objScore.max += weight*cl[dii];
				objScore.min += weight*cu[dii];
			}
		}
	}

}
