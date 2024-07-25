#include "header.h"
#include "hypercube.h"
#include "rentry.h"
#include "rnode.h"
#include "rtree.h"
#include "filemem.h"
#include "tgs.h"
#include "global.h"
#include "skyline.h"
#include "cellTree.h"
#include "param.h"
#include "utk.h"
#include "queryGen.h"
#include "arraY.h"
#include "QuadProg++.h"

#include <sstream>
using namespace std;
/////////////////////////////////

void initContainer(int& insertCount, vector<cell*>& finalResult, vector<long>& newAddSL, multimap<float, int>& heap, 
	unordered_set<long int>& skylines, unordered_set<long int>& removeSL, unordered_set<long int>& singular,
	unordered_map<long int, long int>& RecordIDtoHalfPlaneID, vector<vector<float>>& HalfSpaces, unordered_set<long int>& ignorerecords);

// gobal variables
vector<vector<float>> HalfSpaces; // halfspace 
unordered_map<long long int, long long  int> RecordIDtoHalfPlaneID;  //  record ID to halfplane ID
unordered_map<long int, RtreeNode*> ramTree; // load Rtree to main-memory
unordered_map<long int, cell*> cellID; // cell ID

int objCnt = 0; // # of data objects
double totalIO = 0; // # of IO access
double totalSpaceCost = 0.0; // space cost (MB)
double treeSpacecost = 0.0;
double Spacecost = 0.0;

float U_RANGE = 1000.0;


bool next_arr(vector<int>& curr_arr, int max_steps)
{
	// cout << "curr_arr:";
	// for (auto &v: curr_arr) {
	// 	cout << " " <<  v;
	// }
	// cout << endl;

	int pos_full = -1;
	int sum = 0;
	do {
		pos_full++;
		sum += curr_arr[pos_full];
	} while (sum < max_steps && pos_full < curr_arr.size() - 1);
	// cout << "sum: " << sum << endl;
	// cout << "pos_full: " << pos_full << endl;

	if (pos_full == 0 && sum == max_steps) {
		return false; // which means curr_arr[0] == max_steps, the curr_arr is already the ending
	} else {
		if (sum < max_steps) {
			curr_arr[pos_full]++;
		} else {
			curr_arr[pos_full - 1]++;
			if (pos_full < curr_arr.size()) {
				curr_arr[pos_full] = 0;
			}
		}
		return true;
	}
}

float u_dist(vector<float>& F1, vector<float>& F2)
{
	float d2 = 0.0;
	for (int i = 0; i < F1.size(); i++) {
		d2 += (F1[i] - F2[i]) * (F1[i] - F2[i]);
	}
	return sqrt(d2);
}

void load_domains(int dim, int num_levels, const string& domainfile, vector<vector<int>>& generalization, vector<vector<vector<float>>>& partitionPoint)
{
    // int maxj = 11;
    int maxj = num_levels + 1;
    int maxk = 51;

    // int arr[]={11,11,11,11,11};
    int arr[]={maxj,maxj,maxj,maxj,maxj};
    int b[5] = {-1, -1, -1, -1, -1 };
    c(arr,b,dim,0 ,generalization);
    sort(generalization.begin(), generalization.end(), cmp2);
    partitionPoint.resize(dim);
    
    for (int n = 0; n < dim; n++) {
        partitionPoint[n].resize(maxj);
    }

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < maxj; j++) {
            partitionPoint[i][j].resize(maxk);
        }
    }
    
    ifstream infile(domainfile);
    string s;
    if(!infile) {
        cout<<domainfile<<" not opened"<<endl;
    }
    int i = 0;
    int j = 0;
    float num;
    while(getline(infile, s)) {
        istringstream iss(s);
        int ijk = 0;
        while(iss >> num) {
            partitionPoint[i][j][ijk] = num;
            ijk++;
        }
        j++;
        if(j == arr[i]) {
            j = 0;
            i++;
        }
    }
}

void get_top_k_in_main(vector<float>& utility_func, const vector<long int>& subset, float** PointSet, int dim, int k, vector<long int>& returned_top_k, vector<float>* utilities = nullptr)
{
	vector<pair<float, long int>> tuple_pairs;
	for (const auto & ti: subset) {
	// for (long int ti = 1; ti <= objCnt; ti++) {
		float utility = 0.0;
		for (int di = 0; di < dim; di++) {
			utility += utility_func[di] * (PointSet[ti][di] + PointSet[ti][di + dim]) / 2.0;
		}
		tuple_pairs.push_back(make_pair(utility, ti));
	}
	sort(tuple_pairs.begin(), tuple_pairs.end(), [](pair<float, long int> const& a, pair<float, long int> const& b) -> bool {
		return a.first > b.first;
	});
	for(int ki = 0; ki < k; ki++) {
		returned_top_k.push_back(tuple_pairs[ki].second);

		if (utilities) {
			utilities->push_back(tuple_pairs[ki].first);
		}
	}
}

float get_sp_theory(int k, int numOfKind, vector<int>& kindCount)
{
	int sum = 0;
	for (int i = 0; i < numOfKind; i++) {
		sum += kindCount[i];
	}

	if (numOfKind > 2) {
		return float(k) / float(sum);
	} else {
		float max_sp = -1.0;
		for (int x = 1; x < k; x++) {
			float sp = min(float(x) / kindCount[0], float(k - x) / kindCount[1]);
			if (max_sp < sp) {
				max_sp = sp;
			}
		}
		return max_sp;
	}
}


// void run_brute_force()
// {

// }


int main(const int argc, const char** argv)
{
	// //cout.precision(4);

	int algorithm = 0; // 0 -> exact, 1 -> approx1, 2 -> approx2, 3 -> brute-force, 4 -> measure-fa*ir
	int ablation = 0;
	if (argc > 1) {
		if (string(argv[1]) == "appr1") {
			algorithm = 1;
		} else if (string(argv[1]) == "appr2") {
			algorithm = 2;
		} else if (string(argv[1]) == "bf") {
			algorithm = 3;
		} else if (string(argv[1]) == "measure-fa_ir") {
			algorithm = 4;
		} else if (string(argv[1]) == "nospeedup1") {
			ablation = 1;
		} else if (string(argv[1]) == "nospeedup2") {
			ablation = 2;
		} else if (string(argv[1]) == "nospeedup") {
			ablation = 3;
		}
	}

	if (algorithm == 4 && argc < 4) {
		cerr << "Measure-fa*ir: please supply the two files for fa*ir results" << endl;
		exit(1);
	}

	clock_t at, ad;

	/*
	// parameter parser
	// cout << "Parse Parameters" << endl;
	if (argc == 1)
	{
		helpmsg(argv[0]);
		return -1;
	}
	
	int dim = atoi(Param::read(argc, argv, "-d", ""));
	int k = atoi(Param::read(argc, argv, "-k", ""));
	const char* datafile = Param::read(argc, argv, "-f", "");
	const char* indexfile = Param::read(argc, argv, "-i", "");
	const char* domainfile = Param::read(argc, argv, "-p", "");
	const float sigma = atof(Param::read(argc, argv, "-R", ""));
	const float alpha = atof(Param::read(argc, argv, "-a", ""));
	int numOfKind = atoi(Param::read(argc, argv, "-m", ""));;
	*/

	int dim, k, numOfKind, TotalSize;
	char datafile[255];
	char domainfile[255];
	float alpha, beta;
	const char* indexfile = "index";
	const float sigma = 0.999;
	int theta;
	int num_levels;
	fstream fpConfig;
	fpConfig.open("config.txt", ios::in);
	// fpConfig.open("config_simple.txt", ios::in);
	if (!fpConfig) {
		cout << "config file not found" << endl;
		exit(0);
	}

	fpConfig >> datafile;
	cout << "Dataset file: " << datafile << endl;
	fpConfig >> TotalSize;
	cout << "Size of the dataset (n): " << TotalSize << endl;
	fpConfig >> k;
	cout << "k: " << k << endl;
	fpConfig >> dim;
	cout << "Number of dimensions (d): " << dim << endl;
	fpConfig >> numOfKind;
	cout << "Number of groups (|P|): " << numOfKind << endl;
	fpConfig >> alpha;
	cout << "alpha: " << alpha << endl;
	fpConfig >> beta;
	cout << "beta: " << beta << endl;
	// fpConfig >> theta;
	// cout << "theta: " << theta << endl;
	fpConfig >> num_levels;
	cout << "Number of generalized domains in each dimension (K): " << num_levels << endl;
	fpConfig >> domainfile;
	cout << "Domain file: " << domainfile << endl;

	int resultSize = 0;

	vector<vector<float>>  Queryregions;
	QueryGenerator(Queryregions, dim - 1, sigma);

	// for (auto & qr: Queryregions) {
	// 	for (auto & q: qr) {
	// 		cout << q << " ";
	// 	}
	// 	cout << endl;
	// }
	// exit(0);
	
	// data loading
	cout << "Load data points from file" << endl;
	float** PointSet = new float*[MAXPTS + 1];
	RtreeNodeEntry** rtreeNodeEntry = new RtreeNodeEntry*[MAXPTS];
	fstream fpdata;
	fpdata.open(datafile, ios::in);
	if (!fpdata) {
		cout << "data file not found" << endl;
		exit(0);
	}

	int id = 1;
	vector<int> protectedAttribute;
	srand((unsigned)time(NULL));
	int kind;

	vector<int> kindCount = {0,0,0,0,0,0,0,0,0,0,0,0};

	string fpdata_line;
	while (getline(fpdata, fpdata_line)) {
		stringstream ss(fpdata_line);
		float* cl = new float[dim];
		float* cu = new float[dim];

		PointSet[objCnt + 1] = new float[2 * dim];

		for (int d = 0; d < dim; d++) {
			ss >> cu[d];
			PointSet[objCnt + 1][d+dim] = cu[d];
		}
		for (int d = 0; d < dim; d++) {
			ss >> cl[d];
			PointSet[objCnt + 1][d] = cl[d];
		}
		ss >> kind;
		Hypercube hc(dim, cl, cu);
		rtreeNodeEntry[objCnt] = new RtreeNodeEntry(id, hc, 1);
		
		objCnt++;
		protectedAttribute.push_back(kind);
		kindCount[kind]++;
		id++;

		// log information
		// if (objCnt % 1000 == 0)
		// 	cout << ".";
		if (objCnt % 10000 == 0)
			cout << objCnt << " objects loaded" << endl;
	}

	// cout << "Print data points:" << endl;
	// for (int pi = 1; pi <= objCnt; pi++) {
	// 	for (int di = 0; di < dim; di++) {
	// 		cout << PointSet[pi][di] << " ";
	// 	}
	// 	for (int di = 0; di < dim; di++) {
	// 		cout << PointSet[pi][di + dim] << " ";
	// 	}
	// 	cout << endl;
	// }
	
	// cout << "Num in each protected group";
	// for(int c = 0; c < numOfKind; c++) {
	// 	cout << " " << kindCount[c];
	// }
	// cout << endl;

	// cout << "Total number of objects: " << objCnt << endl;

	double rawSize = dataSize(objCnt, dim);
	totalSpaceCost += rawSize;

	// build rtree
	// cout << "Bulkloading R-tree..." << endl;
	const int maxChild = (PAGESIZE - RtreeNode::size()) / RtreeNodeEntry::size(dim);
	FileMemory mem(PAGESIZE, indexfile, RtreeNodeEntry::fromMem, true);
	Rtree* rtree = TGS::bulkload(mem, dim, maxChild, maxChild, (int) maxChild * 0.3, (int) maxChild * 0.3, rtreeNodeEntry, objCnt, false);
	// cout << "[Rtree build done]" << endl;

	rtreeRAM(*rtree, ramTree);
	totalSpaceCost += ramTree.size() * 4096.00 / MB;

	treeSpacecost = 0.0;

	clock_t alls, alle;
	alls = clock();

	// size = (K+1)^d, which is the total number of global domains, while each member has size = d
	// For example, if K = 2, d = 2, then this list is (0 0) (1 0) (2 0) (0 1) ... (2 2)
	// For another example, if K = 1, d = 3, then this list is (0 0 0) (1 0 0) ... (1 1 1)
    vector<vector<int>> generalization;
    // size = d, each partitionPoint[di] size = (K + 1)
    // For example, if K = 1, d = 3, then each partitionPoint[di] has (K + 1) = 2 lists: <0.5, 1> and <1>, representing [0, 0.5)&[0.5, 1] and [0, 1]
    // Note that the number 0 in generalization[0][di] representing the first list (i.e., [0, 0.5)&[0.5, 1])
    vector<vector<vector<float>>> partitionPoint;
    num_levels -= 1; // to fix a issue that we actually have (K + 1) generalized domains in each dimension
	load_domains(dim, num_levels, domainfile, generalization, partitionPoint);

    // cout << "domains: " << generalization.size() << endl;

    // for (auto & g: generalization) {
    // 	for (auto & gg: g) {
    // 		cout << gg << " ";
    // 	}
    // 	cout << endl;
    // }

    // for (auto & one_g: partitionPoint) {
	//     for (auto & p: one_g) {
	//     	for (auto & pp: p) {
	//     		cout << pp << " ";
	//     	}
	//     	cout << endl;
	//     }
    // }

	/* seems not used
		// unordered_map<int, cell*> utkRet;
		// Point pt;
		// int Noofdominator;
		// int updatek;
		// unordered_set<long int> skylines;
		// unordered_set<long int> ignorerecords;
		// unordered_set<long int> removeSL;
		// unordered_set<long int> singular;
		// multimap<float, int> heap;
		// vector<long> newAddSL;
		// int insertCount = 0;
		// vector<cell*> leaves;
		// vector<cell*> finalResult;
		// utkRet.clear();
	*/

    // Queryregions.size() == 1 must be true
    auto & regions = Queryregions[0];

	UTK *sol = new UTK();
	vector<cell*> exactutk;
	vector<unordered_map<long int, bool>> intersection;

	vector<long int> rskyband;
	vector<long int> skyband;
	vector<long int> klayers;

	kskyband(dim, *rtree, skyband, PointSet, k);
	onionlayer(skyband, PointSet, k, klayers, dim);
	sol->onionToRskyband(regions, dim, *rtree, rskyband, PointSet, k, klayers);

	// run below to replace the above only for small testing data
	// for (long int obj_i = 1; obj_i <= objCnt; obj_i++) {
	// 	rskyband.push_back(obj_i);
	// 	skyband.push_back(obj_i);
	// 	klayers.push_back(obj_i);
	// }

	// cout << "skyband(" << skyband.size() << ") for k = " << k;
	// // for (auto & v: skyband) {
	// // 	cout << " " << v;
	// // }
	// cout << endl;

	// cout << "klayers(" << klayers.size() << ")";
	// // for (auto & v: klayers) {
	// // 	cout << " " << v;
	// // }
	// cout << endl;

	// cout << "# k-onion layer records / # of k-skyband: " << klayers.size() << "/" << skyband.size() << endl;

	// cout << "rskyband(" << rskyband.size() << ")";
	// // for (auto & v: rskyband) {
	// // 	cout << " " << v << " (" << protectedAttribute[v - 1] << ")";
	// // }
	// cout << endl;

	fstream fData;
	fData.open("originalF.txt", ios::in);
	if (!fData) {
		cout << "Original Utility Function not found" << endl;
		exit(0);
	}
	vector<float> originalF;  
	float originalF_numb;
	float originalF_sum = 0;
	for(int di = 0; di < dim-1; di++) {
		fData >> originalF_numb;
		originalF.push_back(originalF_numb * U_RANGE);
		originalF_sum += originalF_numb * U_RANGE;
	}
	originalF.push_back(U_RANGE - originalF_sum);
	fData.close();

	cout << "Original Utility Function:";
	for (auto & v: originalF) {
		cout << " " << (v / U_RANGE);
	}
	cout << endl;

	vector<float> originalF_norm;
	for(int di = 0; di < dim; di++) {
		originalF_norm.push_back(originalF[di] / U_RANGE);
	}

	vector<long int> original_TopKSet;
	vector<float> original_TopKUtilities;
	get_top_k_in_main(originalF_norm, rskyband, PointSet, dim, k, original_TopKSet, &original_TopKUtilities);
	// cout << "Original top-k:";
	// for (int ki = 0; ki < k; ki++) {
	// 	cout << " " << original_TopKSet[ki] << " (" << original_TopKUtilities[ki] << ")";
	// }
	// cout << endl;

	float original_fn = getAlphaFairness(original_TopKSet, PointSet, protectedAttribute, dim, generalization, partitionPoint, alpha, numOfKind, kindCount, num_levels);
	// cout << "Original Fairness: " << original_fn << endl;

	float SPinTheory = get_sp_theory(k, numOfKind, kindCount);
	// cout << "SPinTheory: " << SPinTheory << endl;

	cout << endl;

	if (algorithm == 4) {
		cout << "============================================" << endl;
		cout << "Collecting FA*IR" << endl;

		const string fa_ir_stat = argv[2];
		const string fa_ir_time = argv[3];

		ifstream ifs_fa_ir_stat(fa_ir_stat);
		if (ifs_fa_ir_stat.is_open()) {
			vector<long int> fa_ir_top_k_list;
			vector<int> fa_ir_top_k_grps;
			int vio_count = 0;
			for (int ki = 0; ki < k; ki++) {
				int fa_ir_top_k, fa_ir_top_k_grp;
				ifs_fa_ir_stat >> fa_ir_top_k >> fa_ir_top_k_grp;
				fa_ir_top_k_list.push_back(fa_ir_top_k);
				fa_ir_top_k_grps.push_back(fa_ir_top_k_grp);
			}
			ifs_fa_ir_stat >> vio_count;

			float fa_ir_fn = getAlphaFairness(fa_ir_top_k_list, PointSet, protectedAttribute, dim, generalization, partitionPoint, alpha, numOfKind, kindCount, num_levels);

			ofstream ofs_fa_ir_stat("stat_fa_ir.txt");
			ofs_fa_ir_stat << -100001 << endl;
			ofs_fa_ir_stat << fa_ir_fn << endl;
			ofs_fa_ir_stat << -100001 << endl;
			for (int di = 0; di < dim; di++) {
				ofs_fa_ir_stat << -1 << " ";
			}
			ofs_fa_ir_stat << endl;
			for (int ki = 0; ki < k; ki++) {
				ofs_fa_ir_stat << fa_ir_top_k_list[ki] << " ";
			}
			ofs_fa_ir_stat << endl;
			ofs_fa_ir_stat << vio_count << endl;

			ofs_fa_ir_stat.close();
			ifs_fa_ir_stat.close();
		} else {
			cout << "FA*IR result file " << fa_ir_stat << " not found" << endl;
		}

		ifstream ifs_fa_ir_time(fa_ir_time);
		if (ifs_fa_ir_time.is_open()) {
			float fa_ir_time;
			ifs_fa_ir_time >> fa_ir_time;

			ofstream ofs_fa_ir_time("time_fa_ir.txt");
			ofs_fa_ir_time << fa_ir_time;

			ofs_fa_ir_time.close();
			ifs_fa_ir_time.close();
		} else {
			cout << "FA*IR result file " << fa_ir_time << " not found" << endl;
		}

	}
	else if (algorithm == 3) {
		cout << "============================================" << endl;
		cout << "Brute-force starts" << endl;

		const int num_steps = 1000;
		// const int num_steps = 10;

		float max_SP_bf = -1.0;
		float min_dist_bf = 9999999.99;
		vector<float> best_f_bf(dim, 0.0);

		int test_dim = dim;
		vector<int> arr_bf(test_dim - 1, 0);
		while (next_arr(arr_bf, num_steps)) {
			// cout << "next_arr:";
			// for (auto &v: arr_bf) {
			// 	cout << " " << v;
			// }
			// cout << endl;

			vector<float> testF(test_dim, 0.0);
			int arr_sum = 0;
			for (int di = 0; di < test_dim - 1; di++) {
				testF[di] = float(arr_bf[di]) / float(num_steps);
				arr_sum += arr_bf[di];
			}
			testF[test_dim - 1] = float(num_steps - arr_sum) / float(num_steps);

			// cout << "testF:";
			// for (auto &v: testF) {
			// 	cout << " " << v;
			// }
			// cout << endl;

			// try to restore the top-k set under testF
			// however, I found that the top-k set come from rskyband, which is the intersection of skyband and klayers (in onionlayer)
			// but actually, it *should* be the union of skyband and klayers?
			vector<long int> test_TopKSet;
			get_top_k_in_main(testF, rskyband, PointSet, test_dim, k, test_TopKSet);
			// cout << "test top-k:";
			// for (int ki = 0; ki < k; ki++) {
			// 	cout << " " << test_TopKSet[ki]; // << " (" << tuple_pairs[ki].first << ")";
			// }
			// cout << endl;

			unordered_set<long int> test_TopKSet_set(test_TopKSet.begin(), test_TopKSet.end());

			float SP = getAlphaFairness(test_TopKSet_set, PointSet, protectedAttribute, test_dim , generalization, partitionPoint, alpha, numOfKind, kindCount, num_levels);
			// cout << "SP: " << SP << endl;
			float new_dist = u_dist(testF, originalF_norm);

			if (max_SP_bf < SP) {
				max_SP_bf = SP;
				min_dist_bf = new_dist;
				for (int di = 0; di < dim; di++) {
					best_f_bf[di] = testF[di];
				}
			} else if (max_SP_bf == SP) {
				if (min_dist_bf > new_dist) {
					min_dist_bf = new_dist;
					for (int di = 0; di < dim; di++) {
						best_f_bf[di] = testF[di];
					}
				}
			}
		}

		cout << "max_SP_bf: " << max_SP_bf << endl;
		cout << "min_dist_bf: " << min_dist_bf << endl;
		cout << "best F:";
		for (int di = 0; di < dim; di++) {
			cout << " " << best_f_bf[di];
		}
		cout << endl;

		vector<long int> best_top_k_bf;
		vector<float> best_TopKUtilities;
		get_top_k_in_main(best_f_bf, rskyband, PointSet, dim, k, best_top_k_bf, &best_TopKUtilities);
		cout << "best top-k:";
		for (int ki = 0; ki < k; ki++) {
			cout << " " << best_top_k_bf[ki] << " (" << best_TopKUtilities[ki] << ")";
		}
		cout << endl;
		cout << "============================================" << endl;

		string fn_stat = "stat_ftq_bf_exact.txt";
		// if (algorithm == 1) {
		// 	fn_stat = "stat_ftq_appr1.txt";
		// } else if (algorithm == 2) {
		// 	fn_stat = "stat_ftq_appr2.txt";
		// } else {
		// 	fn_stat = "stat_ftq_bf_exact.txt";
		// }
			
		FILE* stream_stat = fopen(fn_stat.c_str(), "w");

		fprintf(stream_stat, "%f\n", min_dist_bf);
		fprintf(stream_stat, "%f\n", max_SP_bf);
		fprintf(stream_stat, "%f\n", 0.0); // fake loss
		for (int di = 0; di < dim; di++) {
			fprintf(stream_stat, "%f ", best_f_bf[di]);
		}
		fprintf(stream_stat, "\n");
		for (int ki = 0; ki < k; ki++) {
			fprintf(stream_stat, "%d ", best_top_k_bf[ki]);
		}
		fprintf(stream_stat, "\n");

	}
	else {

		if (algorithm == 0) {
			cout << "Run FairTQ-Exact" << endl;
		} else if (algorithm == 1) {
			cout << "Run beta-FairTQ-Approx" << endl;
		}

		if (rskyband.size() == k) {
			cellTree* candidate = new cellTree(2 * (dim - 1));
			exactutk.push_back(candidate->root);
		} else {
			at = clock();

			float maxSP = -100;
			float minLoss = 100;
			float minDistance = 10000.0;

			float* bestF = new float[dim];
			int* bestTopK = new int[k];
			bool relaxed_found = false;

			sol->jaa(regions, exactutk, k, dim, PointSet, *rtree,  intersection, protectedAttribute, generalization, partitionPoint, SPinTheory,
				alpha, beta, numOfKind, originalF, maxSP, minDistance, minLoss, kindCount, num_levels, bestF, bestTopK, algorithm, &relaxed_found); // refinement
			treeSpacecost = sol->treeSize;

			//cout << "Regions: " << exactutk.size() << endl;
			resultSize += exactutk.size();
			ad = clock();

			cout << "============================================" << endl;
			cout << "Total time cost: " << fixed << (ad - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
			//cout << "Total space cost: " << fixed << totalSpaceCost + treeSpacecost << " MB " << endl;

			if (algorithm == 1) {
				if (relaxed_found) {
					unordered_set<long int> bestTopKSet;
					for (int ki = 0; ki < k; ki++) {
						bestTopKSet.insert(bestTopK[ki]);
					}

					maxSP = getAlphaFairness(bestTopKSet, PointSet, protectedAttribute, dim, generalization, partitionPoint, alpha, numOfKind, kindCount, num_levels);
				}
			}

			cout << "Best Utility Function:";
			for (int di = 0; di < dim; di++) {
				cout << " " << bestF[di]/U_RANGE;
			}
			cout << endl;
			cout << "Best Top-k:";
			for (int ki = 0; ki < k; ki++) {
				cout << " " << bestTopK[ki];
			}
			cout << endl;
			cout << "Best alpha-Fairness: " << maxSP << endl;
			cout << "Best Modification Penalty: " << minDistance / U_RANGE <<endl;
			// cout<<"minLoss: "<< minLoss <<endl;

			//Spacecost += treeSpacecost;
			cout << "============================================" << endl;

			string fn_time;
			string fn_stat;
			if (algorithm == 1) {
				fn_time = "time_ftq_appr1.txt";
				fn_stat = "stat_ftq_appr1.txt";
			} else if (algorithm == 2) {
				fn_time = "time_ftq_appr2.txt";
				fn_stat = "stat_ftq_appr2.txt";
			} else {
				fn_time = "time_ftq_exact.txt";
				fn_stat = "stat_ftq_exact.txt";
			}
			
			FILE* stream_time = fopen(fn_time.c_str(), "w");
			FILE* stream_stat = fopen(fn_stat.c_str(), "w");

			fprintf(stream_time, "%f\n", (ad - at)*1.0 / CLOCKS_PER_SEC);
			
			for (int di = 0; di < dim; di++) {
				fprintf(stream_stat, "%f ", bestF[di]/U_RANGE);
			}
			fprintf(stream_stat, "\n");
			for (int ki = 0; ki < k; ki++) {
				fprintf(stream_stat, "%d ", bestTopK[ki]);
			}
			fprintf(stream_stat, "\n");
			fprintf(stream_stat, "%f\n", maxSP);
			fprintf(stream_stat, "%f\n", minDistance/U_RANGE);
			// fprintf(stream_stat, "%f\n", minLoss);

			// // try to restore the top-k set under bestF
			// // however, I found that the top-k set come from rskyband, which is the intersection of skyband and klayers (in onionlayer)
			// vector<pair<float, long int>> tuple_pairs;
			// for (auto & ti: rskyband) {
			// // for (long int ti = 1; ti <= objCnt; ti++) {
			// 	float utility = 0.0;
			// 	for (int di = 0; di < dim ; di++) {
			// 		utility += bestF[di] * (PointSet[ti][di] + PointSet[ti][di + dim]) / 2.0;
			// 	}
			// 	tuple_pairs.push_back(make_pair(utility, ti));
			// }
			// sort(tuple_pairs.begin(), tuple_pairs.end(), [](pair<float, long int> const& a, pair<float, long int> const& b) -> bool {
			// 	return a.first > b.first;
			// });
			// vector<long int> best_TopKSet;
			// for(int ki = 0; ki < k; ki++)
			// {
			// 	best_TopKSet.push_back(tuple_pairs[ki].second);
			// }
			// cout << "best F:";
			// for (int di = 0; di < dim; di++) {
			// 	cout << " " << bestF[di]/U_RANGE;
			// }
			// cout << endl;
			// cout << "best top-k:";
			// for (int ki = 0; ki < k * 2; ki++) {
			// 	cout << " " << best_TopKSet[ki] << " (" << tuple_pairs[ki].first << ")";
			// }
			// cout << endl;


			delete[] bestF;
			delete[] bestTopK;

			fclose(stream_time);
			fclose(stream_stat);
		}
	}

	alle = clock();

	cout << endl;

	return 0;
}



// void initContainer(int& insertCount, vector<cell*>& finalResult, vector<long>& newAddSL, multimap<float, int>& heap, unordered_set<long int>& skylines, unordered_set<long int>& removeSL, unordered_set<long int>& singular,
// 	unordered_map<long int, long int>& RecordIDtoHalfPlaneID, vector<vector<float>>& HalfSpaces, unordered_set<long int>& ignorerecords)
// {
// 	insertCount = 0;
// 	for (int i = 0; i < finalResult.size(); i++)
// 	{
// 		finalResult[i]->release();
// 		delete finalResult[i];
// 	}
// 	finalResult.clear();
// 	vector<cell*>().swap(finalResult);
// 	newAddSL.clear();
// 	vector<long>().swap(newAddSL);
// 	heap.clear();
// 	multimap<float, int>().swap(heap);
// 	skylines.clear();
// 	unordered_set<long int>().swap(skylines);
// 	removeSL.clear();
// 	unordered_set<long int>().swap(removeSL);
// 	singular.clear();
// 	unordered_set<long int>().swap(singular);
// 	RecordIDtoHalfPlaneID.clear();
// 	unordered_map<long int, long int>().swap(RecordIDtoHalfPlaneID);
// 	HalfSpaces.clear();
// 	vector<vector<float>>().swap(HalfSpaces);
// 	ignorerecords.clear();
// 	unordered_set<long int>().swap(ignorerecords);
// }
