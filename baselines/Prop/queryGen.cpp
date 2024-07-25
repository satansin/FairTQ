#include "queryGen.h"

void QueryGenerator(vector < vector<float> >& QryGen, const int dim, const float sigma)
{
	QryGen.clear();
	vector<float> region(2 * dim, 0);
	float query[20];
	srand(0);
	float U_RANGE = 1000.0;
	while (QryGen.size() < 1)
	{
		float sum = 0;
		for (int di = 0; di < dim; di++)
		{
			query[di] = rand() *1.0 / RAND_MAX;
			sum += query[di];
		}

		if (sum < U_RANGE - sigma)
		{
			for (int di = 0; di < dim; di++)
			{
				region[di * 2] = query[di];
				region[di * 2 + 1] = query[di] + sigma*U_RANGE;
			}
			QryGen.push_back(region);
		}
	}
}

void GenerateData(int dim, int size)
{
	srand(NULL);
	vector<float> record(dim, 0);
	for (int i = 0; i < size; i++)
	{
		for (int di = 0; di < dim; di++)
		{
			record[di] = rand() *1.0 / RAND_MAX;
		}
		cout << i+1 << " ";
		for (int di = 0; di < dim; di++)
		{
			printf("%.4f ", record[di] - 0.0001);
		}
		for (int di = 0; di < dim; di++)
		{
			printf("%.4f ", record[di] + 0.0001);
		}
		cout << endl;
	}
	
}