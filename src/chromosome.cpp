#ifdef _WIN32
	#include "headers\\classes.hpp"
#elif _WIN64
	#include "headers\\classes.hpp"
#else
	#include "headers/classes.hpp"
#endif

Chromosome::Chromosome()
{

}

Chromosome::Chromosome(vector<vector<int> > &M, vector<vector<int> > &M_weight, vector<vector<int> > &listIndex)
{
	for(int i = 0; i < M.size(); i++)
	{
		Gene gene = Gene();
		genes.push_back(gene);
	}
	calculateFitness(M, M_weight, listIndex);
}

Chromosome::Chromosome(vector<vector<int> > &M, vector<vector<int> > &M_weight, vector<vector<int> > &listIndex,
	double mut_rate, Chromosome &parent1, Chromosome &parent2, double cross_pointIn, bool flip)
{
	GeneticOperation op;
	cross_point = op.crossover(parent1, parent2, genes, cross_pointIn);
	op.mutate(genes, mut_rate, flip);
	calculateFitness(M, M_weight, listIndex);
}

// Fitness function
void Chromosome::calculateFitness(vector<vector<int> > &M, vector<vector<int> > &M_weight, vector<vector<int> > &listIndex)
{
	int m = M.size();
	int n = M[0].size();

	int num0 = 0;
	int num1 = 0;

	for(int i=0; i < m; i++)
	{
		if(genes[i].getValue() == 0)
			num0++;
		else
			num1++;
	}

	vector<vector<int> > C1(num0, vector<int>(n));
	vector<vector<int> > C2(num1, vector<int>(n));
	
	int count1 = 0;
	int	count2 = 0;
	for(int i=0; i < m; i++)
	{
		if(genes[i].getValue() == 0)
		{
			for(int j=0; j < n; j++)
			{
				C1[count1][j] = M[i][j];
			}
			count1++;
			indexC1.push_back(i);
		}
		else
		{
			for(int j=0; j < n; j++)
			{
				C2[count2][j] = M[i][j];
			}
			count2++;
			indexC2.push_back(i);
		}
	}
	fitness = errorFunction(n, num0, num1, C1, C2, indexC1, indexC2, M_weight, listIndex);
}

double Chromosome::errorFunction(int n, int num0, int num1, vector<vector<int> > &C1, vector<vector<int> > &C2, 
	vector<int> &indexC1, vector<int> &indexC2, vector<vector<int> > &M_weight, vector<vector<int> > &listIndex)
{

	double error1 = 0;
	double error2 = 0;

	h1 = vector<int>(n);
	h2 = vector<int>(n);

	vector<int> N0_1(n);
	vector<int> N1_1(n);

	vector<int> N0_2(n);
	vector<int> N1_2(n);

	if(num0)
	{
		calculate_n0_n1(C1, M_weight, listIndex, indexC1, N0_1, N1_1);
		h1 = calculate_h(N0_1, N1_1);
	}

	if(num1)
	{
		calculate_n0_n1(C2, M_weight, listIndex, indexC2, N0_2, N1_2);
		h2 = calculate_h(N0_2, N1_2);
	}

	if( (num0 > 0) & (num1 > 0) )
	{
		for(int j=0; j < h1.size(); j++)
		{
			if( (h1[j] == h2[j]) & (h1[j] == -1))
			{
				continue;
			}
			else
			{
				if(h1[j] == h2[j])
				{
					if(h1[j] == 0)
					{
						if(N0_1[j] > N0_2[j])
							h2[j] = 1;
						else
							h1[j] = 1;
					}
					else
					{
						if(N1_1[j] > N1_2[j])
							h2[j] = 0;
						else
							h1[j] = 0;
					}
				}
				else
				{
					// no ambigue positions
					if(h1[j] == -1 & h2[j]==0)
						h1[j] = 1;
					else if(h1[j] == -1 & h2[j]==1)
						h1[j] = 0;
					else if(h2[j] == -1 & h1[j]==0)
						h2[j] = 1;
					else if(h2[j] == -1 & h1[j]==1)
						h2[j] = 0;

				}
			}
		}
	}
	
	if(num0)
	{
		for(int j=0; j < C1.size(); j++)
		{
			error1 += distanceFragment(listIndex, C1[j], h1, indexC1[j]);
		}
	}
	if(num1)
	{
		for(int j=0; j < C2.size(); j++)
		{
			error2 += distanceFragment(listIndex, C2[j], h2, indexC2[j]);
		}
	}

	return error1+error2;
}

void Chromosome::calculate_n0_n1(vector<vector<int> > &C, vector<vector<int> > &M_weight, vector<vector<int> > &listIndex, 
	vector<int> &indexC, vector<int> &N0, vector<int> &N1)
{

	int m = C.size();
	int idx = 0;

	for(int i=0; i < m; i++)
	{
		for(int j=0; j < listIndex[indexC[i]][1]; j++)
		{
			if(C[i][j] == 0)
			{
				idx = indexC[i];
				N0[j] += M_weight[idx][j];
			}
			else
			{
				idx = indexC[i];
				N1[j] += M_weight[idx][j];			
			}
		}
	}
}

vector<int> Chromosome::calculate_h(vector<int> &N0, vector<int> &N1)
{
	int n = N0.size();
	vector<int> h(n);

	for(int j=0; j < n; j++)
	{

		if( (N0[j] == 0) & (N1[j] == 0) )
		{
			h[j] = -1;
		}
		else
		{
			if(N0[j] > N1[j])
			{
				h[j] = 0;
			}
			else
			{
				h[j] = 1;
			}
		}
	}

	return h;
}

int Chromosome::distanceFragment(vector<vector<int> > &listIndex, vector<int> &f1, vector<int> &f2, int j)
{
	int sumErr = 0;

	for(int i=0; i < f1.size(); i++)
	{
		sumErr += singleDistance(f1[i], f2[i]);
	}
	return sumErr;
}

int Chromosome::singleDistance(int x, int y)
{
	if( (x == -1) & (y == -1) )
	{
		return 0;
	}
	else
	{
		if( (x != -1) & (y != -1) & (x != y) )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
}

// *************************************************** Metodi Getter ***************************************************
double Chromosome::getFitness() const
{
	return fitness;
}

int Chromosome::getCrossPoint() const
{
	return cross_point;
}

vector<Gene> Chromosome::getGenes() const
{
	return genes;
}

vector<int> Chromosome::getIndexC1() const
{
	return indexC1;
}

vector<int> Chromosome::getIndexC2() const
{
	return indexC2;
}

vector<int> Chromosome::getH1() const
{
	return h1;
}

vector<int> Chromosome::getH2() const
{
	return h2;
}