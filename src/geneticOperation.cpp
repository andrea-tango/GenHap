#ifdef _WIN32
	#include "headers\\classes.hpp"
#elif _WIN64
	#include "headers\\classes.hpp"
#else
	#include "headers/classes.hpp"
#endif


// Mutation operator, implementing both a classic mutation and  bit-flipping mutation in which
// a random number of consecutive elements of the individual is mutated
void GeneticOperation::mutate(vector<Gene> &genes, double rate, bool flip)
{

	random_device rd;
	mt19937 gen{rd()};
	uniform_real_distribution<double> distributionReal(0.0,1.0);
 
	for(int i=0; i < genes.size(); i++)
	{
		if(distributionReal(gen) < rate)
		{
			if(flip)
			{
				if(i < (genes.size()-1))
				{
					uniform_int_distribution<int> distributionInt(i+1, genes.size()-1);
					int endFlip = distributionInt(gen);
					for(int j=i; j < endFlip; j++)
					{
						genes[j].mutatePosition();
					}
					break;
				}
				else
				{
					genes[i].mutatePosition();
				}
			}
			else
			{
				genes[i].mutatePosition();
			}
		}
	}
}

// Crossover operator, implementing a single-point crossover with mixing ratio equal to 0.5
int GeneticOperation::crossover(Chromosome parent1, Chromosome parent2, vector<Gene> &genes, int cross_point)
{
	random_device rd;
	mt19937 gen{rd()};

	vector<Gene> genesP1 = parent1.getGenes();
	vector<Gene> genesP2 = parent2.getGenes();

	int numberGenes = parent1.getGenes().size();
	int half = (int) round((numberGenes*1.0) / 2.0);

	if(cross_point != -1.0)
	{
		if(cross_point >= half)
		{
			for(int i=0; i < (cross_point-half); i++)
				genes.push_back(genesP2[i]);
			
			for(int i=cross_point-half; i < cross_point; i++)
				genes.push_back(genesP1[i]);
			
			for(int i=cross_point; i < numberGenes; i++)
				genes.push_back(genesP2[i]);			
		}
		else
		{
			for(int i=0; i < cross_point; i++)
				genes.push_back(genesP1[i]);
			
			for(int i=cross_point; i < (cross_point+half); i++)
				genes.push_back(genesP2[i]);
			
			for(int i=cross_point+half; i < numberGenes; i++)
				genes.push_back(genesP1[i]);
		}
		return cross_point;
	}	
	else
	{	
		uniform_int_distribution<int> distributionInt(0, numberGenes-1);
		int randNum = distributionInt(gen);

		if(randNum >= half)
		{
			for(int i=0; i < (randNum-half); i++)
				genes.push_back(genesP1[i]);
			
			for(int i=randNum-half; i < randNum; i++)
				genes.push_back(genesP2[i]);
			
			for(int i=randNum; i < numberGenes; i++)
				genes.push_back(genesP1[i]);
		}
		else
		{			
			for(int i=0; i < randNum; i++)
				genes.push_back(genesP2[i]);
			
			for(int i=randNum; i < (randNum+half); i++)
				genes.push_back(genesP1[i]);
			
			for(int i=randNum+half; i < numberGenes; i++)
				genes.push_back(genesP2[i]);
		}
		return randNum;
	}
}