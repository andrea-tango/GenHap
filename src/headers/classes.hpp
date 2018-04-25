#ifndef CLASSES_H_INCLUDED
#define CLASSES_H_INCLUDED

#include "libraries.hpp"

class Gene
{
public:
	Gene();
	void mutatePosition();
	int getValue() const;

private:
	void generatePosition();
	int value;
};


class Chromosome
{
public:
	Chromosome();

	// ~Chromosome();

	Chromosome(vector<vector<int> > &M, vector<vector<int> > &M_weight, vector<vector<int> > &listIndex);

	Chromosome(vector<vector<int> > &M, vector<vector<int> > &M_weight, vector<vector<int> > &listIndex,
	double mut_rate, Chromosome &parent1, Chromosome &parent2, double cross_pointIn, bool flip);

	double getFitness() const;
	vector<Gene> getGenes() const;
	int getCrossPoint() const;
	
	vector<int> getIndexC1() const;
	vector<int> getIndexC2() const;

	vector<int> getH1() const;
	vector<int> getH2() const;

	int getNum0() const;
	int getNum1() const;

	void calculateFitness(vector<vector<int> > &M, vector<vector<int> > &M_weight, vector<vector<int> > &listIndex);


	vector<Gene> genes;
	
private:

	double fitness;
	int cross_point;

	vector<int> indexC1;
	vector<int> indexC2;

	vector<int> h1;
	vector<int> h2;

	double errorFunction(int n, int num0, int num1, vector<vector<int> > &C1, vector<vector<int> > &C2, 
	vector<int> &indexC1, vector<int> &indexC2, vector<vector<int> > &M_weight, vector<vector<int> > &listIndex);

	void calculate_n0_n1(vector<vector<int> > &C, vector<vector<int> > &M_weight, vector<vector<int> > &listIndex, 
	vector<int> &indexC, vector<int> &N0, vector<int> &N1);

	vector<int> calculate_h(vector<int> &N0, vector<int> &N1);

	int distanceFragment(vector<int> &f1, vector<int> &f2);
	int singleDistance(int x, int y);
};


class GeneticOperation
{
public:
	void mutate(vector<Gene> &genes, double rate, bool flip);
	int crossover(Chromosome parent1, Chromosome parent2, vector<Gene> &genes, int cross_pointIn);
};


class GeneticAlgorithm
{
public:
	GeneticAlgorithm(int num_opt, int rank, bool verboseIn, bool savingIn, string settingsIn);
	void startGA(string pathOutput, vector<vector<int> > &MIn, vector<vector<int> > &M_weightIn,
		vector<vector<int> > &listIndexIn, int min, int max, int num_opt);

	void evolve(vector<Chromosome> &pop, int num_opt, int elitism, string method, int numberInd, Chromosome &best);

private:
	void initialize(vector<Chromosome> &pop, int num_opt, Chromosome &best);
	void evolve();
	void getPartitions(int m, int n, vector<Gene> &genes, vector<vector<int> > &M, vector<vector<int> > &C1, vector<vector<int> > &C2);
	void readParameters();
	bool verbose;
	bool saving;
	string OUTPUT_NAME_FIT;
	string OUTPUT_NAME_INFO;
	string OUTPUT_NAME_HAP;

	int POP_SIZE;
	double MUT_RATE;
	double CROSS_RATE;
	int GENERATIONS;
	int CHILDREN_PER_GEN;
	int elitism;
	int numInd;
	string selection;
	vector<vector<int> > M;
	vector<vector<int> > M_weight;
	vector<vector<int> > listIndex;

	string settings;

};


class Reader
{
public:
	// Reader();

	void readFromWif(string pathIn, string pathOut, int gamma, bool saveMatrix,
				vector<vector<int> > &matrixIn, vector<vector<int> > &matrixWeightIn,
				vector<int> &listIndexIn, vector<vector<int> > &listRunLenIn);
};


#endif //CLASSES_H_INCLUDED