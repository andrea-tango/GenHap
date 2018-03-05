#ifdef _WIN32
	#include "headers\\classes.hpp"
	const string slash = "\\";
	const int os = 0;
#elif _WIN64
	#include "headers\\classes.hpp"
	const string slash = "\\";
	const int os = 0;
#else
	#include "headers/classes.hpp"
	const string slash = "/";
	const int os = 1;
#endif

GeneticAlgorithm::GeneticAlgorithm(int num_opt, int rank, bool verboseIn, bool savingIn, string settingsIn):verbose(verboseIn), saving(savingIn), settings(settingsIn)
{
	if(verbose)
		cout << "* Start optimization n. " << num_opt <<  " on rank: " << rank << "\n";
}

// Reading GA parameters (such as, number of generations, crossover rate ...)
void GeneticAlgorithm::readParameters()
{

	// default parameters
	GENERATIONS = 100;
	POP_SIZE    = 100;
	CROSS_RATE  = 0.9;
	MUT_RATE    = 0.05;
	elitism     = 1;
	selection   = "tournament";
	numInd      = (int) round(POP_SIZE*0.1);

	vector<vector<string> > parameters;
	
	try
	{
		string line;
		ifstream infile(settings);
		while (getline(infile, line))
		{
			istringstream buf(line);
			vector<string> sublista;
			for(string token; getline(buf, token, '='); )
			{
				sublista.push_back(token);
			}
			parameters.push_back(sublista);
		}
	}
	catch(const std::exception& e){return;}

	for(int i=0; i < parameters.size(); i++)
	{
		try
		{
			if(strcmp(parameters[i][0].c_str(), "Generations") == 0)
			{
				GENERATIONS = (int) strtod(parameters[i][1].c_str(), NULL);
			}
		}
		catch(const std::exception& e){continue;}
		
		try
		{
			if(strcmp(parameters[i][0].c_str(), "PopSize") == 0)
			{
				POP_SIZE = (int) strtod(parameters[i][1].c_str(), NULL);
			}
		}
		catch(const std::exception& e){continue;}

		try
		{
			if(strcmp(parameters[i][0].c_str(), "Selection") == 0)
			{
				selection = parameters[i][1];
				if(strcmp(selection.c_str(), "tournament") == 0)
				{
					numInd = (int) round(POP_SIZE*0.1);
				}
			}
		}
		catch(const std::exception& e){continue;}

		try
		{
			if(strcmp(parameters[i][0].c_str(), "Mr") == 0)
			{
				MUT_RATE = strtod(parameters[i][1].c_str(), NULL);
			}
		}
		catch(const std::exception& e){continue;}

		try
		{
			if(strcmp(parameters[i][0].c_str(), "Cr") == 0)
			{
				CROSS_RATE = strtod(parameters[i][1].c_str(), NULL);
			}
		}
		catch(const std::exception& e){continue;}

		try
		{
			if(strcmp(parameters[i][0].c_str(), "Elitism") == 0)
			{
				elitism = strtod(parameters[i][1].c_str(), NULL);
			}
		}
		catch(const std::exception& e){continue;}
	}

	CHILDREN_PER_GEN = POP_SIZE - elitism;
}

void GeneticAlgorithm::startGA(string pathOutput, vector<vector<int> > &MIn, vector<vector<int> > &M_weightIn,
		vector<vector<int> > &listIndexIn, int minimo, int massimo, int num_opt)
{

	readParameters();
	M = MIn;
	M_weight = M_weightIn;
	listIndex = listIndexIn;
	FILE * fl;

	string folder1 = pathOutput + slash + "output" + slash + "optimization" + to_string(num_opt);
	struct stat sb;
	if( !(stat(folder1.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)))
	{
		if (os == 0)
		{
			folder1 = "MD " + folder1;
			int sis = system(folder1.c_str());
		}
		else
		{
			folder1 = "mkdir " + folder1;
			int sis = system(folder1.c_str());
		}
	}

	OUTPUT_NAME_INFO = pathOutput + slash + "output" + slash + "optimization" + to_string(num_opt) + slash + "informations";
	OUTPUT_NAME_FIT  = pathOutput + slash + "output" + slash + "optimization" + to_string(num_opt) + slash + "fitness";
	OUTPUT_NAME_HAP  = pathOutput + slash + "output" + slash + "optimization" + to_string(num_opt) + slash;

// ********************************* Writing Genetic Algorithm settings *********************************
	
	if(saving)
	{
		fl=fopen(OUTPUT_NAME_INFO.c_str(), "w");

		if( fl==NULL )
		{
			printf("Error open output file: %s\n", OUTPUT_NAME_INFO.c_str());
			exit(-9);
		}
		fprintf(fl, "%s", "******************************************************\n");
		fprintf(fl, "%s","\t\t\t GA INFORMATION\n");
		fprintf(fl, "%s %d%s", "Number of chromosome:", POP_SIZE, "\n");
		fprintf(fl, "%s %d%s", "Number of elite chromosomes:", elitism, "\n");
		fprintf(fl, "%s %d%s", "Number of genes:", (int) M.size(), "\n");
		fprintf(fl, "%s %d%s", "Number of generations:", GENERATIONS, "\n");
		fprintf(fl, "%s %.2f%s", "Crossover rate:", CROSS_RATE, "\n");
		fprintf(fl, "%s %.2f%s", "Mutation rate:", MUT_RATE, "\n");

		fclose(fl);
	}

// ********************************* Genetic Algorithm execution *********************************
	vector<Chromosome> pop;

	Chromosome best;
	initialize(pop, num_opt, best);

	evolve(pop, num_opt, elitism, selection, numInd, best);

// ********************************* Saving results *********************************

	vector<int> indexC1;
	vector<int> indexC2;

	vector<vector<int> > C1;
	vector<vector<int> > C2;

	vector<Gene> genes = best.getGenes();

	getPartitions(M.size(), M[0].size(), genes, M, C1, C2);

	if(saving)
	{
		string s = OUTPUT_NAME_HAP + slash + "partition0";
		fl=fopen(s.c_str(), "w");

		if( fl==NULL )
		{
			printf("Error open output file: %s\n", OUTPUT_NAME_INFO.c_str());
			exit(-9);
		}

		for(int i=0; i < C1.size(); i++)
		{	
			for(int j=0; j < C1[i].size(); j++)
			{
				if(j < C1[i].size()-1)
				{
					if(C1[i][j]==-1)
						fprintf(fl, "%s ", "-");
					else
						fprintf(fl, "%d ", C1[i][j]);
				}
				else
				{
					if(C1[i][j]==-1)
						fprintf(fl, "%s\n", "-");
					else
						fprintf(fl, "%d\n", C1[i][j]);
				}
			}
		}

		fclose(fl);

		s = OUTPUT_NAME_HAP + slash + "partition1";
		fl=fopen(s.c_str(), "w");

		if( fl==NULL )
		{
			printf("Error open output file: %s\n", OUTPUT_NAME_INFO.c_str());
			exit(-9);
		}

		for(int i=0; i < C2.size(); i++)
		{	
			for(int j=0; j < C2[i].size(); j++)
			{
				if(j < C2[i].size()-1)
				{
					if(C2[i][j]==-1)
						fprintf(fl, "%s ", "-");
					else
						fprintf(fl, "%d ", C2[i][j]);
				}
				else
				{
					if(C2[i][j]==-1)
						fprintf(fl, "%s\n", "-");
					else
						fprintf(fl, "%d\n", C2[i][j]);
				}
			}
		}

		fclose(fl);
	}

	string path = OUTPUT_NAME_HAP + slash + "haplotypes";
	fl=fopen(path.c_str(), "w");

	if( fl==NULL )
	{
		printf("Error open output file: %s\n", OUTPUT_NAME_INFO.c_str());
		exit(-9);
	}

	vector<int> h1 = best.getH1();
	vector<int> h2 = best.getH2();

	for(int j=0; j < h1.size(); j++)
	{
		if(h1[j]==0 || h1[j]==1)
		{
			int num0 = 0;
			int num1 = 0;
			for(int i=0; i < C1.size(); i++)
			{
				if(C1[i][j] == 0)
				{
					indexC1 = best.getIndexC1();
					int idx = indexC1[i];
					num0 += M_weight[idx][j];
				}
				else
				{
					if(C1[i][j] == 1)
					{
						indexC1 = best.getIndexC1();
						int idx = indexC1[i];
						num1 += M_weight[idx][j];
					}
				}
			}
			if(num0 > num1)
				h1[j] = 0;
			else if(num1 > num0)
				h1[j] = 1;
		}
	}

	for(int j=0; j < h2.size(); j++)
	{
		if(h2[j]==0 || h2[j]==1)
		{
			int num0 = 0;
			int num1 = 0;
			for(int i=0; i < C2.size(); i++)
			{
				int num0 = 0;
				int num1 = 0;
				if(C2[i][j] == 0)
				{
					indexC2 = best.getIndexC2();
					int idx = indexC2[i];
					num0 += M_weight[idx][j];
				}
				else
				{
					if(C2[i][j] == 1)
					{
						indexC2 = best.getIndexC2();
						int idx = indexC2[i];
						num1 += M_weight[idx][j];
					}
				}
			}
			if(num0 > num1)
				h2[j] = 0;
			else if(num1 > num0)
				h2[j] = 1;
		}
	}

	for(int i=0; i < h1.size(); i++)
	{
		if(h1[i]==0)
		{
			int count=0;
			int countNoHole = 0;

			for(int j=0; j < C1.size(); j++)
			{	
				if(C1[j][i]!=-1)
					countNoHole++;

				if(C1[j][i]==1)
				{
					count++;
				}
			}
			if(countNoHole==count && count!=0)
				h1[i]=1;
		}
		else if(h1[i]==1)
		{
			int count=0;
			int countNoHole = 0;

			for(int j=0; j < C1.size(); j++)
			{	
				if(C1[j][i]!=-1)
					countNoHole++;

				if(C1[j][i]==0)
				{
					count++;
				}
			}
			if(countNoHole==count && count!=0)
				h1[i]=0;
		}
	}

	for(int i=0; i < h2.size(); i++)
	{
		if(h2[i]==0)
		{
			int count=0;
			int countNoHole = 0;

			for(int j=0; j < C2.size(); j++)
			{	
				if(C2[j][i]!=-1)
					countNoHole++;

				if(C2[j][i]==1)
				{
					count++;
				}
			}
			if(countNoHole==count && count!=0)
				h2[i]=1;
		}
		else if(h2[i]==1)
		{
			int count=0;
			int countNoHole = 0;

			for(int j=0; j < C2.size(); j++)
			{	
				if(C2[j][i]!=-1)
					countNoHole++;

				if(C2[j][i]==0)
				{
					count++;
				}
			}
			if(countNoHole==count && count!=0)
				h2[i]=0;
		}
	}

	if(C1.size() > 0 & C2.size() > 0)
	{
		for(int i=0; i < minimo; i++)
			fprintf(fl, "%s", "-");
		
		for(int i=0; i < h1.size(); i++)
		{
			if(h1[i] == -1)
				fprintf(fl, "%s", "-");
			else
				fprintf(fl, "%d", h1[i]);
		}

		for(int i=h1.size()+minimo; i < massimo; i++)
			fprintf(fl, "%s", "-");
		fprintf(fl, "%s", "\n");

		for(int i=0; i < minimo; i++)
			fprintf(fl, "%s", "-");

		for(int i=0; i < h2.size(); i++)
		{
			if(h2[i] == -1)
				fprintf(fl, "%s", "-");
			else
				fprintf(fl, "%d", h2[i]);
		}

		for(int i=h2.size()+minimo; i < massimo; i++)
			fprintf(fl, "%s", "-");
	}
	else
	{
		if(C1.size()>0)
		{
			for(int i=0; i < minimo; i++)
				fprintf(fl, "%s", "-");

			for(int i=0; i < h1.size(); i++)
			{
				if(h1[i] == -1)
					fprintf(fl, "%s", "-");
				else
					fprintf(fl, "%d", h1[i]);
			}

			for(int i=h1.size()+minimo; i < massimo; i++)
				fprintf(fl, "%s", "-");
			fprintf(fl, "%s", "\n");

			for(int i=0; i < massimo; i++)
				fprintf(fl, "%s", "-");
		}

		else
		{
			for(int i=0; i < massimo; i++)
				fprintf(fl, "%s", "-");

			fprintf(fl, "%s", "\n");

			for(int i=0; i < minimo; i++)
				fprintf(fl, "%s", "-");
			
			for(int i=0; i < h2.size(); i++)
			{
				if(h2[i] == -1)
					fprintf(fl, "%s", "-");
				else
					fprintf(fl, "%d", h2[i]);
			}

			for(int i=h2.size()+minimo; i < massimo; i++)
				fprintf(fl, "%s", "-");
		}
	}

	fclose(fl);
}

void GeneticAlgorithm::getPartitions(int m, int n, vector<Gene> &genes, vector<vector<int> > &M, vector<vector<int> > &C1, vector<vector<int> > &C2)
{
	int num0 = 0;
	int num1 = 0;

	for(int i=0; i < m; i++)
	{
		if(genes[i].getValue() == 0)
			num0++;
		else
			num1++;
	}

	C1 = vector<vector<int> >(num0, vector<int>(n));
	C2 = vector<vector<int> >(num1, vector<int>(n));
	
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
		}
		else
		{
			for(int j=0; j < n; j++)
			{
				C2[count2][j] = M[i][j];
			}
			count2++;
		}
	}
}

// Initialization of the population
void GeneticAlgorithm::initialize(vector<Chromosome> &pop, int num_opt, Chromosome &best)
{

	for(int i=0; i < POP_SIZE; i++)
	{
		Chromosome c = Chromosome(M, M_weight, listIndex);
		pop.push_back(c);
	}

	// riordino della popolazione in base alla fitness
	sort(pop.begin(), pop.end(),
		[]( const Chromosome &lhs, const Chromosome &rhs )
		{
			return lhs.getFitness() < rhs.getFitness();
		}
	);

	best = pop[0];

	if(verbose)
		cout << "Best fitness at generation 0 of optimization n. " << num_opt << " = " << best.getFitness() << "\n";

	FILE * fl;
	fl=fopen(OUTPUT_NAME_FIT.c_str(), "w");

	if( fl==NULL )
	{
		printf("Error open output file: %s\n", OUTPUT_NAME_INFO.c_str());
		exit(-9);
	}

	fprintf(fl, "%.2f%s", best.getFitness(), "\n");

	fclose(fl);

}

// Evolution of the population
void GeneticAlgorithm::evolve(vector<Chromosome> &pop, int num_opt, int elitism, string method, int numberInd, Chromosome &best)
{
	random_device rd;
	mt19937 gen{rd()};
	uniform_int_distribution<int> distributionInt(0, POP_SIZE-1);
	uniform_real_distribution<double> distributionReal(0.0,1.0);

	FILE * fl;
	
	if(saving)
	{
		fl=fopen(OUTPUT_NAME_INFO.c_str(), "a");

		if( fl==NULL )
		{
			printf("Error open output file: %s\n", OUTPUT_NAME_INFO.c_str());
			exit(-9);
		}

		if(strcmp(method.c_str(), "wheel") == 0)
		{
			fprintf(fl, "%s", "Selection: wheel roulette selection\n");
		}
		else
		{
			if(strcmp(method.c_str(), "rank") == 0)
			{
				fprintf(fl, "%s", "Selection: rank selection \n");
			}
			else
			{
				fprintf(fl, "%s %d %s", "Selection: tournament selection with", numberInd, "individuals\n");
			}
		}
		fclose(fl);
	}

	double oldFitness = best.getFitness();
	bool flippa = false;
	int counterEqualFitness = 0;

	Chromosome parent_1;
	Chromosome parent_2;

	GeneticOperation op;

	bool noImprovement = false;
	int prematureExit = (int) round(GENERATIONS*0.25);

	for(int i=0; i < GENERATIONS; i++)
	{

		if(counterEqualFitness == prematureExit-1)
		{
			noImprovement = true;
		}

		if(best.getFitness() == 0.0 or noImprovement)
		{
			if(verbose)
				cout << "Best fitness at generation " << i << " of optimization n. " << num_opt << " = " << best.getFitness() << "\n";
			
			break;
		}
		vector<double> probabilities;

		// roulette wheel selection
		if(strcmp(method.c_str(), "wheel") == 0)
		{
			
			double sum_fit = 0;
			for(int j=0; j < POP_SIZE; j++)
			{
				probabilities.push_back(pop[j].getFitness());
				sum_fit += pop[j].getFitness();
			}

			double sumProb = 0;
			for(int j=0; j < POP_SIZE; j++)
			{
				probabilities[j] = probabilities[j] / sum_fit;
				probabilities[j] = 1 - probabilities[j];
				sumProb += probabilities[j];
			}
			for(int j=0; j < POP_SIZE; j++)
			{
				probabilities[j] = probabilities[j] / sumProb;
			}
		}
		// ranking selection
		else
		{
			if(strcmp(method.c_str(), "rank") == 0)
			{
				vector<double> probabilities;
				double sumRank = 0;
				
				for(int j=0; j < POP_SIZE; j++)
				{
					probabilities.push_back(POP_SIZE-j);
					sumRank += (POP_SIZE-j);
				}
				for(int j=0; j < POP_SIZE; j++)
				{
					probabilities[j] = probabilities[j] / sumRank;
				}
			}
		}

		int countWhile = CHILDREN_PER_GEN;
		vector<Chromosome> children;

		if(counterEqualFitness >= 2)
			flippa = true;
		
		while(countWhile > 0)
		{
			// tournament selection
			if(strcmp(method.c_str(), "tournament") == 0)
			{	
				vector<int> dist1;
				vector<int> dist2;
				for(int j=0; j < numberInd; j++)
				{
					dist1.push_back(distributionInt(gen));
					dist2.push_back(distributionInt(gen));
				}

				vector<Chromosome> individuals1;
				vector<Chromosome> individuals2;

				for(int k=0; k < numberInd; k++)
				{
					individuals1.push_back(pop[dist1[k]]);
					individuals2.push_back(pop[dist2[k]]);
				}

				// sorting the individuals
				sort(individuals1.begin(), individuals1.end(),
					[]( const Chromosome &lhs, const Chromosome &rhs )
					{
						return lhs.getFitness() < rhs.getFitness();
					}
				);

				sort(individuals2.begin(), individuals2.end(),
					[]( const Chromosome &lhs, const Chromosome &rhs )
					{
						return lhs.getFitness() < rhs.getFitness();
					}
				);

				parent_1 = individuals1[0];
				parent_2 = individuals2[0];
			}

			else
			{
				discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
				int value1 = distribution(gen);
				int value2 = distribution(gen);
				parent_1 = pop[value1];
				parent_2 = pop[value2];				
			}

			if(countWhile == 1)
			{
				Chromosome child0 = Chromosome(M, M_weight, listIndex, MUT_RATE, parent_1, parent_2, -1.0, flippa);
				Chromosome child1 = Chromosome(M, M_weight, listIndex, MUT_RATE, parent_1, parent_2, child0.getCrossPoint(), flippa);

				if(child0.getFitness() < child1.getFitness())
					children.push_back(child0);
				else
					children.push_back(child1);

				countWhile -= 1;
			}

			else
			{
				if(distributionReal(gen) < CROSS_RATE)
				{
					Chromosome child0 = Chromosome(M, M_weight, listIndex, MUT_RATE, parent_1, parent_2, -1.0, flippa);
					Chromosome child1 = Chromosome(M, M_weight, listIndex, MUT_RATE, parent_1, parent_2, child0.getCrossPoint(), flippa);
					children.push_back(child0);
					children.push_back(child1);
					countWhile -= 2;
				}
				else
				{
					op.mutate(parent_1.genes, MUT_RATE, flippa);
					op.mutate(parent_2.genes, MUT_RATE, flippa);

					parent_1.calculateFitness(M, M_weight, listIndex);
					parent_2.calculateFitness(M, M_weight, listIndex);

					children.push_back(parent_1);
					children.push_back(parent_2);
					countWhile -= 2;
				}
			} 		
		}

		for(int j=elitism; j < POP_SIZE; j++)
		{
			pop[j] = children[j-elitism];
		}

		sort(pop.begin(), pop.end(),
			[]( const Chromosome &lhs, const Chromosome &rhs )
			{
				return lhs.getFitness() < rhs.getFitness();
			}
		);

		if(pop[0].getFitness() < best.getFitness())
		{
			best = pop[0];
		}

		if(best.getFitness() == oldFitness)
		{
			counterEqualFitness += 1;
		}
		else
		{
			counterEqualFitness = 0;
			oldFitness = best.getFitness();
			flippa = false;
		}

		if(verbose)
			cout << "Best fitness at generation " << i << " of optimization n. " << num_opt << " = " << best.getFitness() << "\n";

		if(verbose & (i == GENERATIONS - 1))
			cout << "Best fitness at last generation of optimization n. " << num_opt << " = " << best.getFitness() << "\n";
	
		fl=fopen(OUTPUT_NAME_FIT.c_str(), "a");

		if( fl==NULL )
		{
			printf("Error open output file: %s\n", OUTPUT_NAME_INFO.c_str());
			exit(-9);
		}

		fprintf(fl, "%.2f%s", best.getFitness(), "\n");

		fclose(fl);

	}
}