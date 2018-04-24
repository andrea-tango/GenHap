/**
 *
 *                              GenHap
 * A Novel Computational Method Based on Genetic Algorithms for Haplotype Assembly
 *
 * Copyright (C) 2018  Andrea Tangherloni
 *
 * Distributed under the terms of the GNU General Public License (GPL)
 *
 * This file is part of GenHap.
 *
 * HapCol is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License v3.0 as published by
 * the Free Software Foundation.
 * 
 * GenHap is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the
 * GNU General Public License for more details.
 *
 **/

#ifdef _WIN32
	#include "headers\\libraries.hpp"
	#include "headers\\classes.hpp"
	const string slash = "\\";
	const int os = 0;
#elif _WIN64
	#include "headers\\libraries.hpp"
	#include "headers\\classes.hpp"
	const string slash = "\\";
	const int os = 0;
#else
	#include "headers/libraries.hpp"
	#include "headers/classes.hpp"
	const string slash = "/";
	const int os = 1;
#endif


const int WORKTAG = 0;
const int DIETAG  = 1;

const int END_S  = 11;


int slave(int rank, string pathOutput, bool verbose, bool saving, string settings);
void master(vector<vector<int> > &matrix, vector<vector<int> > &matrixWeight,
	vector<vector<int> > &listRunLen, vector<int> &listIndex, int size);

void findMinMax(int start, int end, vector<vector<int> > &v, int &min, int &max);
void saveHaplotypes(string pathOut, string name, int num, vector<vector<int> > &matrix, bool ambiguous);
int calculateSimilarity(string oldSubHap, vector<string> newSubHap, int length);

void saveHaplotypes(string pathOut, string name, int num, vector<vector<int> > &matrix, bool ambiguous)
{
	vector<vector<string> > subHaplotypes;

	string line, path;
	int length;
	
	for(int i=0; i < num; i++)
	{
		path = pathOut+slash+"output"+slash+"optimization"+to_string(i)+slash+"haplotypes";
		ifstream infile(path.c_str());
		vector<string> sublista;
		while (getline(infile, line))
		{
			istringstream buf(line);
			for(string token; getline(buf, token, ' '); )
			{
				sublista.push_back(token);
			}
		}
		subHaplotypes.push_back(sublista);
	}

	vector<string> haplotypes;

	haplotypes.push_back(subHaplotypes[0][0]);
	haplotypes.push_back(subHaplotypes[0][1]);
	
	for(int i=0; i < num-1; i++)
	{
		vector<string> oldSubHap;
		vector<string> newSubHap;
		
		length = subHaplotypes[i][0].length();	
		
		oldSubHap.push_back(subHaplotypes[i][0]);
		oldSubHap.push_back(subHaplotypes[i][1]);
		newSubHap.push_back(subHaplotypes[i+1][0]);
		newSubHap.push_back(subHaplotypes[i+1][1]);

		vector<int> indexOld;
		vector<int> indexNew;

		for(int h=0; h < 2; h++)
		{
			int indexN = calculateSimilarity(oldSubHap[h], newSubHap, length);
			int indexO = calculateSimilarity(haplotypes[h], oldSubHap, length);

			indexOld.push_back(indexO);
			indexNew.push_back(indexN);
		}

		if(indexNew[0] == indexNew[1] && i!=num-1)
		{
			indexNew.clear();
			vector<string> tmp;
			tmp.push_back(subHaplotypes[i+1][0]);
			tmp.push_back(subHaplotypes[i+1][1]);
			for(int h=0; h < 2; h++)
			{
				int indexN = calculateSimilarity(newSubHap[h], tmp, length);
				indexNew.push_back(indexN);
			}
		}

		for(int h=0; h < 2; h++)
		{
			int idxN = indexNew[h];
			int idxO = indexOld[h];

			string s = "";
			for(int j=0; j < length; j++)
			{
				if(newSubHap[idxN][j] == '-')
					s += haplotypes[idxO][j];
				else
					s += newSubHap[idxN][j];
			}

			haplotypes[idxO] = s;			
		}	
	}


	for(int i=0; i < haplotypes[0].size(); i++)
	{	
		if (haplotypes[0][i]=='1')
		{
			bool trovato = false;
			for(int j=0; j < matrix.size(); j++)
			{
				if(matrix[j][i]==1)
				{
					trovato = true;
					break;
				}
			}
			if(!trovato)
			{
				if(!ambiguous)
					haplotypes[0][i]='0';
				else
					haplotypes[0][i]='X';
			}
		}
		else if (haplotypes[0][i]=='0')
		{
			bool trovato = false;
			for(int j=0; j < matrix.size(); j++)
			{
				if(matrix[j][i]==0)
				{
					trovato = true;
					break;
				}
			}
			if(!trovato)
			{
				if(!ambiguous)
					haplotypes[0][i]='1';
				else
					haplotypes[0][i]='X';
			}
		}
	}

	for(int i=0; i < haplotypes[1].size(); i++)
	{	
		if (haplotypes[1][i]=='1')
		{
			bool trovato = false;
			for(int j=0; j < matrix.size(); j++)
			{
				if(matrix[j][i]==1)
				{
					trovato = true;
					break;
				}
			}
			if(!trovato)
			{
				if(!ambiguous)
					haplotypes[1][i]='0';
				else
					haplotypes[1][i]='X';
			}
		}
		else if (haplotypes[1][i]=='0')
		{
			bool trovato = false;
			for(int j=0; j < matrix.size(); j++)
			{
				if(matrix[j][i]==0)
				{
					trovato = true;
					break;
				}
			}
			if(!trovato)
			{
				if(!ambiguous)
					haplotypes[1][i]='1';
				else
					haplotypes[1][i]='X';
			}
		}
	}

	FILE * fl;
	string file = pathOut+slash+name;
	fl=fopen(file.c_str(), "w");

	if( fl==NULL )
	{
		printf("Error open output file: %s\n", file.c_str());
		exit(-9);
	}

	fprintf(fl, "%s\n", haplotypes[0].c_str());
	fprintf(fl, "%s",   haplotypes[1].c_str());
	fclose(fl);	
}

int calculateSimilarity(string oldSubHap, vector<string> newSubHap, int length)
{

	int count1 = 0;
	int count2 = 0;

	for(int j=0; j < length; j++)
	{
		if(newSubHap[0][j] != '-' && oldSubHap[j] == newSubHap[0][j])
			count1 += 1;
		if(newSubHap[1][j] != '-' && oldSubHap[j] == newSubHap[1][j])
			count2 += 1;
	}

	if(count1 > count2)
		return 0;
	else
	 	return 1;
}

void findMinMax(int start, int end, vector<vector<int> > &v, int &min, int &max)
{
	min = INT_MAX;
	max = INT_MIN;
	
	for(int i=start; i < end; i++)
	{
		if(v[i][0] <= min)
			min = v[i][0];

		if(v[i][1] >= max)
			max = v[i][1];
	}
}

// Master process that orchestrates the slave processes
void master(vector<vector<int> > &matrix, vector<vector<int> > &matrixWeight,
	vector<vector<int> > &listRunLen, vector<int> &listIndex, int size)
{
	MPI_Status status;
	int min, max, numRow, numCol;

	int shape = matrix[0].size();

	vector<int> v;

	if(listIndex.size() < size)
	{
		for(int i=1; i < listIndex.size(); i++)
		{
			int start = listIndex[i-1];
			int end   = listIndex[i];

			numRow = end - start;
			findMinMax(start, end, listRunLen, min, max);
			numCol = max - min;

			int idx = i-1;
			MPI_Send(&idx, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
			MPI_Send(&numRow, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
			MPI_Send(&numCol, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);

			for(int j=start; j < end; j++)
			{
				vector<int> v(numCol);
				for(int k=min; k < max; k++)
				{
					v[k-min] = matrix[j][k];
				}

				if(j != end-1)
					MPI_Send(&v[0], v.size(), MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
				else
					MPI_Send(&v[0], v.size(), MPI_INT, i, END_S, MPI_COMM_WORLD);
			}

			for(int j=start; j < end; j++)
			{
				vector<int> v(numCol);
				for(int k=min; k < max; k++)
				{
					v[k-min] = matrixWeight[j][k];
				}

				if(j != end-1)
					MPI_Send(&v[0], v.size(), MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
				else
					MPI_Send(&v[0], v.size(), MPI_INT, i, END_S, MPI_COMM_WORLD);
			}

			for(int j=start; j < end; j++)
			{
				vector<int> v(2);
				v[0] = listRunLen[j][0] - min;
				v[1] = listRunLen[j][1] - min;

				if(j != end-1)
					MPI_Send(&v[0], v.size(), MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
				else
					MPI_Send(&v[0], v.size(), MPI_INT, i, END_S, MPI_COMM_WORLD);

			}

			MPI_Send(&min, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
			MPI_Send(&shape, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
		}

	}
	else
	{

		for(int i=1; i < size; i++)
		{
			int start = listIndex[i-1];
			int end   = listIndex[i];

			numRow = end - start;
			findMinMax(start, end, listRunLen, min, max);
			numCol = max - min;

			int idx = i-1;
			MPI_Send(&idx, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
			MPI_Send(&numRow, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
			MPI_Send(&numCol, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);

			for(int j=start; j < end; j++)
			{
				vector<int> v(numCol);
				for(int k=min; k < max; k++)
				{
					v[k-min] = matrix[j][k];
				}

				if(j != end-1)
					MPI_Send(&v[0], v.size(), MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
				else
					MPI_Send(&v[0], v.size(), MPI_INT, i, END_S, MPI_COMM_WORLD);
			}

			for(int j=start; j < end; j++)
			{
				vector<int> v(numCol);
				for(int k=min; k < max; k++)
				{
					v[k-min] = matrixWeight[j][k];
				}

				if(j != end-1)
					MPI_Send(&v[0], v.size(), MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
				else
					MPI_Send(&v[0], v.size(), MPI_INT, i, END_S, MPI_COMM_WORLD);
			}

			for(int j=start; j < end; j++)
			{
				vector<int> v(2);
				v[0] = listRunLen[j][0] - min;
				v[1] = listRunLen[j][1] - min;

				if(j != end-1)
					MPI_Send(&v[0], v.size(), MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
				else
					MPI_Send(&v[0], v.size(), MPI_INT, i, END_S, MPI_COMM_WORLD);

			}

			MPI_Send(&min, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
			MPI_Send(&shape, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
		}

		for(int i=size; i < listIndex.size(); i++)
		{
			int im_free;
			MPI_Recv(&im_free, 1, MPI_INT, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD, &status);

			int start = listIndex[i-1];
			int end   = listIndex[i];

			numRow = end - start;
			findMinMax(start, end, listRunLen, min, max);
			numCol = max - min;

			int idx = i-1;
			MPI_Send(&idx, 1, MPI_INT, im_free, WORKTAG, MPI_COMM_WORLD);
			MPI_Send(&numRow, 1, MPI_INT, im_free, WORKTAG, MPI_COMM_WORLD);
			MPI_Send(&numCol, 1, MPI_INT, im_free, WORKTAG, MPI_COMM_WORLD);

			for(int j=start; j < end; j++)
			{
				vector<int> v(numCol);
				for(int k=min; k < max; k++)
				{
					v[k-min] = matrix[j][k];
				}

				if(j != end-1)
					MPI_Send(&v[0], v.size(), MPI_INT, im_free, WORKTAG, MPI_COMM_WORLD);
				else
					MPI_Send(&v[0], v.size(), MPI_INT, im_free, END_S, MPI_COMM_WORLD);
			}

			for(int j=start; j < end; j++)
			{
				vector<int> v(numCol);
				for(int k=min; k < max; k++)
				{
					v[k-min] = matrixWeight[j][k];
				}

				if(j != end-1)
					MPI_Send(&v[0], v.size(), MPI_INT, im_free, WORKTAG, MPI_COMM_WORLD);
				else
					MPI_Send(&v[0], v.size(), MPI_INT, im_free, END_S, MPI_COMM_WORLD);
			}

			for(int j=start; j < end; j++)
			{
				vector<int> v(2);
				v[0] = listRunLen[j][0] - min;
				v[1] = listRunLen[j][1] - min;

				if(j != end-1)
					MPI_Send(&v[0], v.size(), MPI_INT, im_free, WORKTAG, MPI_COMM_WORLD);
				else
					MPI_Send(&v[0], v.size(), MPI_INT, im_free, END_S, MPI_COMM_WORLD);

			}

			MPI_Send(&min, 1, MPI_INT, im_free, WORKTAG, MPI_COMM_WORLD);
			MPI_Send(&shape, 1, MPI_INT, im_free, WORKTAG, MPI_COMM_WORLD);
		}
	}

	for(int i=1; i < size; i++)
	{
		int len = -1;
		MPI_Send(&len, 1, MPI_INT, i, DIETAG, MPI_COMM_WORLD);
	}
}


// Slave processes that executes Genetic Algorithms for solving the Haplotype Assembly Problem
int slave(int rank, string pathOutput, bool verbose, bool saving, string settings)
{
	int c = 0;
	int opt, numRow, numCol, min, max;

	MPI_Status status;
	while(true)
	{
		
		MPI_Recv(&opt, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		if( status.MPI_TAG )
			break;

		MPI_Recv(&numRow, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&numCol, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		vector<vector<int> > matrix(numRow, vector<int>(numCol));
		vector<vector<int> > matrixWeight(numRow, vector<int>(numCol));
		vector<vector<int> > listRunLen(numRow, vector<int>(2));

		// receiving matrix
		int count = 0;
		while(true)
		{
			vector<int> v(numCol);
			MPI_Recv(&v[0], numCol, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			matrix[count] = v;
			count++;

			if( status.MPI_TAG == 11)
				break;

		}

		// receiving matrix weight
		count = 0;
		while(true)
		{
			vector<int> v(numCol);
			MPI_Recv(&v[0], numCol, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			matrixWeight[count] = v;
			count++;

			if( status.MPI_TAG == 11)
				break;

		}

		// receiving listRunLen list
		count = 0;
		while(true)
		{
			vector<int> v(2);
			MPI_Recv(&v[0], 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			listRunLen[count] = v;
			count++;

			if( status.MPI_TAG == 11)
				break;

		}

		MPI_Recv(&min, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&max, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		GeneticAlgorithm GA = GeneticAlgorithm(opt, rank, verbose, saving, settings);
		GA.startGA(pathOutput, matrix, matrixWeight, listRunLen, min, max, opt);

		// GA computation finished		
		c += 1;
		MPI_Send(&rank, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
	}
	return c;
}

int main(int argc, char* argv[])
{

	string pathInput  = "";
	string pathOutput = "output";
	string name       = "haplotypes";
	int gamma         = 30;
	bool verbose      = false;
	bool saving       = false;
	bool ambiguous    = true;
	string settings	  = "";
	vector<vector<int> > matrix;

	int p, value, len;
	double t1, t2;

	while ((p = getopt (argc, argv, "i:o:n:g:v:s:x:p:")) != -1)
	{
		switch (p)
		{
			case 'i':
				pathInput = optarg;
				break;
			case 'o':
				pathOutput = optarg;
				break;
			case 'p':
				settings   = optarg;
				break;
			case 'n':
				name       = optarg;
				break;
			case 'g':
				gamma      = (int) strtod(optarg, NULL);
				break;
			case 'v':
				value  = (int) strtod(optarg, NULL);
				if(value==1)
					verbose = true;
				else
					verbose = false;
				break;
			case 's':
				value   = (int) strtod(optarg, NULL);
				if(value==1)
					saving = true;
				else
					saving = false;
				break;
			case 'x':
				value   = (int) strtod(optarg, NULL);
				if(value==1)
					ambiguous = true;
				else
					ambiguous = false;
				break;
			}
	}

	// Initialize the MPI environment
	MPI_Init(NULL, NULL);

	// Get the number of processes
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();

	if(rank == 0)
	{
		if(argc == 1)
		{
			cout << "************************************************************************************************************************\n";
			cout << "GenHap: Genetic Algorithm for Haplotype Assembly, parameters" << endl;
			cout << " * -i: path to the input wif file (mandatory)" << endl;
			cout << " * -o: output folder (optional, default value: " << pathOutput << ")" <<endl;
			cout << " * -n: file name containing the estimated haplotypes (optional, default value " << name << ")" << endl;
			cout << " * -g: number of reads composing each sub-matrix (optional, default value " << gamma << ")" << endl;
			cout << " * -v: verbose modality (optional, default value " << 0 << ")" << endl;
			cout << " * -s: save partitions and fragment matrix (optional, default value " << 1 << ")" << endl;
			cout << " * -x: mask ambiguous positions in the output haplotypes with a X (optional, default value " << 1 << ")" << endl;
			cout << " * -p: path to the file of GA settings" << endl;
			cout << "************************************************************************************************************************\n";
			
			for(int i=1; i < size; i++)
			{
				int len = -1;
				MPI_Send(&len, 1, MPI_INT, i, DIETAG, MPI_COMM_WORLD);
			}
			exit(-1);
		}

		if(strcmp(pathInput.c_str(), "") == 0)
		{
			cout << " * The path to the input wif file is mandatory" << endl;
			cout << " * Please, run GenHap without parameters to obtain information" << endl;

			for(int i=1; i < size; i++)
			{
				int len = -1;
				MPI_Send(&len, 1, MPI_INT, i, DIETAG, MPI_COMM_WORLD);
			}
			exit(-2);
		}

		// printing parameters
		cout << "************************************************************************************************************************\n";
		cout << "GenHap: Genetic Algorithm for Haplotype Assembly\n" << endl;
		cout << "* Using " << pathInput  << " as input wif file" << endl;
		cout << "* Using " << pathOutput << " as output folder" << endl;
		cout << "* Using " << name       << " as file name containing the estimated haplotypes" << endl;
		cout << "* Using " << gamma      << " reads for each sub-matrix" << endl;

		if(verbose)
			cout << "* Verbose modality enabled" << endl;
		else
			cout << "* Verbose modality disabled" << endl;
		if(saving)
			cout << "* Saving partitions and fragment matrix" << endl;
		else
			cout << "* Not saving partitions and fragment matrix" << endl;
		if(ambiguous)
			cout << "* Masking ambiguous positions" << endl;
		else
			cout << "* Not masking ambiguous positions" << endl;
		
		//creation output folder
		string folder1;
		struct stat sb;
		if( !(stat(pathOutput.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)))
		{
			if (os == 0)
			{
				folder1 = "MD " + pathOutput;
				int sis = system(folder1.c_str());
			}
			else
			{
				folder1 = "mkdir " + pathOutput;
				int sis = system(folder1.c_str());
			}
		}

		folder1 = pathOutput + slash + "output";

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

		vector<vector<int> > matrixWeight;
		vector<int> listIndex;
		vector<vector<int> > listRunLen;

		Reader reader;
		reader.readFromWif(pathInput, pathOutput, gamma, saving, matrix, matrixWeight, listIndex, listRunLen);

		len = listIndex.size()-1;

		cout << "* Number of reads: " << matrix.size() << endl;
		cout << "* Number of columns:  " << matrix[0].size() << endl;
		cout << "* Start " << listIndex.size()-1 << " optimizations" << endl;
		cout << "************************************************************************************************************************\n";
		
		master(matrix, matrixWeight, listRunLen, listIndex, size);

	}
	else
	{
		int c = slave(rank, pathOutput, verbose, saving, settings);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0)
	{	

		saveHaplotypes(pathOutput, name, len, matrix, ambiguous);
		t2 = MPI_Wtime();
		double elapsed = t2-t1;

		printf("Elapsed time %5.2f seconds\n", elapsed);
		cout << "************************************************************************************************************************\n";
	}

	// Finalize the MPI environment.
    MPI_Finalize();

}