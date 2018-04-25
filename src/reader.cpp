#ifdef _WIN32
	#include "headers\\classes.hpp"
	const string slash = "\\";
#elif _WIN64
	#include "headers\\classes.hpp"
	const string slash = "\\";
#else
	#include "headers/classes.hpp"
	const string slash = "/";
#endif

void Reader::readFromWif(string pathIn, string pathOut, int gamma, bool saveMatrix,
				vector<vector<int> > &matrixIn, vector<vector<int> > &matrixWeightIn,
				vector<int> &listIndexIn, vector<vector<int> > &listRunLenIn)
{
// ********************************* Generating matrices from WIF file *********************************

	vector<vector<string> > lista;
	string line;
	ifstream infile(pathIn.c_str());
	while (getline(infile, line))
	{
		istringstream buf(line);
		vector<string> sublista;
		for(string token; getline(buf, token, ' '); )
		{
			sublista.push_back(token);
		}
		lista.push_back(sublista);
	}

	vector<vector<vector<int> > > SNPValues;
	vector<int> SNPColums;

	for(int i=0; i<lista.size(); i++)
	{
		vector<vector<int> > values;
		for(int j=0; j<lista[i].size(); j++)
		{
			vector<int> l;
			if(j % 5 == 0)
			{
				if(lista[i][j].compare("#") != 0)
				{
					l.push_back((int) strtod(lista[i][j].c_str(), NULL));
					SNPColums.push_back((int) strtod(lista[i][j].c_str(), NULL));
					l.push_back((int) strtod(lista[i][j+2].c_str(), NULL));
					l.push_back((int) strtod(lista[i][j+3].c_str(), NULL));
					values.push_back(l);
				}
			}
		}
		SNPValues.push_back(values);
	}

	sort(SNPColums.begin(), SNPColums.end());

	SNPColums.erase(unique(SNPColums.begin(), SNPColums.end()), SNPColums.end());

	sort(SNPValues.begin(), SNPValues.end(),
		[](const vector<vector<int> > &lhs, const vector<vector<int> > &rhs)
		{
			return lhs[0][0] < rhs[0][0];
		}
	);

	vector<vector<int> > matrix(SNPValues.size(), vector<int>(SNPColums.size()));
	vector<vector<int> > matrixWeight(SNPValues.size(), vector<int>(SNPColums.size()));

	for(int i=0; i < SNPValues.size(); i++)
	{
		for(int j=0; j < SNPColums.size(); j++)
		{
			matrix[i][j] = -1;
			matrixWeight[i][j] = 0;
		}

	}

	int k = 0;
	for(int i=0; i < SNPValues.size(); i++)
	{
		for(int j=0; j < SNPValues[i].size(); j++)
		{
			for(int l=0; l < SNPColums.size(); l++)
			{
				if(SNPColums[l] == SNPValues[i][j][0])
				{
					k = l;
					break;
				}
			}
			matrix[i][k] = SNPValues[i][j][1];
			matrixWeight[i][k] = SNPValues[i][j][2];
		}
	}

	if(saveMatrix)
	{
		FILE * fl;
		string s = pathOut + slash + "output" + slash + "matrix";
		fl=fopen(s.c_str(), "w");

		if( fl==NULL )
		{
			printf("Error open output file: %s\n", s.c_str());
			exit(-9);
		}

		for(int i=0; i < matrix.size(); i++)
		{
			for(int j=0; j < matrix[i].size(); j++)
			{
				if(matrix[i][j]==-1)
				{
					if(j == matrixWeight[i].size()-1)
						fprintf(fl, "%s", "-");
					else
						fprintf(fl, "%s ", "-");
				}
				else
				{
					if(j == matrix[i].size()-1)
						fprintf(fl, "%d", matrix[i][j]);
					else
						fprintf(fl, "%d ", matrix[i][j]);
				}
			}
			fprintf(fl, "\n");
		}
		fclose(fl);
	}

	int maxCov = 0;
	int cov    = 0;
	for(int j=0; j < matrix[0].size(); j++)
	{
		cov = 0;
		for(int i=0; i < matrix.size(); i++)
		{
			if(matrix[i][j] != -1)
				cov += 1;
		}
		if(maxCov < cov)
			maxCov = cov;
	}

	FILE * fl;
	string s = pathOut + slash + "informations";
	fl=fopen(s.c_str(), "w");

	if( fl==NULL )
	{
		printf("Error open output file: %s\n", s.c_str());
		exit(-9);
	}

	fprintf(fl, "Number of reads: %d\n", (int) matrix.size());
	fprintf(fl, "Number of columns: %d\n", (int) matrix[0].size());
	fprintf(fl, "Maximum coverage: %d", maxCov);
	fclose(fl);

// ********************************* Generating structures for parallelization *********************************
	vector<vector<int> > listRunLen;

	for(int i=0; i < matrix.size(); i++)
	{
		vector<int> pos;
		vector<int> tmp;
		for(int j=0; j < matrix[i].size(); j++)
		{
			if(matrix[i][j] != -1)
				pos.push_back(j);
		}

		tmp.push_back(pos[0]);
		tmp.push_back(pos[pos.size()-1]+1);
		listRunLen.push_back(tmp);
	}

	vector<int> listIndex;

	for(int i=0; i < matrix.size(); i+=gamma)
		listIndex.push_back(i);

	if(listIndex[listIndex.size()-1] != (matrix.size()))
		listIndex.push_back(matrix.size());

	matrixIn = matrix;
	matrixWeightIn = matrixWeight;
	listIndexIn = listIndex;
	listRunLenIn = listRunLen;
}