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

int Reader::get_first(vector<int> &row)
{
	int i=0;
	for(int j=0; j < row.size(); j++)
	{
		if(row[j]!=-1)
		{
			i=j;
			break;
		}
	}
	return i;
}

int Reader::get_last(vector<int> &row)
{
	int i=0;
	for(int j=0; j < row.size(); j++)
	{
		if(row[j]!=-1)
		{
			i=j;
		}
	}
	return i;
}

void Reader::readFromWif(string pathIn, string pathOut, int gamma, bool saveMatrix,
				vector<vector<int> > &matrixIn, vector<vector<int> > &matrixWeightIn,
				vector<vector<int> > &listIndexIn, vector<vector<int> > &listRunLenIn,
				bool verbose)
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
		string s = pathOut + slash + "matrix";
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
					if(j == matrix[i].size()-1)
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

	cout << " * Number of reads: " << matrix.size() << endl;
	cout << " * Number of SNPs:  " << matrix[0].size() << endl;

// ********************************* Finding haplotype blocks *************************************************

	int minCov    = numeric_limits<int>::max();
	double avgCov = 0.0;
	int maxCov    = 0;
	int cov       = 0;

	for(int j=0; j < matrix[0].size(); j++)
	{
		cov = 0;
		for(int i=0; i < matrix.size(); i++)
		{
			if(matrix[i][j] != -1)
			{
				cov += 1;
				avgCov += 1;
			}
		}

		if(maxCov < cov)
			maxCov = cov;
		if(minCov > cov)
			minCov = cov;
	}

	avgCov = round(avgCov/(double) matrix[0].size());

	vector<int> listBlocks;
	listBlocks.push_back(0);

	int begin = get_first(matrix[0]);
	int end   = get_last(matrix[0]);
	int pos   = 1;

	while(pos < matrix.size())
	{
		int new_begin = get_first(matrix[pos]);
		int new_end   = get_last(matrix[pos]);

		if(new_begin>end)
		{
			begin = new_begin;
			end   = new_end;
			listBlocks.push_back(pos);
			if(verbose)
				printf(" \t* New haplotype block found at column %d, %d is the last read before the block\n", new_begin, pos-1);

		}
		else
		{
			if(new_begin>=begin)
			{
				if(new_end>end)
					end = new_end;
			}
		}
		pos++;
	}
	listBlocks.push_back(pos);

	if(listBlocks.size()-1 == 1)
		printf("\n * Detected 1 haplotype block\n");
	else
		printf("\n * Detected %d haplotype blocks\n", (int) listBlocks.size()-1);

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

	vector<int> steps;
	if(gamma != 0)
	{
		for(int k=0; k < listBlocks.size(); k++)
			steps.push_back(gamma);
	}
	else
	{
		for(int k=0; k < listBlocks.size()-1; k++)
		{
			begin = listBlocks[k];
			end   = listBlocks[k+1];

			int firstCol = get_first(matrix[begin]);
			int lastCol  = get_last(matrix[end-1]);

			double avgCov = 0;
			for(int j=firstCol; j <= lastCol; j++)
			{
				for(int i=begin; i < end; i++)
				{
					if(matrix[i][j] != -1)
						avgCov++;
				}
			}
			if ((lastCol-firstCol) != 0)
				avgCov = round(avgCov/(double) (lastCol-firstCol));
			else
				avgCov = round(avgCov);

			steps.push_back((int)avgCov);
		}
	}

	vector<vector<int> > listIndex;

	for(int j=0; j < listBlocks.size()-1; j++)
	{
		begin = listBlocks[j];
		end   = listBlocks[j+1];

		int step = steps[j];

		vector<int> pos;
		for(int i=begin; i <= end; i+=step)
			pos.push_back(i);

		if(pos.size() == 0)
			pos.push_back(end);
		else
		{
			if(pos[pos.size()-1] != end)
				pos.push_back(end);
		}
			
		listIndex.push_back(pos);
	}

	FILE * fl;
	string s = pathOut + slash + "informations";
	fl=fopen(s.c_str(), "w");

	if( fl==NULL )
	{
		printf("Error open output file: %s\n", s.c_str());
		exit(-9);
	}

	if(verbose)
	{
		cout << " * Maximum coverage detected " << maxCov << endl;
		cout << " * Minimum coverage detected " << minCov << endl;
		cout << " * Average coverage detected " << (int) avgCov << endl;
		for(int j=0; j < listBlocks.size()-1; j++) 
			printf(" * Block (%d) gamma equal to %d for each sub-matrix\n", j, (int) steps[j]);
	}

	fprintf(fl, "Number of reads: %d\n", (int) matrix.size());
	fprintf(fl, "Number of columns: %d\n", (int) matrix[0].size());
	fprintf(fl, "Maximum coverage: %d\n", maxCov);
	fprintf(fl, "Minimum coverage: %d\n", minCov);
	fprintf(fl, "Average coverage: %d\n", (int) avgCov);
	fprintf(fl, "Number of haplotype blocks: %d\n", (int) listIndex.size());
	for(int j=0; j < listBlocks.size()-1; j++) 
		fprintf(fl, "Block (%d) gamma equal to %d for each sub-matrix\n", j, (int) steps[j]);
	fclose(fl);

	matrixIn = matrix;
	matrixWeightIn = matrixWeight;
	listIndexIn = listIndex;
	listRunLenIn = listRunLen;

}