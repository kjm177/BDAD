#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <fstream>
#include <mpi.h>
using namespace std;

vector<vector<float>> coefficients;
vector<int> unknowns;
//vector<int> constants;

float maxError;
int numOfUnknowns;


int convertStringToInteger(string s)
{
    stringstream convert(s);
    int ans;
    convert >> ans;
    return(ans);
}

float convertStringToFloat(string s)
{
    stringstream convert(s);
    float ans;
    convert >> ans;
    return(ans);
}

vector<string> tokenizer(string input)
{
    vector <string> tokens;
    stringstream check1(input);
    string intermediate;
    while(getline(check1, intermediate, ' '))
    {
        tokens.push_back(intermediate);
    }
    return tokens;
}

void readInput(string inputFileName)
{
    ifstream inputFile;
    inputFile.open(inputFileName);
    if(!inputFile.is_open())
        cout<<"ERROR! Unable to open file!";
    else
    {
        string token;
        getline(inputFile, token);
        numOfUnknowns = convertStringToInteger(token);
        getline(inputFile, token);
        maxError = convertStringToFloat(token);
        getline(inputFile, token);
        vector<string> temp = tokenizer(token);
        for(int i = 0 ; i < numOfUnknowns ; i++)
            unknowns.push_back(convertStringToFloat(temp[i]));

        for(int i = 0 ; i < numOfUnknowns ; i++)
        {
            vector<float> currentLine;
            getline(inputFile, token);
            vector<string> tokens = tokenizer(token);
            for(int j = 0 ; j < tokens.size() ; j++)
                currentLine.push_back(convertStringToFloat(tokens[j]));
            coefficients.push_back(currentLine);
        }
    }
    return;
}

bool validError(const vector<float>& current)
{
    for(int i = 0 ; i < numOfUnknowns ; i++)
    {
        if( fabs(current[i] - unknowns[i])/current[i] > maxError)
            return false;
    }
    return true;
}

void solveEquation(vector<float>& y, int myRank)
{
    float sum = 0;
    for(int i = 0 ; i < numOfUnknowns ; i++)
    {
        if(i != myRank)
        {
            sum += coefficients[myRank][i] * unknowns[i];
        }
    }
    y[myRank] = (coefficients[myRank][numOfUnknowns] - sum)/coefficients[myRank][myRank];
}

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        cout<<"ERROR! Input file missing!"<<endl;
        return 0;
    }

    int comm_sz, my_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Status stat;

	string filename = "input.txt";

	int numOfProcesses = comm_sz;
	int iterations = 0;

	if(numOfProcesses == 1)
    {
        vector<float> y;
        for(int i = 0 ; i < numOfUnknowns ; i++)
        {
            y.push_back(0);
        }

        do
        {
            for(int i = 0 ; i < numOfUnknowns ; i++)
                solveEquation(y, i);

            iterations++;

            for(int i = 0 ; i < numOfUnknowns ; i++)
                unknowns[i] = y[i];
        }
        while((!validError(y)));

    }
    else
    {
        cout<<"not implemented"<<endl;
    }

	MPI_Finalize();

    return 1;
}
