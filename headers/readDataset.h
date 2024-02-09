#ifndef READ_DATASET_H
#define READ_DATASET_H

#include <fstream>	// file io
#include <vector>
#include <sstream>	// stringstream

using namespace std;

// function to read a dataset
vector<vector<double> > readDataset(string pathDataset, int n);

// function to read the # of lines (instances) into a file
int numLines(ifstream& dataset);

// function to take a n random lines (instances) from a file without rep
vector<vector<double> > sampleInstances(int n, int numLines, ifstream& dataset);

// function to split a string
vector<double> splitString(int index, string line);

// function to check if a string is a digit:
bool isStringDigit(string str);

#endif // READ_DATASET_H
