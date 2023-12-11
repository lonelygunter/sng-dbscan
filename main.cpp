#include <iostream>
#include <fstream>	// file io
#include <vector>
#include <sstream>	// stringstream
#include <cstdlib>	// rand()
#include <string>
using namespace std;

// declare functions:
vector<vector<string> > readDataset(string pathDataset, int s);
int numLines(ifstream& dataset);
vector<int> sampleInstances(int n, int numLines);
vector<string> splitString(string line);

int main() {
	// sampling dataset
	vector<vector<string> > instances;
	instances = readDataset("datasets/iris/iris.data", 40);

	// inizialization matrix for graph
	vector<vector<int> > matrix(instances.size(), vector<int>(instances.size(), 0));
	
	return 0;
}

// functions:

/* function to read a dataset:
	pathDataset:	path del dataset
	s:				dimension of the sample
	instances:		vector of vector to store the istances and their datas
	dataset:		dataset file
	numInstances:	# of lines (instances) of the dataset
	vecInstances:	vector of # of the instances sampled
	currentLine:	current line pointer from the pointer of the file
	ss:				string stream of the current line
	element:		single data of the readed line
	instance:		lincked vector of single instance
*/
vector<vector<string> > readDataset(string pathDataset, int s){
	vector<vector<string> > instances;
	ifstream dataset(pathDataset);
	int numInstances;
	vector<int> vecInstances;
	string currentLine;

	// file check
	if (!dataset.is_open()){
		cout << "Unable to open file";
	}
	
	// take the number of the total instances
	numInstances = numLines(dataset);

	// random sample numbers wothout rep
	vecInstances = sampleInstances(s, numInstances);
	sort(vecInstances.begin(), vecInstances.end());

	// reset file pointer
	dataset.seekg(0, ios::beg);

	// random sample the dataset:
	// jump to first sampled line
	for (int i = 0; i < vecInstances[0]; i++) {
		if (!getline(dataset, currentLine)) {
			cerr << "Error reading the file or the file is too short." << endl;
			dataset.close();
		}
	}

	// add the instance element by element
	instances.push_back(splitString(currentLine));

	// jump to other sampled line
	for (int i = 1; i < vecInstances.size(); i++){
		for (int j = 0; j < vecInstances[i]-vecInstances[i-1]; j++) {
			if (!getline(dataset, currentLine)) {
				cerr << "Error reading the file or the file is too short." << endl;
				dataset.close();
			}
		}

		// add the instance element by element
		instances.push_back(splitString(currentLine));
	}

	dataset.close();
	return instances;
}

/* function to read the # of lines (instances) into a file
	numLines:	# of lines into the file
*/
int numLines(ifstream& dataset){
	int numLines = 1; // init to 1 cos last line in file not have '\n'

	numLines += count(istreambuf_iterator<char>(dataset), istreambuf_iterator<char>(), '\n');

	return numLines;
}

/* function to take a n random lines (instances) from a file without rep
	n:				# of lines (instances) to take
	numLines:		# of lines into the file
	lines:			vector with all # of lines
	linesSize:		original # of the lines
*/
vector<int> sampleInstances(int n, int numLines){
	// create a vector to semplify generation of random number without rep
	vector<int> lines;
	for (int i = 0; i < numLines; i++){
		lines.push_back(i);
	}
	
	// take random elements from "lines"
	vector<int> numLineInstances;
	int linesSize = lines.size();
	while (lines.size() != linesSize - n){
		int randline = rand() % lines.size();
		numLineInstances.push_back(lines[randline]);

		lines.erase(lines.begin() + randline);
	}
	
	return numLineInstances;
}

/* function to split a string
	line:	string to split
*/
vector<string> splitString(string line){
	stringstream ss(line);
	string element;
	vector<string> instance;

	// put all datas of an instance in a  vector
	while (getline(ss, element, ',')){
		instance.push_back(element);
	}

	return instance;
}