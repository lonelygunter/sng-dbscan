#include <iostream>
#include <fstream>	// file io
#include <vector>
#include <sstream>	// stringstream
#include <cstdlib>	// rand()
#include <string>
#include <cmath>

#include <chrono>

using namespace std;

// declare functions:
vector<vector<double> > readDataset(string pathDataset, int n);
int numLines(ifstream& dataset);
vector<int> sampleInstances(int n, int numLines);
vector<double> splitString(string line);
bool isStringDigit(string str);
void compareWsnSample(vector<vector<double> > instances, double s, double n, double epsilon[], vector<vector<int> >& graph);
double euclDist(vector<double> point1, vector<double> point2);

int main(){
	auto t_start = chrono::high_resolution_clock::now();

	// sampling dataset
	vector<vector<double> > instances;
	int n = 40;
	double s = 0.4; // sample S
	double epsilon[2] = {0.2, 2.5};

	instances = readDataset("datasets/iris/iris.data", n);

	// inizialization matrix for graph
	vector<vector<int> > graph(instances.size(), vector<int>(instances.size(), 0));

	// check if a point is in epsilon range
	compareWsnSample(instances, s, n, epsilon, graph);

	auto t_end = chrono::high_resolution_clock::now();
	cout << "Total time required = " << chrono::duration<double, milli>(t_end-t_start).count() << endl;

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
vector<vector<double> > readDataset(string pathDataset, int n){
	vector<vector<double> > instances;
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
	vecInstances = sampleInstances(n, numInstances);
	sort(vecInstances.begin(), vecInstances.end());

	// reset file pointer
	dataset.seekg(0, ios::beg);

	// random sample the dataset:
	// jump to first sampled line
	for (int i = 0; i <= vecInstances[0]; i++) {
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
vector<double> splitString(string line){
	stringstream ss(line);
	string element;
	vector<double> instance;

	// put all datas of an instance in a  vector
	while (getline(ss, element, ',')){
		if (isStringDigit(element)){
			instance.push_back(stod(element));
		}
	}

	return instance;
}

/* function to check if a string is a digit:
	str:	string to check
	ss:		string stream
	digit:	digit variable
*/
bool isStringDigit(string str){
    stringstream ss(str);
    double digit;
    ss >> digit;

    return ss.eof() && !ss.fail();
}

/* function that take one point of the "n" sample to compare it with sampled sn points:
*/
void compareWsnSample(vector<vector<double> > instances, double s, double n, double epsilon[], vector<vector<int> >& graph){
	for (int i = 0; i < instances.size(); i++){
		// sample of s
		vector<vector<double> > instancesCpy = instances;
		int instCpySize = instancesCpy.size();

		while (instancesCpy.size() != instCpySize - s*n){
			int randline = rand() % instancesCpy.size();

			// to not randomly take the same instance
			while (instances[i] == instancesCpy[randline]){
				randline = rand() % instancesCpy.size();
			}

			// calculate Euclidean distance for all instance infos
			double euclDistij = euclDist(instances[i], instancesCpy[randline]);

			// check if point are in range epsilon
			if (euclDistij >= epsilon[0] && euclDistij <= epsilon[1]){
				graph[i][randline] = 1;
			}

			instancesCpy.erase(instancesCpy.begin() + randline);
		}
	}
}

/* function that calculate the Euclidean distance between two points
	point1:			1st point of the distance
	point2: 		2nd point of the distance
	euclDistij:		Euclidean distance
*/
double euclDist(vector<double> point1, vector<double> point2){
	double euclDistij = 0;
	for (int j = 0; j < point1.size(); j++){
		euclDistij += pow(point1[j] - point2[j], 2);
	}
	
	euclDistij = sqrt(euclDistij);

	return euclDistij;
} 