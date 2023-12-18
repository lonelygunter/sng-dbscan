#include <iostream>
#include <fstream>	// file io
#include <vector>
#include <sstream>	// stringstream
#include <cstdlib>	// rand()
#include <string>
#include <cmath>	// math
#include <getopt.h>	// getopt
#include <numeric>

#include <chrono>

using namespace std;

// declare functions:
void getOpt(int& n, double& s, double *epsilon, int& l, int& minpts, int argc, char **argv);
vector<vector<double> > readDataset(string pathDataset, int n);
int numLines(ifstream& dataset);
vector<int> sampleInstances(int n, int numLines);
vector<double> splitString(string line);
bool isStringDigit(string str);
void compareWsnSample(vector<vector<double> > instances, double s, double n, double epsilon[], vector<vector<int> >& graph);
double euclDist(vector<double> point1, vector<double> point2);
void connectedComponents(vector<vector<int> > graph, int minpts, int l, int totInstances, vector<vector<vector<int> > >& k);
int calcLinks(int instance, vector<vector<int> >& graph);
void findConnComp(int i, vector<vector<int> >& graph, vector<vector<int> > graphCpy, vector<int>& minPtsInstances, int totInstances, vector<vector<int> >& k, vector<int>& usedMinPtsInstances);
bool isInMinPts(int instance, vector<int> minPtsInstances);
void setAllToZero(vector<vector<int> >& graph, vector<int> minPtsInstances, int instance);
int findMinPtsInstances(int instance, vector<int> minPtsInstances);

/*
	n:			dimension of the 1st sample
	s:			dimension of the 2nd sample
	epsilon:	range for a point
	l:			# of classes
	minpts:		min # of points in the epsilon range of a point
	k:			subgraph induced by vertices of degree at least MinPts
*/
int main(int argc, char **argv){
	chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();

	vector<vector<double> > instances;
	int n = 30;
	double s = 0.5;
	double epsilon[2] = {0.2, 2.4};
	int l = 8;
	int minpts = 1;

	// getopt:
	getOpt(n, s, epsilon, l, minpts, argc, argv);

	// 1. sampling dataset
	instances = readDataset("datasets/iris/iris.data", n);

	// 2. inizialization of a matrix for graph
	vector<vector<int> > graph(instances.size(), vector<int>(instances.size(), 0));

	// 3. check if a point is in epsilon range
	compareWsnSample(instances, s, n, epsilon, graph);
	
	// 4. inizialization of an array of l items
	vector<vector<vector<int> > > k;

	// print graph
	// for (int i = 0; i < instances.size(); ++i) {
    //     for (int j = 0; j < instances.size(); ++j) {
    //         cout << graph[i][j] << " ";
    //     }
    //     cout << endl;
    // }

	// 5. create the connected components
	k.resize(l);
	connectedComponents(graph, minpts, l, instances.size(), k);

	// print k
	// for (size_t i = 0; i < k.size(); ++i) {
    //     for (size_t j = 0; j < k[i].size(); ++j) {
    //         for (size_t m = 0; m < k[i][j].size(); ++m) {
    //             cout << k[i][j][m] << " ";
    //         }
    //         cout << "\n";
    //     }
    //     cout << endl;
    // }

	// 6. inizialization of an array of l items
	vector<vector<vector<int> > > c;

	// 7. create the clusters
	c.resize(l);
	

	chrono::high_resolution_clock::time_point t_end = chrono::high_resolution_clock::now();
	cout << "Total time required = " << chrono::duration<double, milli>(t_end-t_start).count() << endl;

	return 0;
}

// functions:

/* function to get options for command line:
*/
void getOpt(int& n, double& s, double *epsilon, int& l, int& minpts, int argc, char **argv){
	int opt;

    while ((opt = getopt(argc, argv, "n:s:e:l:m:")) != -1) {
        switch (opt) {
        case 'n':
            n = stod(optarg);
			cout << "n: " << n << endl;
            break;
        case 's':
            s = stod(optarg);
			cout << "s: " << s << endl;
            break;
        case 'e': {
			vector<double> eps = splitString(optarg);
			epsilon[0] = eps[0];
			epsilon[1] = eps [1];

			cout << "epsilon: [" << epsilon[0] << ", " << epsilon[1] << "]" << endl;
			}
            break;
        case 'l':
            l = stod(optarg);
			cout << "l: " << l << endl;
            break;
        case 'm':
            minpts = stod(optarg);
			cout << "minpts: " << minpts << endl;
            break;
        case '?':
			cout << "This is not an implemented option" << endl;
            break;
        default:
            abort ();
        }
    }
}

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
		int instCpysnSize = instancesCpy.size() - s*n;

		while (instancesCpy.size() != instCpysnSize){
			int randline = rand() % instancesCpy.size();

			// to not randomly take the same instance
			while (i == randline){
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

/* function to create the connected components:
	minPtsInstances:		vector of istances that have minPts vectors connected
	usedMinPtsInstances:	vector to track which instance was used
	graphCPy:				graph where can delete 1s without influence k items
*/
void connectedComponents(vector<vector<int> > graph, int minpts, int l, int totInstances, vector<vector<vector<int> > >& k){
	vector<int> minPtsInstances;
	vector<int> usedMinPtsInstances;
	vector<vector<int> > graphCpy = graph;

	// check all instances with MinPts vartices
	for (int i = 0; i < totInstances; i++){
		int calcololinks = calcLinks(i, graph);
		if (calcololinks >= minpts){
			minPtsInstances.push_back(i);
		}
	}

	// fill k with all k_i subgraphs
	for (int i = 0; i < l; i++){
		if (minPtsInstances.size() != 0){
			// add manually the 1st instance
			k[i].push_back(graphCpy[minPtsInstances[0]]);

			findConnComp(0, graph, graphCpy, minPtsInstances, totInstances, k[i], usedMinPtsInstances);

			// track all used minPtsInstances in this k_i
			for (int j = 0; j < usedMinPtsInstances.size(); j++){
				int eraseMinPts = findMinPtsInstances(usedMinPtsInstances[j], minPtsInstances);
				if (eraseMinPts >= 0){
					minPtsInstances.erase(minPtsInstances.begin() + eraseMinPts);
				}
			}
		}
	}
}

/* function to calculate # of links into an istance in the graph:
	sum:		summarize of the 1s of a graph row
*/
int calcLinks(int instance, vector<vector<int> >& graph){
	int sum = 0;

	for(int i = 0; i < graph[instance].size(); i++){
		sum += graph[instance][i];
	}

	return sum;
}

/* function to find the connected components:
	i:						index of the current instance
	graph:					graph where can delete 1s without influence k items
	graphCPy:				original graph
	minPtsInstances:		vector of all minPts instances
	totInstances:			size of minPtsInstances
	k:						vector where put all instaces
	usedMinPtsInstances:	vector with used minPts instances
	zero:					bool to check if current instance was used or not
*/
void findConnComp(int i, vector<vector<int> >& graph, vector<vector<int> > graphCpy, vector<int>& minPtsInstances, int totInstances, vector<vector<int> >& k, vector<int>& usedMinPtsInstances){
	bool zero = false;

	for (int j = 0; j < totInstances; j++){
		if (graph[minPtsInstances[i]][j] == 1){
			k.push_back(graphCpy[j]);

			// step into only in instaces that have minPts vertices
			if (isInMinPts(j, minPtsInstances)){
				setAllToZero(graph, minPtsInstances, j);

				// to shorten the recursive sequence
				setAllToZero(graph, minPtsInstances, minPtsInstances[i]);
			
				// notice that all instances connection with the corrent instance was set to zero
				zero = true;

				// track the used minPtsIntance
				if (!isInMinPts(minPtsInstances[i], usedMinPtsInstances)){
					usedMinPtsInstances.push_back(minPtsInstances[i]);
				}
				
				findConnComp(findMinPtsInstances(j, minPtsInstances), graph, graphCpy, minPtsInstances, totInstances, k, usedMinPtsInstances);
			}

			setAllToZero(graph, minPtsInstances, j);
		}
	}

	// set to zero all instances connection with the corrent instance
	if (!zero){
		// track the used minPtsIntance
		if (!isInMinPts(minPtsInstances[i], usedMinPtsInstances)){
			usedMinPtsInstances.push_back(minPtsInstances[i]);
		}
		
		setAllToZero(graph, minPtsInstances, minPtsInstances[i]);
	}
}

/* function to check if an instance is in minPts range:
	instance:		number of the instance to check
*/
bool isInMinPts(int instance, vector<int> minPtsInstances){
	for (int i = 0; i < minPtsInstances.size(); i++){
		if (minPtsInstances[i] == instance){
			return true;
		}
	}
	
	return false;
}

/* function to shorten the recursive function fincConnComp:
	instance:		number of the column instance of the graph
*/
void setAllToZero(vector<vector<int> >& graph, vector<int> minPtsInstances, int instance){
	for (int i = 0; i < minPtsInstances.size(); i++){
		graph[minPtsInstances[i]][instance] = 0;
	}
}

/* function to find the id of an instance into minPtsInstances:
	instance:		number of ther instance to find
*/
int findMinPtsInstances(int instance, vector<int> minPtsInstances){
	for (int i = 0; i < minPtsInstances.size(); i++){
		if (minPtsInstances[i] == instance){
			return i;
		}
	}
	
	return -1;
}
