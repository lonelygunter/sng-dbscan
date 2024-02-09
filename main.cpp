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
void getOpt(string& filePath, int& n, double& s, double *epsilon, int& l, int& minpts, int argc, char **argv);
vector<vector<double> > readDataset(string pathDataset, int n);
int numLines(ifstream& dataset);
vector<vector<double> > sampleInstances(int n, int numLines, ifstream& dataset);
vector<double> splitString(int index, string line);
bool isStringDigit(string str);
void compareWsnSample(vector<vector<double> > instances, double s, int n, double epsilon[], vector<vector<int> >& graph, int minPts, vector<int>& minPtsInstances);
double euclDist(vector<double> point1, vector<double> point2);
void connectedComponents(vector<vector<int> > graph, int l, vector<vector<vector<double> > >& k, vector<vector<double> > instances, vector<int>& minPtsInstances, vector<vector<vector<double> > >& c);
void findConnComp(int i, vector<vector<int> >& graph, vector<vector<double> > instances, vector<int>& minPtsInstances, vector<int> originalMinPtsInstances, vector<vector<double> >& k, vector<int>& usedMinPtsInstances, vector<vector<vector<double> > >& c, vector<vector<int> >& originalGraph, int l);
bool isInMinPts(int instance, vector<int> minPtsInstances);
void eraseAll(vector<vector<int> >& graph, vector<int> minPtsInstances, int instance);
int findMinPtsInstances(int instance, vector<int> minPtsInstances);
void findClusters(vector<vector<double> > instances, int l, vector<vector<vector<double> > >& c, vector<vector<int> >& graph, int ki);
void printLoadingBar(int index, int totalIterations);

/*
	instances:	istances randomly sampled
	n:			dimension of the 1st sample
	s:			dimension of the 2nd sample
	epsilon:	range for a point
	l:			# of classes
	minpts:		min # of points in the epsilon range of a point
	filePath:	path of the dataset
	graph:		graph of the dataset
	k:			subgraph induced by vertices of degree at least MinPts
	c:			vector of the clusters
*/
int main(int argc, char **argv){
	vector<vector<double> > instances;
	vector<int> minPtsInstances;

	ofstream outFile("new.txt");

	// iris150
	int n = 150;
	double s = 0.5;
	double epsilon[2] = {0.2, 2.9};
	int l = 10;
	int minPts = 2;
	string filePath = "datasets/iris150.data";

	// eshop100k
	// int n = 3000;
	// double s = 0.4;
	// double epsilon[2] = {15, 21};
	// int l = 10;
	// int minPts = 2;
	// string filePath = "datasets/eshop100k.data";

	// getopt:
	getOpt(filePath, n, s, epsilon, l, minPts, argc, argv);

	chrono::high_resolution_clock::time_point t_start_total = chrono::high_resolution_clock::now();

	// sampling n instances from dataset
	chrono::high_resolution_clock::time_point t_start_instances = chrono::high_resolution_clock::now();

	cout << "\nSampling " << n << " (n) instances..." << endl;
	instances = readDataset(filePath, n);

	chrono::high_resolution_clock::time_point t_end_instances = chrono::high_resolution_clock::now();
	cout << "\n" << chrono::duration<double, milli>(t_end_instances-t_start_instances).count()/1000 << " sec\n" << endl;

	// inizialization of a matrix for graph
	vector<vector<int> > graph;
	graph.resize(instances.size());

	// sampling sn instances from dataset and check minPts
	chrono::high_resolution_clock::time_point t_start_sn = chrono::high_resolution_clock::now();

	cout << "\nSampling " << s*n << " (sn) instances..." << endl;
	compareWsnSample(instances, s, n, epsilon, graph, minPts, minPtsInstances);

	chrono::high_resolution_clock::time_point t_end_sn = chrono::high_resolution_clock::now();
	cout << "\n" << chrono::duration<double, milli>(t_end_sn-t_start_sn).count()/1000 << " sec\n" << endl;
	
	// inizialization of an vectors of l items
	vector<vector<vector<double> > > k;
	vector<vector<vector<double> > > c;
	k.resize(l);
	c.resize(l);

	// create clusters with connected components
	chrono::high_resolution_clock::time_point t_start_conncomp = chrono::high_resolution_clock::now();

	cout << "\nCreating " << l << " (l) clusters using connected components..." << endl;
	connectedComponents(graph, l, k, instances, minPtsInstances, c);

	chrono::high_resolution_clock::time_point t_end_conncomp = chrono::high_resolution_clock::now();
	cout << "\n" << chrono::duration<double, milli>(t_end_conncomp-t_start_conncomp).count()/1000 << " sec\n" << endl;

	// print k
	for (size_t i = 0; i < k.size(); ++i) {
		outFile << "k" << i << ":" << endl;
		for (size_t j = 0; j < k[i].size(); ++j) {
			for (size_t m = 0; m < k[i][j].size(); ++m) {
				outFile << k[i][j][m] << " ";
			}
			outFile << "\n";
		}
		outFile << endl;
	}

	// // 6. inizialization of an array of l items
	// vector<vector<vector<double> > > c;
	// c.resize(l);

	// // 7. create the clusters
	// chrono::high_resolution_clock::time_point t_start_clust = chrono::high_resolution_clock::now();

	// cout << "\nCreating " << l << " (l) clusters..." << endl;
	// findClusters(instances, l, k, c, graph);

	// chrono::high_resolution_clock::time_point t_end_clust = chrono::high_resolution_clock::now();
	// cout << "\n" << chrono::duration<double, milli>(t_end_clust-t_start_clust).count()/1000 << " sec\n" << endl;

	// print c
	outFile << "\n\n\n";
	for (size_t i = 0; i < c.size(); ++i) {
		outFile << "c" << i << ":" << endl;
        for (size_t j = 0; j < c[i].size(); ++j) {
            for (size_t m = 0; m < c[i][j].size(); ++m) {
                outFile << c[i][j][m] << " ";
            }
            outFile << "\n";
        }
        outFile << endl;
    }

	outFile.close();

	chrono::high_resolution_clock::time_point t_end_total = chrono::high_resolution_clock::now();
	cout << "\nTotal computation in " << chrono::duration<double, milli>(t_end_total-t_start_total).count()/1000 << " sec" << endl;

	return 0;
}

// FUNCTIONS:

/* function to get options for command line:
*/
void getOpt(string& filePath, int& n, double& s, double *epsilon, int& l, int& minpts, int argc, char **argv){
	int opt;

    while ((opt = getopt(argc, argv, "f:n:s:e:l:m:")) != -1) {
        switch (opt) {
		case 'f':
			filePath = optarg;
			cout << "filePath: " << filePath << endl;
			break;
        case 'n':
            n = stod(optarg);
			cout << "n: " << n << endl;
            break;
        case 's':
            s = stod(optarg);
			cout << "s: " << s << endl;
            break;
        case 'e': {
			vector<double> eps = splitString(-1, optarg);
			epsilon[0] = eps[0];
			epsilon[1] = eps[1];

			cout << "epsilon: [" << eps[0] << ", " << eps[1] << "]" << endl;
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
	instances:		vector of vector to store the istances and their datas
	dataset:		dataset file
	numInstances:	# of lines (instances) of the dataset
	vecInstances:	vector of # of the instances sampled
	currentLine:	current line pointer from the pointer of the file
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

	// reset file pointer
	dataset.seekg(0, ios::beg);

	// random sample numbers wothout rep
	instances = sampleInstances(n, numInstances, dataset);
	
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
	dataset:		dataset to use
	currentLine:	line taked in that time
	instances:		istances randomly sampled
	oldRandLine:	previous rand line generated
	range:			range of the sampling
	minRange:		min lim
	maxRange:		max lim
	loopIter:		distance from the previous rand to the next
	index:			index of the instance
*/
vector<vector<double> > sampleInstances(int n, int numLines, ifstream& dataset){
	// create a vector to semplify generation of random number without rep
	string currentLine;
	vector<vector<double> > instances;
	int oldRandLine = 0;
	int range = numLines/n;
	int minRange = 0;
	int maxRange = range;
	int index = 0;
	
	// take random elements from "lines"
	for (int i = 0; i < numLines; i += range){
		// printing the loading bar
		printLoadingBar(i, numLines);

		int randline = (rand() % (maxRange - minRange)) + minRange;
		int loopIter = randline - oldRandLine;

		// check if have many instances
		if (instances.size() == n){
			return instances;
		}

		// check if want to take the first instance
		if (randline == 0 || loopIter == 0){
			getline(dataset, currentLine);
		}

		// random sample into the dataset after current position:
		for (int j = 0; j < loopIter; j++){
			if (!getline(dataset, currentLine)) {
				cerr << "Error reading the file or the file is too short." << endl;

				// check if the current pointer of the file is in EOF
				if (dataset.eof()){
					return instances;
				}

				dataset.close();
			}
		}

		instances.push_back(splitString(index++, currentLine));

		oldRandLine = randline;
		minRange += range;
		maxRange += range;

		// check if the max limit of the range if > of the total of the lines in the dataset
		if (maxRange > numLines){
			maxRange = numLines;
		}
	}
	
	return instances;
}

/* function to split a string
	index:			index of the instance (-1 to use the function to split random string)
	line:			string to split
	ss:				string stream of the current line
	element:		single data of the readed line
	instance:		lincked vector of single instance
*/
vector<double> splitString(int index, string line){
	stringstream ss(line);
	string element;
	vector<double> instance;

	// add the index of the instance
	if (index >= 0){
		instance.push_back(index);
	}

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
	instances:	istances randomly sampled
	s:			dimension of the 2nd sample
	n:			dimension of the 1st sample
	epsilon:	range for a point
	graph:		graph of the dataset
*/
void compareWsnSample(vector<vector<double> > instances, double s, int n, double epsilon[], vector<vector<int> >& graph, int minPts, vector<int>& minPtsInstances){
	int instacesSize = instances.size();

	for (int i = 0; i < instacesSize; i++){
		// printing the loading bar
		printLoadingBar(i, instacesSize);

		int range = n/(s*n);
		int minRange = 0;
		int maxRange = range;

		// sample of s
		for (int j = 0; j < s*n; j++){
			int randline = (rand() % (maxRange - minRange)) + minRange;

			// to not randomly take the same instance
			while (i == randline){
				randline = (rand() % (maxRange - minRange)) + minRange;
			}

			// to disregard a comparison between two instances 
			if (instances[i].size() != instances[randline].size()){
				j--;
				continue;
			}

			double euclDistij = euclDist(instances[i], instances[randline]);

			// check if point are in range epsilon
			if (euclDistij >= epsilon[0] && euclDistij <= epsilon[1]){
				graph[i].insert(graph[i].begin(), randline);
			}

			minRange += range;
			maxRange += range;
			if (maxRange > n){
				maxRange = n;
			}
		}

		// check if is MinPts
		if (graph[i].size() >= minPts){
			minPtsInstances.push_back(i);
			graph[i].insert(graph[i].begin(), 1);
		} else {
			graph[i].insert(graph[i].begin(), 0);
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
*/
void connectedComponents(vector<vector<int> > graph, int l, vector<vector<vector<double> > >& k, vector<vector<double> > instances, vector<int>& minPtsInstances, vector<vector<vector<double> > >& c){
	vector<int> usedMinPtsInstances;
	vector<vector<int> > originalGraph(graph);

	// fill k with all k_i subgraphs
	for (int i = 0; i < l; i++){
		// printing the loading bar
		printLoadingBar(i, l);

		if (minPtsInstances.size() != 0){
			// add manually the 1st instance
			k[i].push_back(instances[minPtsInstances[0]]);

			// to shorten the recursive sequence
			eraseAll(graph, minPtsInstances, minPtsInstances[0]);

			// create a copy of minPtsInstance that will not be modified
			vector<int> originalMinPtsInstances(minPtsInstances);

			// add instances to the i-th cluster
			findClusters(instances, i, c, originalGraph, minPtsInstances[0]);

			findConnComp(minPtsInstances[0], graph, instances, minPtsInstances, originalMinPtsInstances, k[i], usedMinPtsInstances, c, originalGraph, i);
		}
	}
}

/* function to find the connected components:
	i:						index of the current instance
	graph:					graph where can delete 1s without influence k items
	graphCPy:				original graph
	minPtsInstances:		vector of all minPts instances
	totInstances:			size of minPtsInstances
	k:						vector where put all instaces
	usedMinPtsInstances:	vector with used minPts instances
	erased:					bool to check if current instance was used or not
	nextInst:				next instance to check after "i" (need to have like a variable because cange with deletion)
	elseMinPts:				to check if "i" joined into the else of the if
*/
void findConnComp(int i, vector<vector<int> >& graph, vector<vector<double> > instances, vector<int>& minPtsInstances, vector<int> originalMinPtsInstances, vector<vector<double> >& k, vector<int>& usedMinPtsInstances, vector<vector<vector<double> > >& c, vector<vector<int> >& originalGraph, int l){
	bool erased = false;
	int nextInst = 0;
	int j = graph[i].size()-1; // initialized to have right indexing

	while (graph[i].size() > 1){
		nextInst = graph[i][j];
		k.push_back(instances[nextInst]);

		// add instances to the i-th cluster
		findClusters(instances, l, c, originalGraph, nextInst);

		// step into only in instaces that have minPts vertices
		if (graph[nextInst][0] == 1){

			eraseAll(graph, originalMinPtsInstances, nextInst);

			// notice that all instances connection with the corrent instance was erased
			erased = true;

			// track the used minPtsIntance
			if (!isInMinPts(i, usedMinPtsInstances)){
				usedMinPtsInstances.push_back(i);

				// erase the instance from glogal minPtsInstances
				int eraseMinPts = findMinPtsInstances(i, minPtsInstances);
				minPtsInstances.erase(minPtsInstances.begin() + eraseMinPts);
			}

			findConnComp(nextInst, graph, instances, minPtsInstances, originalMinPtsInstances, k, usedMinPtsInstances, c, originalGraph, l);
		} else {
			// block the loop when the istance isn't a minPts
			eraseAll(graph, originalMinPtsInstances, nextInst);
		}

		j--;
	}

	// erase all instances connection with the current instance
	if (!erased){
		// track the used minPtsIntance
		if (!isInMinPts(i, usedMinPtsInstances)){
			usedMinPtsInstances.push_back(i);
			
			// erase the instance from glogal minPtsInstances
			int eraseMinPts = findMinPtsInstances(i, minPtsInstances);
			minPtsInstances.erase(minPtsInstances.begin() + eraseMinPts);
		}
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
void eraseAll(vector<vector<int> >& graph, vector<int> minPtsInstances, int instance){
	for (int minPtsInstance : minPtsInstances){
		for (int j = 1; j < graph[minPtsInstance].size(); j++){
			if (graph[minPtsInstance][j] == instance){
				graph[minPtsInstance].erase(graph[minPtsInstance].begin() + j);
			} else if (graph[minPtsInstance][j] < instance){
				break;
			}
		}
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

/* function to determine the clusters:
	l:		# of the l-th loop
	ki:		id of the k-th instance
*/
void findClusters(vector<vector<double> > instances, int l, vector<vector<vector<double> > >& c, vector<vector<int> >& graph, int ki){
	vector<int> allGraph(graph.size());
	iota(allGraph.begin(), allGraph.end(), 0);

	for (vector<double> instance : instances){
		if (instance[0] == ki){
			c[l].push_back(instance);
			eraseAll(graph, allGraph, instance[0]);
		} else {
			for (int i = 1; i < graph[ki].size(); i++){
				if (instance[0] == graph[ki][i]){
					c[l].push_back(instance);
					eraseAll(graph, allGraph, instance[0]);
				}
			}
		}
	}
}

/* function to print progression of the computation
	index:				index of the loop (+1 to get 100%)
	totalIterations:	total iterations of the loop
	barWidth:			width of the loading bar
	percentage:			percentage of the loading bar
	barLength:			# of "#" in the loading bar
*/
void printLoadingBar(int index, int totalIterations){
	int barWidth = 50;
	float percentage = float(index+1) / totalIterations;
	int barLength = int(percentage * barWidth);

	cout << "[";

	for (int i = 0; i < barWidth; ++i) {
		if (i < barLength) {
			cout << "#";
		} else {
			cout << ".";
		}
	}

	cout << "] " << int(percentage * 100.0) << "%" << "\r";
	cout.flush();
}