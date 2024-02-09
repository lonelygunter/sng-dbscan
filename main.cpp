#include <iostream>
#include <fstream>	// file io
#include <string>
#include <getopt.h>	// getopt
#include <chrono>

#include "headers/readDataset.h"
#include "headers/compareWsnSample.h"
#include "headers/connectedComponents.h"

using namespace std;

// declare functions:
void getOpt(string& filePath, int& n, double& s, double *epsilon, int& l, int& minpts, int argc, char **argv);


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

	// australian-1M
	// int n = 3000;
	// double s = 0.4;
	// double epsilon[2] = {15, 21};
	// int l = 10;
	// int minPts = 2;
	// string filePath = "datasets/australian-1M.arff";

	// banck-8k
	// int n = 3000;
	// double s = 0.4;
	// double epsilon[2] = {15, 21};
	// int l = 10;
	// int minPts = 2;
	// string filePath = "datasets/banck-8k.arff";

	// ionosphere-351
	// int n = 3000;
	// double s = 0.4;
	// double epsilon[2] = {15, 21};
	// int l = 10;
	// int minPts = 2;
	// string filePath = "datasets/ionosphere-351.data";

	// iris-150
	int n = 150;
	double s = 0.5;
	double epsilon[2] = {0.2, 2.9};
	int l = 10;
	int minPts = 2;
	string filePath = "datasets/iris-150.data";

	// kc2-522
	// int n = 3000;
	// double s = 0.4;
	// double epsilon[2] = {15, 21};
	// int l = 10;
	// int minPts = 2;
	// string filePath = "datasets/kc2-522.arff";

	// libras-360
	// int n = 3000;
	// double s = 0.4;
	// double epsilon[2] = {15, 21};
	// int l = 10;
	// int minPts = 2;
	// string filePath = "datasets/libras-360.data";

	// mozilla-15k
	// int n = 3000;
	// double s = 0.4;
	// double epsilon[2] = {15, 21};
	// int l = 10;
	// int minPts = 2;
	// string filePath = "datasets/mozilla-15k.arff";

	// pageblocks-5k
	// int n = 3000;
	// double s = 0.4;
	// double epsilon[2] = {15, 21};
	// int l = 10;
	// int minPts = 2;
	// string filePath = "datasets/pageblocks-5k.data";

	// satimage-1M
	// int n = 3000;
	// double s = 0.4;
	// double epsilon[2] = {15, 21};
	// int l = 10;
	// int minPts = 2;
	// string filePath = "datasets/satimage-1M.arff";

	// tokyo-959
	// int n = 3000;
	// double s = 0.4;
	// double epsilon[2] = {15, 21};
	// int l = 10;
	// int minPts = 2;
	// string filePath = "datasets/tokyo-959.arff";

	// vehicle-846
	// int n = 3000;
	// double s = 0.4;
	// double epsilon[2] = {15, 21};
	// int l = 10;
	// int minPts = 2;
	// string filePath = "datasets/vehicle-846.arff";

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
