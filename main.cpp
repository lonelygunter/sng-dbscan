#include <iostream>
#include <fstream>	// file io
#include <string>
#include <getopt.h>	// getopt
#include <chrono>

#include "headers/readDataset.h"
#include "headers/compareWsnSample.h"
#include "headers/connectedComponents.h"
#include "headers/conn_comp_igraph.h"

using namespace std;

// declare functions:
void getOpt(string& filePath, int& n, double& s, double *epsilon, int& minpts, int argc, char **argv);


/*
	instances:	istances randomly sampled
	n:			dimension of the 1st sample
	s:			dimension of the 2nd sample
	epsilon:	range for a point
	minpts:		min # of points in the epsilon range of a point
	filePath:	path of the dataset
	graph:		graph of the dataset
	k:			subgraph induced by vertices of degree at least MinPts
	c:			vector of the clusters
*/
int main(int argc, char **argv){
	vector<vector<double> > instances;
	vector<int> minPtsInstances;

	// inizialize parameters
	int n = 0;
	double s = 0;
	double epsilon[2];
	int minPts = 0;
	string filePath = "datasets/";

	ofstream outFile("clustering.txt");

	// getopt:
	getOpt(filePath, n, s, epsilon, minPts, argc, argv);

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

	cout << "\nSampling " << static_cast<int>(s*n) << " (sn) instances..." << endl;
	compareWsnSample(instances, s, n, epsilon, graph, minPts, minPtsInstances);

	chrono::high_resolution_clock::time_point t_end_sn = chrono::high_resolution_clock::now();
	cout << "\n" << chrono::duration<double, milli>(t_end_sn-t_start_sn).count()/1000 << " sec\n" << endl;

	// create clusters with connected components
	chrono::high_resolution_clock::time_point t_start_conncomp = chrono::high_resolution_clock::now();

	cout << "\nCreating clusters using connected components..." << endl;
	// connectedComponents(graph, l, k, instances, minPtsInstances, c);

	// calculate connected components:
	vector<vector<int> > conncomps;
	igraph_integer_t num_components = 2;

	conncomps = conn_comp_igraph(graph, n, num_components);

	// adact graph for clustering
	for (int i = 0; i < graph.size(); i++) {
		graph[i].erase(graph[i].begin()); // Erase the first item that indicate if is >= Minpts
	}

	// create clusters
	vector<vector<vector<double> > > c;
	c.resize(num_components);

	    vector<int> allConnComps(conncomps.size());
    for (int i = 0; i < conncomps.size(); ++i) {
        allConnComps.push_back(i);
    }

    for (int component = 0; component < num_components; component++){
        // cout << "component: " << component << endl;

        for (int i = 0; i < instances.size(); i++){
            // cout << "instance: " << i << endl;

            for (int j = 0; j < conncomps.size(); j++){
                // cout << "conncomps: " << j << endl;

                if (!conncomps[j].empty() && conncomps[j][0] == component){
                    // cout << "if " << conncomps[j][0] << " == " << component << endl;

                    if (graph[i] == conncomps[j] && graph[i].size() > 1){
                        // cout << "if conncomps[j] == conncomps[j]" << endl;
                        // cout << "added instance " << i << " in cluster " << component << endl;
                        c[component].push_back(instances[i]);
                        eraseAll(conncomps, allConnComps, i);

                    } else {

                        for (int w = 1; w < conncomps[j].size(); w++){
                            // cout << "conncomps " << j << " link: " << conncomps[j].size() << endl;

                            if (i == conncomps[j][w]){
                                // cout << "if " << i << " == " << conncomps[j][w] << endl;
                                // cout << "added instance " << i << " in cluster " << component << endl;
                                c[component].push_back(instances[i]);
                                eraseAll(conncomps, allConnComps, i);
                            }
                        }
                    }
                }
            }
        }
    }

	chrono::high_resolution_clock::time_point t_end_conncomp = chrono::high_resolution_clock::now();
	cout << "\n" << chrono::duration<double, milli>(t_end_conncomp-t_start_conncomp).count()/1000 << " sec\n" << endl;

	chrono::high_resolution_clock::time_point t_end_total = chrono::high_resolution_clock::now();
	cout << "\nTotal computation in " << chrono::duration<double, milli>(t_end_total-t_start_total).count()/1000 << " sec" << endl;

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

	return 0;
}

// FUNCTIONS:

/* function to get options for command line:
*/
void getOpt(string& filePath, int& n, double& s, double *epsilon, int& minpts, int argc, char **argv){

	int opt;

    while ((opt = getopt(argc, argv, "f:n:s:e:m:")) != -1) {
        switch (opt) {
		case 'f':
			filePath += optarg;
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
