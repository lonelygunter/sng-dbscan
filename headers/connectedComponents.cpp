#include "connectedComponents.h"

#include "loadingBar.h"

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

/* function to shorten the recursive function fincConnComp:
	instance:		number of the column instance of the graph
*/
void eraseAll(vector<vector<int> >& graph, vector<int> minPtsInstances, int instance){
	for (int minPtsInstance : minPtsInstances){
		for (size_t j = 1; j < graph[minPtsInstance].size(); j++){
			if (graph[minPtsInstance][j] == instance){
				graph[minPtsInstance].erase(graph[minPtsInstance].begin() + j);
			} else if (graph[minPtsInstance][j] < instance){
				break;
			}
		}
	}
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
			for (size_t i = 1; i < graph[ki].size(); i++){
				if (instance[0] == graph[ki][i]){
					c[l].push_back(instance);
					eraseAll(graph, allGraph, instance[0]);
				}
			}
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
	for (size_t i = 0; i < minPtsInstances.size(); i++){
		if (minPtsInstances[i] == instance){
			return true;
		}
	}
	
	return false;
}

/* function to find the id of an instance into minPtsInstances:
	instance:		number of ther instance to find
*/
int findMinPtsInstances(int instance, vector<int> minPtsInstances){
	for (size_t i = 0; i < minPtsInstances.size(); i++){
		if (minPtsInstances[i] == instance){
			return i;
		}
	}
	
	return -1;
}