#include "compareWsnSample.h"

#include "loadingBar.h"

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
		if (graph[i].size() >= static_cast<size_t>(minPts)){
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

	for (size_t j = 0; j < point1.size(); j++){
		euclDistij += pow(point1[j] - point2[j], 2);
	}
	
	euclDistij = sqrt(euclDistij);

	return euclDistij;
}