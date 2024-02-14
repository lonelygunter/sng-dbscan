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
		// cout << "i " << i << " instacesSize " << instacesSize << endl;
		// printing the loading bar
		printLoadingBar(i, instacesSize);

		int range = n/(s*n);
		int minRange = 0;
		int maxRange = range;

		// to not have infinite loop into the while
		if (range <= 1){
			range = 2;
			maxRange = range;
		}

		// sample of s
		for (int j = 0; j < static_cast<int>(s*n); j++){
			// cout << "j " << j << " s*n " << static_cast<int>(s*n) << endl;
			int randline = (rand() % (maxRange - minRange)) + minRange;
			// cout << "maxRange " << maxRange << " minRange " << minRange << endl;
			// cout << "randline " << randline << endl;
			

			// to not randomly take the same instance
			while (i == randline){
				randline = (rand() % (maxRange - minRange)) + minRange;
			}

			// to disregard a comparison between two instances 
			if (instances[i].size() != instances[randline].size()){
				// cout << instances[i].size() << " != " << instances[randline].size() << endl;
				continue;
			}

			double calcDistij = calcDist(instances[i], instances[randline]);

			// check if point are in range epsilon
			if (calcDistij >= epsilon[0] && calcDistij <= epsilon[1]){
				graph[i].insert(graph[i].begin(), randline);
			}

			minRange += range;
			maxRange += range;
			// cout << "POST: maxRange " << maxRange << " minRange " << minRange << endl;
			if (maxRange > n){
				maxRange = n;
				minRange = n - range;
				// cout << "IF: maxRange " << maxRange << " minRange " << minRange << endl;
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
	calcDistij:		Euclidean distance
*/
// double calcDist(vector<double> point1, vector<double> point2){
// 	double calcDistij = 0;

// 	for (size_t j = 1; j < point1.size(); j++){
// 		calcDistij += pow(point1[j] - point2[j], 2);
// 	}
	
// 	calcDistij = sqrt(calcDistij);

// 	return calcDistij;
// }
double calcDist(vector<double> point1, vector<double> point2){
	double calcDistij = 0;

	for (size_t j = 1; j < point1.size(); j++){
		calcDistij += pow(point1[j ] - point2[j], 2);
	}

	return calcDistij;
}