#ifndef COMPARE_WSN_SAMPLE_H
#define COMPARE_WSN_SAMPLE_H

#include <vector>
#include <cmath>	// math

using namespace std;


// function that take one point of the "n" sample to compare it with sampled sn points
void compareWsnSample(vector<vector<double> > instances, double s, int n, double epsilon[], vector<vector<int> >& graph, int minPts, vector<int>& minPtsInstances);

// function that calculate the Euclidean distance between two points
double euclDist(vector<double> point1, vector<double> point2);

#endif // COMPARE_WSN_SAMPLE_H
