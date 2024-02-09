#ifndef CONNECTED_COMPONENTS_H
#define CONNECTED_COMPONENTS_H

#include <vector>
#include <numeric>

using namespace std;

// function to create the connected components
void connectedComponents(vector<vector<int> > graph, int l, vector<vector<vector<double> > >& k, vector<vector<double> > instances, vector<int>& minPtsInstances, vector<vector<vector<double> > >& c);

// function to shorten the recursive function fincConnComp
void eraseAll(vector<vector<int> >& graph, vector<int> minPtsInstances, int instance);

// function to determine the clusters
void findClusters(vector<vector<double> > instances, int l, vector<vector<vector<double> > >& c, vector<vector<int> >& graph, int ki);

// function to find the connected components
void findConnComp(int i, vector<vector<int> >& graph, vector<vector<double> > instances, vector<int>& minPtsInstances, vector<int> originalMinPtsInstances, vector<vector<double> >& k, vector<int>& usedMinPtsInstances, vector<vector<vector<double> > >& c, vector<vector<int> >& originalGraph, int l);

// function to check if an instance is in minPts range
bool isInMinPts(int instance, vector<int> minPtsInstances);

// function to find the id of an instance into minPtsInstances
int findMinPtsInstances(int instance, vector<int> minPtsInstances);

#endif // CONNECTED_COMPONENTS_H
