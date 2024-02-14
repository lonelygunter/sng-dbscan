#ifndef CONN_COMP_IGRAPH_H
#define CONN_COMP_IGRAPH_H

#include <igraph/igraph.h>
#include <iostream>
#include <vector>
#include <utility>
#include <fstream>

using namespace std;

// function to compute connected components of an igraph
vector<vector<int> > conn_comp_igraph(vector<vector<int> > dataset, int n, igraph_integer_t &num_components);

// function to adapt the graph in a proper data structure
bool createGraphFromDataset(const vector<vector<int> >& dataset, igraph_t* graph, int n);

#endif // CONN_COMP_IGRAPH_H
