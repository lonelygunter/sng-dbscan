#include "readDataset.h"

#include "loadingBar.h"

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
		if (instances.size() == static_cast<size_t>(n)){
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