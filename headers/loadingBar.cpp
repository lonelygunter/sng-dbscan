#include "loadingBar.h"

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