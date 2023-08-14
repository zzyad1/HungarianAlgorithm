#include <iostream>
#include "Hungarian.h"
#include "testMain.h"


int main(void)
{
	//initializing costmatrix variable with the type vector of vectors matrix
	vector< vector<double> > costMatrix = { { 82, 83, 69, 92 }, 
										  { 77, 37, 49, 92}, 
										  { 11, 69,5, 86}, 
										  { 8, 9, 98, 23 } };

	HungarianAlgorithm HungAlgo; // initializing hungalgo object from hungarianalgorithm class
	vector<int> assignment;
	//initializing assignment variable with the type vector matrix 
	double cost = HungAlgo.Solve(costMatrix, assignment);  //calling solve function with passing costmatrix and assignment parameters to it

	for (unsigned int x = 0; x < costMatrix.size(); x++)
		std::cout << x << "," << assignment[x] << "\t"; // printing the the positions of zeros

	std::cout << "\ncost: " << cost << std::endl; // printing the cost

	return 0;
}
