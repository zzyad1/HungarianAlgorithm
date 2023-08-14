/*
Step 0: Create a cost matrix called n × m n\times mnXm Matrix, where each element represents the n nn One of the workers is assigned to m mmThe cost of one job.

Rotate the matrix so that there are at least as many columns as there are rows, and makek = min ⁡ (n, m) k=\min(n,m)k=min ( n ,m ).

Step 1: For each row of the matrix, find the smallest element and subtract it from each element in that row. Go to step 2.

Step 2: Find zero in the result matrix (z zz). If there is no starred zero in its row and column, it willz zzStar. Repeat this operation for each element in the matrix. Go to step 3.

Step 3: Cover the columns that contain starred zeros.

If covered in total k kkColumn, the asterisked zero represents a complete unique allocation set. In this case, please go to "Finish";

Otherwise, go to step 4.

Step 4:

Find an uncovered zero and prime it.

If there are no asterisked zeros in the row containing the primed zero, go to step 5.

Otherwise, please overwrite the row and uncover the column containing the starred zero.

Continue in this way until there are no zeros left.

Save the smallest uncovered value, and then go to step 6.

Step 5: Follow the steps below to construct a series of alternating primed and starred zeros.

make z0 ​ Indicates the uncovered primed zeros found in step 4. 

make z1 exepress z0	​ The starred zero in the column (if any).​

  make z2 Express z1​ Add prime zero in the row (there will always be one).

 Continue the sequence until the zero hour that is not starred in the column where the primed zero is added.

 Cancel the asterisks of zeros in the sequence, and add stars to the apostrophes, erase all apostrophes and uncover every row in the matrix.

 Go back to step 3.

Step 6: Add the value found in step 4 to each element of each covered row, and subtract it from each element of each uncovered column.

Go back to step 4 without changing any stars, apostrophes, or overlay lines.

Done: The allocation pair is indicated by the position of an asterisked zero in the cost matrix.

ifC (i, j) C(i,j)C ( i ,j ) Is starred, it will be the same as the line i ii The associated element is assigned to the column j jj Associated elements.

In step 4, the possible situation is that there is an uncovered and primed zero.

If there is no asterisked zero in its row, the procedure goes to step 5;

if there is no uncovered zero at all, the procedure goes to step 6. 

Step 5 is the augmented path algorithm, and step 6 is to modify the cost matrix.


				HungarianAlgorithm::Solve
							|
							|
							|
							v
			HungarianAlgorithm::assignmentoptimal


*/
#include <stdlib.h>
#include <cfloat> // for DBL_MAX
#include <cmath>  // for fabs()
#include "Hungarian.h"


HungarianAlgorithm::HungarianAlgorithm(){}
HungarianAlgorithm::~HungarianAlgorithm(){}
//Vectors are same as dynamic arrays with the ability to resize itself automatically when an element is inserted or deleted,
//with their storage being handled automatically by the container.Vector elements are placed in contiguous storage
//so that they can be accessedand traversed using iterators

//Iterators            An iterator is an object (like a pointer) that points to an element inside the container.

// vectors are also containers

//begin() – Returns an iterator pointing to the first element in the vector
//end() – Returns an iterator pointing to the theoretical element that follows the last element in the vector
//rbegin() – Returns a reverse iterator pointing to the last element in the vector(reverse beginning).It moves from last to first element
//rend() – Returns a reverse iterator pointing to the theoretical element preceding the first element in the vector(considered as reverse end)
//cbegin() – Returns a constant iterator pointing to the first element in the vector.
//cend() – Returns a constant iterator pointing to the theoretical element that follows the last element in the vector.
//crbegin() – Returns a constant reverse iterator pointing to the last element in the vector(reverse beginning).It moves from last to first element
//crend() – Returns a constant reverse iterator pointing to the theoretical element preceding the first element in the vector(considered as reverse end)


//********************************************************//
// A single function wrapper for solving assignment problem.
//********************************************************//


double HungarianAlgorithm::Solve(vector <vector<double> >& DistMatrix, vector<int>& Assignment)//
{
	//Vector of Vectors is a two-dimensional vector with a variable number of rows where each row is vector.
	//Each index of vector stores a vector which can be traversed and accessed using iterators.
	//It is similar to an Array of Vectors but with dynamic properties.
	unsigned int nRows = DistMatrix.size();

	unsigned int nCols = DistMatrix[0].size();



	//The keyword unsigned is a data type specifier,
	//which only represents non-negative integers i.e. positive numbers and zero

	// DistMatrix construction distMatrixIn, the input matrix is ​​first stored in MATLAB column.

	// DistMatrix Compute a symmetric matrix of distances (or similarities) between the rows or columns of a matrix

	double *distMatrixIn = new double[nRows * nCols];

	int *assignment = new int[nRows];

	//construct a new array
	// full array and array for assignment each one has a pointer

	double cost = 0.0;
	// initialization for the variable cost


	// Fill in the distMatrixIn. Mind the index is "i + nRows * j".
	// Here the cost matrix of size MxN is defined as a double precision array of N*M elements. 
	// In the solving functions matrices are seem  to be saved MATLAB-internally in row-order.
	// (i.e. the matrix [1 2; 3 4] will be stored as a vector [1 3 2 4], NOT [1 2 3 4]).


	//HungarianAlgorithm::assignmentoptimal is the main implementation.
	//	The calculated result is given Assignment , so we are directly passing it
	for (unsigned int i = 0; i < nRows; i++)
		for (unsigned int j = 0; j < nCols; j++)
			distMatrixIn[i + nRows * j] = DistMatrix[i][j];
	// numbering the blocks of the  given matrix i.e[0,1 0,2 0,3]
	


	// call solving function

	assignmentoptimal(assignment, &cost, distMatrixIn, nRows, nCols);
	//The parameter declaration int & cost defines a as a reference to a int

	Assignment.clear();
	//clear() – It is used to remove all the elements of the vector container
	

	  // Inserting elements into vector
	for (unsigned int r = 0; r < nRows; r++)
		Assignment.push_back(assignment[r]);


	//where r refers to the element
	//to be added in the back of the vector
	delete[] distMatrixIn;
	delete[] assignment;
	return cost;
	//erasing vectors from the memory and returning the cost
}


//********************************************************//
// Solve optimal solution for assignment problem using Munkres algorithm, also known as Hungarian Algorithm.
//********************************************************//



void HungarianAlgorithm::assignmentoptimal(int *assignment, double *cost, double *distMatrixIn, int nOfRows, int nOfColumns)
{

	//the * is Declaration of a pointer
	// so that *assignment is a pointer to int...
	double *distMatrix, *distMatrixTemp, *distMatrixEnd, *columnEnd, value, minValue;
	bool *coveredColumns, *coveredRows, *starMatrix, *newStarMatrix, *primeMatrix;
	int nOfElements, minDim, row, col;
	//assignment is The matching result for each trajectory.

	/* initialization */
	*cost = 0;

	// Generate distance matrix 
		// and check matrix elements positiveness :)
	for (row = 0; row<nOfRows; row++)
		assignment[row] = -1;

	/* generate working copy of distance Matrix */
	/* check if all matrix elements are positive */



	// Total elements number
	nOfElements = nOfRows * nOfColumns;

	// Memory allocation
	distMatrix = (double *)malloc(nOfElements * sizeof(double));

	//The malloc function allocates a memory block of at least size bytes.
	//The block may be larger than size bytes because of the space that's required for alignment and maintenance information.
	
	
	distMatrixEnd = distMatrix + nOfElements;

	// distMatrixEndIs the end of the matrix.  simply its a Pointer to last element in the matrix
	

	// It is distMatrixIn constructed distMatrix to prevent data from being overwritten.Frustrated,
	// just newStarMatrixpass in HungarianAlgorithm::step4 all the way to use it.
		
	{
		value = distMatrixIn[row];
		if (value < 0)
			cerr << "All matrix elements have to be non-negative." << endl;
		distMatrix[row] = value;
	}


	/* memory allocation */
	coveredColumns = (bool *)calloc(nOfColumns, sizeof(bool));
	coveredRows = (bool *)calloc(nOfRows, sizeof(bool));
	starMatrix = (bool *)calloc(nOfElements, sizeof(bool));
	primeMatrix = (bool *)calloc(nOfElements, sizeof(bool));
	newStarMatrix = (bool *)calloc(nOfElements, sizeof(bool)); /* used in step4 */

	//The calloc function allocates storage space for an array of number elements,
	// each of length size bytes. Each element is initialized to 0.


	/* preliminary steps */

	//	Start with a small number of directions and look for the minimum distance.

	if (nOfRows <= nOfColumns)
	{
		minDim = nOfRows;

		for (row = 0; row<nOfRows; row++)
		{
			/* find the smallest element in the row */
			distMatrixTemp = distMatrix + row;
			minValue = *distMatrixTemp;
			// assuming that the element we pointed to is the minimum value
			distMatrixTemp += nOfRows;
			// moving the pointer to the next one
			while (distMatrixTemp < distMatrixEnd) 	// first comparison
			{

				value = *distMatrixTemp;
				
				if (value < minValue) // second comparison
					minValue = value;
				distMatrixTemp += nOfRows;

				// looping until it finds the smallest value 
			}

			/* subtract the smallest element from each element of the row */
			distMatrixTemp = distMatrix + row;
			while (distMatrixTemp < distMatrixEnd)
			{
				*distMatrixTemp -= minValue;
				distMatrixTemp += nOfRows;
				// looping and subtracting smallest element from all elements in this row
			}
		}

		/* Steps 1 and 2a */
	//Modify directly distMatrix, subtract the minimum value from each row(column).
		for (row = 0; row<nOfRows; row++)
			for (col = 0; col<nOfColumns; col++)
				if (fabs(distMatrix[row + nOfRows*col]) < DBL_EPSILON)
					if (!coveredColumns[col])
					{
						//starMatrix : Records the current match.
						starMatrix[row + nOfRows*col] = true;
						coveredColumns[col] = true;
						break;
					
						//	If there is no conflict(one - to - many), set the coveredColumns corresponding position to true,
						

						//The fabs() function  returns the absolute value of the argument.
					
						// fabs() Return Value : the absolute value of distMatrix[row + nOfRows*col i.e. |distMatrix[row + nOfRows*col|
					
						
					}
	}

	else /* if(nOfRows > nOfColumns) */
	{

		//	Search by column will face the problem of discontinuous access.
		//so we will access a two-dimensional array
		minDim = nOfColumns;

		
		
		
		for (col = 0; col<nOfColumns; col++)
		{
			/* find the smallest element in the column */

			//distMatrixTemp : The name is a little weird, just a pointer to the element.
			distMatrixTemp = distMatrix + nOfRows*col;
			columnEnd = distMatrixTemp + nOfRows;

			minValue = *distMatrixTemp++;
			while (distMatrixTemp < columnEnd)
			{
				value = *distMatrixTemp++;
				if (value < minValue)
					minValue = value;
			}

			/* subtract the smallest element from each element of the column */
			distMatrixTemp = distMatrix + nOfRows*col;
			while (distMatrixTemp < columnEnd)
				*distMatrixTemp++ -= minValue;
		}
		/* Steps 1 and 2a */

		
		for (col = 0; col<nOfColumns; col++)
			for (row = 0; row<nOfRows; row++)
				if (fabs(distMatrix[row + nOfRows*col]) < DBL_EPSILON)
					if (!coveredRows[row])
					{
						starMatrix[row + nOfRows*col] = true;
						coveredColumns[col] = true;
						coveredRows[row] = true;
						break;
						




						
					}
		
		for (row = 0; row<nOfRows; row++)
			coveredRows[row] = false;
		// we set coveredRows[row] to when there are more rows false.

	}
	

	/* move to step 2b */
// counts the number of covered columns.


	step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

	/* 
			remove invalid assignments
		and also calculates the cost of the allocation plan.

 */
	computeassignmentcost(assignment, cost, distMatrixIn, nOfRows);

	/* free allocated memory */
	free(distMatrix);
	free(coveredColumns);
	free(coveredRows);
	free(starMatrix);
	free(primeMatrix);
	free(newStarMatrix);

	return;
}

/********************************************************/
void HungarianAlgorithm::buildassignmentvector(int *assignment, bool *starMatrix, int nOfRows, int nOfColumns)
{
//buildassignmentvector gets the result by asterisk.
	int row, col;

	for (row = 0; row<nOfRows; row++)
		for (col = 0; col<nOfColumns; col++)
			if (starMatrix[row + nOfRows*col])
			{
#ifdef ONE_INDEXING
				assignment[row] = col + 1; /* MATLAB-Indexing */
#else
				assignment[row] = col;
#endif
				break;
			}
}
//If the starMatrix corresponding position is marked, the matching result is recorded in the assignmentarray.
/********************************************************/
void HungarianAlgorithm::computeassignmentcost(int *assignment, double *cost, double *distMatrix, int nOfRows)
{
	int row, col;

	for (row = 0; row<nOfRows; row++)
	{
		col = assignment[row];
		if (col >= 0)
			*cost += distMatrix[row + nOfRows*col];
	}

}

/********************************************************/
void HungarianAlgorithm::step2a(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	bool *starMatrixTemp, *columnEnd;
	int col;

	/* cover every column containing a starred zero */
	for (col = 0; col<nOfColumns; col++)
	{
		starMatrixTemp = starMatrix + nOfRows*col;
		columnEnd = starMatrixTemp + nOfRows;
		while (starMatrixTemp < columnEnd){
			if (*starMatrixTemp++)
			{
				coveredColumns[col] = true;
				break;
			}
		}
	}

	/* move to step 3 */
	step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/
/*
			 ____ step2b ________
			|					 |
			|					 |
		finished				 |
			|					 |
			v				   	 v
	buildassignmentvector		step3


*/
/********************************************************/

//Count the number of covered columns. we just Determine whether it is complete or not

void HungarianAlgorithm::step2b(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	int col, nOfCoveredColumns;

	/* count covered columns */
	nOfCoveredColumns = 0;
	for (col = 0; col<nOfColumns; col++)
		if (coveredColumns[col])
			nOfCoveredColumns++;

	if (nOfCoveredColumns == minDim)
	{
		/* algorithm finished */
		buildassignmentvector(assignment, starMatrix, nOfRows, nOfColumns);
	}
	else
	{
		/* move to step 3 */
		step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
	}

}
/********************************************************/
/*
			 ____ step3 ________
			|					 |
			|					 |
			|					 |
			|					 |
			v				   	 v
		  step4					step5

//The indentation structure is too deep.

*/
/********************************************************/
void HungarianAlgorithm::step3(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	bool zerosFound;
	int row, col, starCol;

	zerosFound = true;
	while (zerosFound)
	{
		//If there is no asterisked zero in the line containing the primed zero, call HungarianAlgorithm::step4
		zerosFound = false;
		for (col = 0; col<nOfColumns; col++)
			if (!coveredColumns[col])
				for (row = 0; row<nOfRows; row++)
					if ((!coveredRows[row]) && (fabs(distMatrix[row + nOfRows*col]) < DBL_EPSILON))
					{
						//Find an uncovered zero and prime it. 
						/* prime zero */
						primeMatrix[row + nOfRows*col] = true;

						/* find starred zero in current row */
						for (starCol = 0; starCol<nOfColumns; starCol++)
							if (starMatrix[row + nOfRows*starCol])
								break;

						if (starCol == nOfColumns) /* no starred zero found */
						{

						
							/* move to step 4 */
							step4(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim, row, col);
							return;
						}
						else
						{
							coveredRows[row] = true;
							coveredColumns[starCol] = false;
							zerosFound = true;
							break;
						}
					}
	}
//Otherwise, overwrite the row and uncover the column containing the starred zero. Continue in this way until there are no zeros left.
//Save the smallest uncovered value, and then call HungarianAlgorithm::step5 .

	/* move to step 5 */
	step5(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/
/*
The rank index is specified.
newStarMatrix Copy the whole first starMatrix, which is slightly inefficient.
*/
void HungarianAlgorithm::step4(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim, int row, int col)
{
	int n, starRow, starCol, primeRow, primeCol;
	int nOfElements = nOfRows*nOfColumns;

	/* generate temporary copy of starMatrix */
	for (n = 0; n<nOfElements; n++)
		newStarMatrix[n] = starMatrix[n];

	/* star current zero */
	newStarMatrix[row + nOfRows*col] = true;

	/* find starred zero in current column */
	starCol = col;
	for (starRow = 0; starRow<nOfRows; starRow++)
		if (starMatrix[starRow + nOfRows*starCol])
			break;

	while (starRow<nOfRows)
	{
		/* unstar the starred zero */
		newStarMatrix[starRow + nOfRows*starCol] = false;

		/* find primed zero in current row */
		primeRow = starRow;
		for (primeCol = 0; primeCol<nOfColumns; primeCol++)
			if (primeMatrix[primeRow + nOfRows*primeCol])
				break;

		/* star the primed zero */
		newStarMatrix[primeRow + nOfRows*primeCol] = true;

		/* find starred zero in current column */
		starCol = primeCol;
		for (starRow = 0; starRow<nOfRows; starRow++)
			if (starMatrix[starRow + nOfRows*starCol])
				break;
	}

	/* use temporary copy as new starMatrix */
	/* delete all primes, uncover all rows */
	for (n = 0; n<nOfElements; n++)
	{
		primeMatrix[n] = false;
		starMatrix[n] = newStarMatrix[n];
	}
	for (n = 0; n<nOfRows; n++)
		coveredRows[n] = false;

	/* move to step 2a */
	step2a(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/

void HungarianAlgorithm::step5(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	double h, value;
	int row, col;

	/* find smallest uncovered element h */
	h = DBL_MAX;
	for (row = 0; row<nOfRows; row++)
		if (!coveredRows[row])
			for (col = 0; col<nOfColumns; col++)
				if (!coveredColumns[col])
				{
					value = distMatrix[row + nOfRows*col];
					if (value < h)
						h = value;
				}

	/* add h to each covered row */
	for (row = 0; row<nOfRows; row++)
		if (coveredRows[row])
			for (col = 0; col<nOfColumns; col++)
				distMatrix[row + nOfRows*col] += h;

	/* subtract h from each uncovered column */
	for (col = 0; col<nOfColumns; col++)
		if (!coveredColumns[col])
			for (row = 0; row<nOfRows; row++)
				distMatrix[row + nOfRows*col] -= h;

	/* move to step 3 */
	step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}
/*


				HungarianAlgorithm::step2a
							|
							|
							|
							v
				HungarianAlgorithm::step2b

Iterate through each column and cover the column if there are starred zeros in the element


*/