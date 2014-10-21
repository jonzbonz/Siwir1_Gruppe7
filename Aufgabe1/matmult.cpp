/*==========================================================================
 * Performance so far:														|
 * Naive Implementation: 0.344938s(512x512), 2.93396s(1024x1024)			|
 * TODO why does it take over 6 seconds for the make perf1024 command to	|
 * 		finish?																|
 *==========================================================================
 */

#include "Timer.h"

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>


/* 
#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}
#endif
*/


int main(int argc, char* argv[])
{
	if(argc != 4){
		std::cerr << "Usage: ./matmult A.in B.in C.out" << std::endl;
		return 1;
	}

	//Reading input / creating matrices
	double* mat1 = 0;
	double* mat2 = 0;
	int n1, m1, n2, m2;

	std::string line;
	
	std::ifstream firstMatFile(argv[1]);
	std::getline(firstMatFile, line);
	m1 = atoi(line.c_str());
	line.erase(0, line.find(' ') + 1);
	n1 = atoi(line.c_str());
	mat1 = new double[n1*m1];

	double* cur = mat1;
	while(std::getline(firstMatFile, line)){
		(cur++)[0] = atof(line.c_str());
	}
	firstMatFile.close();

	std::ifstream secondMatFile(argv[2]);
	std::getline(secondMatFile, line);
	m2 = atoi(line.c_str());
	line.erase(0, line.find(' ') + 1);
	n2 = atoi(line.c_str());
	mat2 = new double[n2*m2];
	
	cur = mat2;
	while(std::getline(secondMatFile, line)){
		(cur++)[0] = atof(line.c_str());
	}
	secondMatFile.close();

	if(n1 != m2){ 
		std::cerr << "Die Spaltenanzahl der ersten Matrix (hier " << n1 << ") muss gleich der Zeilenanzahl der zweiten Matrix (hier " << m2 << ") sein!" << std::endl;
		return 1;
	}

	int n3 = n2;
	int m3 = m1;
	double* ergmat = new double[n3*m3];

/* 
#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("matmult");
#endif
*/

	siwir::Timer timer;

	// bitte hier koten
	for(int i = 0; i  < m3; ++i){					// ueber die Zeilen von ergmat
		for(int j = 0; j  < n3; ++j){				// ueber die Spalten von ergmat
			for(int k = 0;  k < n1/*bzw m2*/; ++k){	// ueber die Spalten von mat1 / Zeilen von mat2
				ergmat[i*n3 + j] += mat1[i*n1 + k] * mat2[k*n2 + j];
			}
		}
	}

   	double time = timer.elapsed();
	

/*   
#ifdef USE_LIKWID
	likwid_markerStopRegion( "matmult" );
#endif  
*/

	
	std::ofstream outfile(argv[3]);
	outfile << m3 << " " << n3 << std::endl;

	for(int i = 0; i < n3*m3; ++i){
		outfile << ergmat[i] << std::endl;
	}

	outfile.close();
   
    std::cout << "Calculation took " << time << " seconds\n";  
}
