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
	n1 = atoi(line.c_str());
	line.erase(0, line.find(' ') + 1);
	m1 = atoi(line.c_str());
	mat1 = new double[n1*m1];

	double* cur = mat1;
	while(std::getline(firstMatFile, line)){
		(cur++)[0] = atof(line.c_str());
	}
	firstMatFile.close();

	std::ifstream secondMatFile(argv[2]);
	std::getline(secondMatFile, line);
	n2 = atoi(line.c_str());
	line.erase(0, line.find(' ') + 1);
	m2 = atoi(line.c_str());
	mat2 = new double[n2*m2];
	
	cur = mat2;
	while(std::getline(secondMatFile, line)){
		(cur++)[0] = atof(line.c_str());
	}
	secondMatFile.close();

	if(n1 != m2){
		std::cerr << "Die Zeilenanzahl der ersten Matrix (hier " << n1 << ") muss gleich der Spaltenanzahl der zweiten Matrix (hier " << m2 << ") sein!" << std::endl;
		return 1;
	}


   siwir::Timer timer;

/*
#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("matmult");
#endif
*/

	// bitte hier koten
   std::cout << "Hallo :)\n";

/*  
#ifdef USE_LIKWID
   likwid_markerStopRegion( "matmult" );
#endif  
*/

   double time = timer.elapsed();

   
    std::cout << "Calculation took " << time << " seconds\n";  
}
