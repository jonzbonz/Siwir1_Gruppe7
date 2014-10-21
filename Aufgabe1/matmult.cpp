/*==========================================================================
 * TODO Enter todo's														|
 *==========================================================================
 */

#include "Timer.h"

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


 
#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}
#endif



int main(int argc, char* argv[])
{
	if(argc != 4){
		std::cerr << "Usage: ./matmult A.in B.in C.out" << std::endl;
		return 1;
	}

	//Reading input / creating matrices
	double* matA = 0;
	double* matB = 0;
	int na, ma, nb, mb;

	std::string line;
	
	std::ifstream firstMatFile(argv[1]);
	std::getline(firstMatFile, line);
	ma = atoi(line.c_str());
	line.erase(0, line.find(' ') + 1);
	na = atoi(line.c_str());

	std::ifstream secondMatFile(argv[2]);
	std::getline(secondMatFile, line);
	mb = atoi(line.c_str());
	line.erase(0, line.find(' ') + 1);
	nb = atoi(line.c_str());
	
	bool useStrassen = false;
	int nn = 1;
	if(na == ma && nb == mb && n2*m1 >= 0){ //TODO gute grenze finden
		std::cout << "using strassen" << std::endl;
		useStrassen = true;
		
		while(nn < na) nn = nn << 1;

		matA = new double[nn*nn];
		memset(matA, 0, nn*nn*sizeof(double));
		matB = new double[nn*nn];
		memset(matB, 0, nn*nn*sizeof(double));
	} else {
		matA = new double[na*ma];
		std::cout << na << "x" << ma << " Matrix erstellt" << std::endl;
		matB = new double[nb*mb];
		std::cout << nb << "x" << mb << " Matrix erstellt" << std::endl;
	}
	

	for(int i = 0; i < ma; ++i){
		for(int j = 0; j < na; ++j){
			if(!std::getline(firstMatFile, line)){
				std::cerr << "Unexpected end of file!" << std::endl;
				return 1;
			}
			matA[i*na + j] = atof(line.c_str());
		}
	}
	firstMatFile.close();

	for(int i = 0; i < mb; ++i){
		for(int j = 0; j < nb; ++j){
			if(!std::getline(secondMatFile, line)){
				std::cerr << "Unexpected end of file!" << std::endl;
				return 1;
			}
			matB[i*nb + j] = atof(line.c_str());
		}
	}
	secondMatFile.close();

	if(na != mb){ 
		std::cerr << "Die Spaltenanzahl der ersten Matrix (hier " << na << ") muss gleich der Zeilenanzahl der zweiten Matrix (hier " << mb << ") sein!" << std::endl;
		return 1;
	}

	int nc = nb;
	int mc = ma;
	double* matC = new double[nc*mc];

 
#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("matmult");
#endif


	siwir::Timer timer;


	if(!useStrassen){
		std::cout << na << ma << nb << mb << nc << mc << std::endl;


		for(int j = 0; j  < nc; ++j){						// ueber die Spalten von matC
			for(int i = 0; i  < mc; ++i){					// ueber die Zeilen von matC
				for(int k = 0;  k < na/*bzw mb*/; ++k){		// ueber die Spalten von matA / Zeilen von matB
					matC[i*nc + j] += matA[i*na + k] * matB[k*nb + j];
				}
			}
		}
	} else {
		
		//Strassen-Algorithmus - nur fuer grosse und quadratische matrizen
		
		std::cout << "Strassen not implemented yet!!!" << std::endl;	
		
		int offset = nn*nn/4;
		double* M = new double[offset*7];

		for(int j = 0; j  < nn/2; ++j){
			for(int i = 0; i  < nn/2; ++i){
				
				for(int k = 0; k < nn/2; k++){
					//TODO
				}

			}
		}


	}




   	double time = timer.elapsed();
	

   
#ifdef USE_LIKWID
	likwid_markerStopRegion("matmult");
	likwid_markerClose();
#endif  


	
	std::ofstream outfile(argv[3]);
	outfile << mc << " " << nc << std::endl;

	for(int i = 0; i < nc*mc; ++i){
		outfile << matC[i] << std::endl;
	}

	outfile.close();
   
    std::cout << "Calculation took " << time << " seconds\n";  
}
