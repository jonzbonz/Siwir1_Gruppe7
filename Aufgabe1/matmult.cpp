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
	if(na == ma && nb == mb && nb*ma >= 0){ //TODO gute grenze finden
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
		
		int nnh = nn/2;

		int offset = nnh*nnh;
		double* M = new double[offset*7];
		memset(M, 0, offset*7*sizeof(double));

		double* M1 = M;
		double* M2 = M1 + offset;
		double* M3 = M2 + offset;
		double* M4 = M3 + offset;
		double* M5 = M4 + offset;
		double* M6 = M5 + offset;
		double* M7 = M6 + offset;

		for(int j = 0; j < nnh; ++j){
			for(int i = 0; i < nnh; ++i){
				
				int pos = i*nnh + j;

				for(int k = 0; k < nnh; ++k){
					//(A11 + A22) * (B11 + B22)
					M1[pos] += (matA[i*nn + k] + matA[(i+nnh)*nn + k + nnh]) * (matB[k*nn + j] + matB[(k+nnh)*nn + j + nnh]);
					//(A21 + A22) * B11
					M2[pos] += (matA[(i+nnh)*nn + k] + matA[(i+nnh)*nn + k + nnh]) * matB[k*nn + j];
					//A11 * (B12 - B22)
					M3[pos] += matA[i*nn + k] * (matB[k*nn + j + nnh] - matB[(k+nnh)*nn + j + nnh]);
					//A22 * (B21 - B11)
					M4[pos] += matA[(i+nnh)*nn + k + nnh] * (matB[(k+nnh)*nn + j] - matB[k*nn + j]);
					//(A11 + A12) * B22
					M5[pos] += (matA[i*nn + k] + matA[i*nn + k + nnh]) *  matB[(k+nnh)*nn + j + nnh];
					//(A21 - A11) * (B11 + B12)
					M6[pos] += (matA[(i+nnh)*nn + k] - matA[i*nn + k]) * (matB[k*nn + j] + matB[k*nn + j + nnh]);
					//(A12 - A22) * (B21 + B22)
					M7[pos] += (matA[i*nn + k + nnh] - matA[(i+nnh)*nn + k + nnh]) * (matB[(k+nnh)*nn + j] + matB[(k+nnh)*nn + j + nnh]);
				}

				matC[i*nn + j] = M1[pos] + M4[pos] - M5[pos] + M7[pos];
				
				matC[i*nn + j + nnh] = M3[pos] + M5[pos];
				
				matC[(i+nnh)*nn + j] = M2[pos] + M3[pos];

				matC[(i+nnh)*nn + j + nnh] = M1[pos] - M2[pos] + M3[pos] + M6[pos];
			}
		}
		
		
		for(int i = 1; i < 8; ++i){
			std::cout << "M" << i << ":" << std::endl;
			for(int j = 0; j < nnh; ++j){
				for(int k = 0; k < nnh; ++k){
					std::cout << M[i*nnh*nnh + j*nnh + k] << " ";
				}
				std::cout << std::endl;
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
