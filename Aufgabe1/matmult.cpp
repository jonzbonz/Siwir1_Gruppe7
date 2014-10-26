/*==========================================================================
 * TODO Enter todo's														|
 * ueberpruefen, ob man loop umbauen kann --> nur eine Laufvariable			|
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

#include "immintrin.h"

using namespace std;
 
#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}
#endif

#define THRESHOLD 128

inline void naiveMatmult(double* matA, double* matB, double* matC, int nc, int mc, int na);

inline void naiveMatmultQ(double* matA, double* matB, double* matC, int nn);

inline void naiveMatmultTQ(double* matA, double* matTB, double* matC, int nn);

void strassenMult(double* matA, double* matB, double* matC, int nn);

//void align(void** pointer){
//	
//}

int main(int argc, char* argv[])
{
	if(argc != 4){
		cerr << "Usage: ./matmult A.in B.in C.out" << endl;
		return 1;
	}
	cout << "with Threshold " << THRESHOLD << endl;
	//Reading input / creating matrices
	double* matA = NULL;
	double* matB = NULL;
	int na, ma, nb, mb;

	string line;
	
	ifstream firstMatFile(argv[1]);
	getline(firstMatFile, line);
	ma = atoi(line.c_str());
	line.erase(0, line.find(' ') + 1);
	na = atoi(line.c_str());

	ifstream secondMatFile(argv[2]);
	getline(secondMatFile, line);
	mb = atoi(line.c_str());
	line.erase(0, line.find(' ') + 1);
	nb = atoi(line.c_str());
	
	bool useStrassen = false;
	int nn = 1;
	if(na == ma && nb == mb && nb*ma >= 0){ //TODO gute grenze finden
		useStrassen = true;
		
		nn=na; //WORKAROUND!
		/*	
		while(nn < na) nn = nn << 1;
		//TODO THIS IS BUGGED! (Ask joni)
		
		matA = new double[nn*nn];
		memset(matA, 0, nn*nn*sizeof(double));
		matB = new double[nn*nn];
		memset(matB, 0, nn*nn*sizeof(double));
		*/
	} else {
		matA = new double[na*ma];
		matB = new double[nb*mb];
	}
	
	//matA = new double[na*ma];
	matA = NULL;
	posix_memalign((void**)&matA, 32, na*ma*sizeof(double));
	//matB = new double[nb*mb];
	matB = NULL;
	posix_memalign((void**)&matB, 32, nb*mb*sizeof(double));

	for(int i = 0; i < ma; ++i){
		for(int j = 0; j < na; ++j){
			if(!getline(firstMatFile, line)){
				cerr << "Unexpected end of file!" << endl;
				return 1;
			}
			matA[i*na + j] = atof(line.c_str());
		}
	}
	firstMatFile.close();

	for(int i = 0; i < mb; ++i){
		for(int j = 0; j < nb; ++j){
			if(!getline(secondMatFile, line)){
				cerr << "Unexpected end of file!" << endl;
				return 1;
			}
			matB[i*nb + j] = atof(line.c_str());
		}
	}
	secondMatFile.close();

 

	int nc = nb;
	int mc = ma;
	//double* matC = new double[nc*mc];
	double* matC = NULL;
	posix_memalign((void**)&matC, 32, nc*mc*sizeof(double));

	memset(matC, 0, nc*mc*sizeof(double));

#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("matmult");
#endif

	

	siwir::Timer timer;


	if(!useStrassen){ 
	
		naiveMatmult(matA, matB, matC, nc, mc, na);

	} else {
		//cout << "using strassen" << endl;	
		//Strassen-Algorithmus - nur fuer grosse und quadratische matrizen

		//Transpose B
		for(int i = 0; i < nn-1; ++i){
			for(int j = i+1; j < nn; ++j){
				double tmp = matB[i*nn + j];
				matB[i*nn + j] = matB[j*nn + i];
				matB[j*nn + i] = tmp;
			}
		}

		strassenMult(matA, matB, matC, nn);
	}

   	double time = timer.elapsed();

#ifdef USE_LIKWID
	likwid_markerStopRegion("matmult");
	likwid_markerClose();
#endif  

	
	ofstream outfile(argv[3]);
	outfile << mc << " " << nc << endl;

	for(int i = 0; i < nc*mc; ++i){ 
		outfile << matC[i] << endl;
	}

	matC = NULL;
	posix_memalign((void**)&matC, 32, nc*mc*sizeof(double));
	outfile.close();
   
	delete[] matC;

    cout << "Calculation took " << time << " seconds" << endl;  
}

inline void naiveMatmult(double* matA, double* matB, double* matC, int nc, int mc, int na){
	for(int j = 0; j < nc; ++j){							// ueber die Spalten von matC
		for(int i = 0; i  < mc; ++i){						// ueber die Zeilen von matC
			for(int k = 0;  k  < na/*bzw mb*/; ++k){		// ueber die Spalten von matA / Zeilen von matB
				matC[i*nc + j] += matA[i*na + k] * matB[k*nc + j];
			}
		} 
	}
}

inline void naiveMatmultQ(double* matA, double* matB, double* matC, int nn){
	
	//transpose matB
	double* matBT = new double[nn*nn];
	for(int i = 0; i < nn; ++i){
		for(int j = 0; j < nn; ++j){
			matBT[i*nn + j] = matB[j*nn +i];
		}
	}

/*
	for(int j = 0; j < nn; ++j){							// ueber die Spalten von matC
		for(int i = 0; i  < nn; ++i){						// ueber die Zeilen von matC
			for(int k = 0;  k  < nn; ++k){					// ueber die Spalten von matA / Zeilen von matB
				matC[i*nn + j] += matA[i*nn + k] * matBT[j*nn + k];
			}
		} 
	}
*/

	double* temp = new double[4];
	for(int j = 0; j < nn; ++j){
		for(int i = 0; i < nn; ++i){							// ueber die Zeilen und spalten von matC
			__m256d sum = _mm256_setzero_pd();
			for(int k = 0; k < nn; k+=4){						// ueber die Spalten von matA / Zeilen von matB
				__m256d A = _mm256_loadu_pd(&matA [i*nn + k]);
				__m256d B = _mm256_loadu_pd(&matBT[j*nn + k]);
				__m256d C = _mm256_mul_pd(A, B);
				
				sum = _mm256_add_pd(sum, C);
			}
			_mm256_storeu_pd(temp, sum);
			matC[i*nn + j] = temp[0] + temp[1] + temp[2] + temp[3];
		}
	}
}

inline void naiveMatmultTQ(double* matA, double* matBT, double* matC, int nn){
	//matA and matB have to be 32 aligned!!!
	//double* temp = new double[4];
	double* temp = NULL;
	posix_memalign((void**)&temp, 32, 4*sizeof(double));

	for(int j = 0; j < nn; ++j){
		for(int i = 0; i < nn; ++i){							// ueber die Zeilen und spalten von matC
			__m256d sum = _mm256_setzero_pd();
			for(int k = 0; k < nn; k+=4){						// ueber die Spalten von matA / Zeilen von matB
				__m256d A = _mm256_load_pd(&matA [i*nn + k]);
				__m256d B = _mm256_load_pd(&matBT[j*nn + k]);
				__m256d C = _mm256_mul_pd(A, B);
				
				sum = _mm256_add_pd(sum, C);
			}
			_mm256_store_pd(temp, sum);
			matC[i*nn + j] = temp[0] + temp[1] + temp[2] + temp[3];
		}
	}
}

void strassenMult(double* matA, double* matB, double* matC, int nn){
	if(nn <= THRESHOLD){
		naiveMatmultTQ(matA, matB, matC, nn);
		return;
	}
	
	int nnh = nn/2;

	int nnhq = nnh*nnh; // = nn halb quadrat


	
	//double* M = new double[nnhq*7];
	double* M = NULL;
	posix_memalign((void**)&M, 32, nnhq*7*sizeof(double));

	double* M1 = M + nnhq*0;
	double* M2 = M + nnhq*1;
	double* M3 = M + nnhq*2;
	double* M4 = M + nnhq*3;
	double* M5 = M + nnhq*4;
	double* M6 = M + nnhq*5;
	double* M7 = M + nnhq*6;
		
	// Fuer den naechsen schritt brauchen wir folgende Matritzen:
	// A11+A22, A21+A22, A11,     A22,     A11+A12, A21-A11, A12-A22
	// B11+B22, B11,     B12-B22, B21-B11, B22,     B11+B12, B21+B22
	//    0   ,    2   ,    4   ,    6   ,    8   ,    10  ,    12 
	//    1   ,    3   ,    5   ,    7   ,    9   ,    11  ,    13
	// Dafuer sollten untereinanderstehende Matritzen auch nacheinander im Cache stehen (Meint Joni)
	// Diese werden naemlich anschließend miteinander multipliziert, um M zu erhalten
	// Die entstehende Matrix nennen wir matT (für temp)
	

	//double* matT = new double[7*2*nnhq]; //nn = größe z.b. von matA, nnh = hälfte davon, nnhq = quadrat von nnh --> Größe von A11 z.b.
	double* matT = NULL;
	posix_memalign((void**)&matT, 32, 7*2*nnhq*sizeof(double));
	
	double* matT0 = matT + nnhq*0;
	double* matT1 = matT + nnhq*2;
	double* matT2 = matT + nnhq*4;
	double* matT3 = matT + nnhq*6;
	double* matT4 = matT + nnhq*8;
	double* matT5 = matT + nnhq*10;
	double* matT6 = matT + nnhq*12;

	
	//cout << "Matrix ist " << (7*2*nnhq*sizeof(double)/(1024*1024.0f)) << "MByte groß!" << endl;
	//memset(matT, 0, 7*2*nnhq*sizeof(double));


	if(nnh%4 != 0){
		cerr << "nnh ist nicht durch 4 teilbar!!! Threshold zu klein?" << endl;
		exit(-1);
	}
	
	int pos = -4;
	for(int i = 0; i < nnh; ++i){
		for(int j = 0; j < nnh; j+=4){
			pos+=4;

			__m256d A11 = _mm256_load_pd(&matA[i*nn + j]);
			__m256d A12 = _mm256_load_pd(&matA[i*nn + j + nnh]);
			__m256d A21 = _mm256_load_pd(&matA[(i+nnh)*nn + j]);
			__m256d A22 = _mm256_load_pd(&matA[(i+nnh)*nn + j + nnh]);

			_mm256_store_pd(&matT0[pos], _mm256_add_pd(A11, A22));
			_mm256_store_pd(&matT1[pos], _mm256_add_pd(A21, A22));
			_mm256_store_pd(&matT2[pos], A11);
			_mm256_store_pd(&matT3[pos], A22);
			_mm256_store_pd(&matT4[pos], _mm256_add_pd(A11, A12));
			_mm256_store_pd(&matT5[pos], _mm256_sub_pd(A21, A11));
			_mm256_store_pd(&matT6[pos], _mm256_sub_pd(A12, A22));
		}	
	}
	
	pos = -4 + nnhq;
	for(int i = 0; i < nnh; ++i){
		for(int j = 0; j < nnh; j+=4){
			pos+=4;
			__m256d B11 = _mm256_load_pd(&matB[i*nn + j]);
			__m256d B21 = _mm256_load_pd(&matB[i*nn + j + nnh]);
			__m256d B12 = _mm256_load_pd(&matB[(i+nnh)*nn + j]);
			__m256d B22 = _mm256_load_pd(&matB[(i+nnh)*nn + j + nnh]);
			
			_mm256_store_pd(&matT0[pos], _mm256_add_pd(B11, B22));
			_mm256_store_pd(&matT1[pos], B11);
			_mm256_store_pd(&matT2[pos], _mm256_sub_pd(B12, B22));
			_mm256_store_pd(&matT3[pos], _mm256_sub_pd(B21, B11));
			_mm256_store_pd(&matT4[pos], B22);
			_mm256_store_pd(&matT5[pos], _mm256_add_pd(B11, B12));
			_mm256_store_pd(&matT6[pos], _mm256_add_pd(B21, B22));		
		}
	}


	//Not needed any more
	//delete[] matA;
	//delete[] matB;
	
	strassenMult(matT0, matT0 + nnhq, M + 0*nnhq, nnh);
	strassenMult(matT1, matT1 + nnhq, M + 1*nnhq, nnh);
	strassenMult(matT2, matT2 + nnhq, M + 2*nnhq, nnh);
	strassenMult(matT3, matT3 + nnhq, M + 3*nnhq, nnh);
	strassenMult(matT4, matT4 + nnhq, M + 4*nnhq, nnh);
	strassenMult(matT5, matT5 + nnhq, M + 5*nnhq, nnh);
	strassenMult(matT6, matT6 + nnhq, M + 6*nnhq, nnh);
	

	double* tmp1 = NULL;
	posix_memalign((void**)&tmp1, 32, 4*sizeof(double));
	double* tmp2 = NULL;
	posix_memalign((void**)&tmp2, 32, 4*sizeof(double));
	double* tmp3 = NULL;
	posix_memalign((void**)&tmp3, 32, 4*sizeof(double));
	double* tmp4 = NULL;
	posix_memalign((void**)&tmp4, 32, 4*sizeof(double));
	
	pos = -4;
	for(int i = 0; i < nnh;  ++i){
		for(int j = 0; j < nnh; j+=4){
			pos+=4;
			
			__m256d m1 = _mm256_load_pd(&M1[pos]);
			__m256d m2 = _mm256_load_pd(&M2[pos]);
			__m256d m3 = _mm256_load_pd(&M3[pos]);
			__m256d m4 = _mm256_load_pd(&M4[pos]);
			__m256d m5 = _mm256_load_pd(&M5[pos]);
			__m256d m6 = _mm256_load_pd(&M6[pos]);
			__m256d m7 = _mm256_load_pd(&M7[pos]);

			_mm256_store_pd(&matC[i*nn + j], _mm256_add_pd(_mm256_add_pd(m1, m4), _mm256_sub_pd(m7, m5)));

			_mm256_store_pd(&matC[i*nn + j + nnh], _mm256_add_pd(m3, m5));

			_mm256_store_pd(&matC[(i+nnh)*nn + j], _mm256_add_pd(m2, m4));
				
			_mm256_store_pd(&matC[(i+nnh)*nn + j + nnh],_mm256_add_pd(_mm256_sub_pd(m1, m2), _mm256_add_pd(m3, m6)));

		}
	}
}
