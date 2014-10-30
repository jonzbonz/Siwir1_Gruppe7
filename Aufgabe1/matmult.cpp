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

#define THRESHOLD 512
#define ALIGNMENT 32 //needed for AVX

inline void naiveMatmult(double* matA, double* matB, double* matC, int nc, int mc, int na);

inline void naiveMatmultTQ(double* matA, double* matTB, double* matC, int nn);

inline void naiveMatmultTQ_fast(double* matA, double* matTB, double* matC);

void strassenMult(double* matA, double* matB, double* matC, int nn);

//double matMultTime = 0; //TODO delete

int main(int argc, char* argv[])
{
	if(argc != 4){
		cerr << "Usage: ./matmult A.in B.in C.out" << endl;
		exit(1);
	}

//	cout << "with Threshold " << THRESHOLD << endl;

	//Reading input / creating matrices
	double* matA;// = NULL;
	double* matB;// = NULL;
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
	if(na == ma && nb == mb && nb >= THRESHOLD){
		int tmp = na;
		useStrassen = true;
		while(tmp > 1){ // Pruefen auf zweierpotenz
			if(tmp%2 == 1){
				useStrassen = false;
				break;
			}
			tmp = tmp >> 1;
		}	
	}
	
	if(useStrassen){
		posix_memalign((void**)&matA, ALIGNMENT, na*na*sizeof(double));
		posix_memalign((void**)&matB, ALIGNMENT, na*na*sizeof(double));
	} else {
		matA = new double[na*ma];
		matB = new double[nb*mb];
	}

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
	double* matC;// = NULL;
	if(useStrassen){
		posix_memalign((void**)&matC, ALIGNMENT, nc*mc*sizeof(double));
	} else {
		matC = new double[nc*mc];
	}

	memset(matC, 0, nc*mc*sizeof(double)); //TODO pruefen ob man das braucht

#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("matmult");
#endif

	

	siwir::Timer timer;


	if(!useStrassen){ 
	
		naiveMatmult(matA, matB, matC, nc, mc, na);

	} else {
		//Transpose B
		for(int i = 0; i < na-1; ++i){
			for(int j = i+1; j < na; ++j){
				double tmp = matB[i*na + j];
				matB[i*na + j] = matB[j*na + i];
				matB[j*na + i] = tmp;
			}
		}

		strassenMult(matA, matB, matC, na);
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
	
	free(matA);
	free(matB);
	free(matC);

    cout << THRESHOLD << "  " << time << endl;
//	cout << "It spent " <<  matMultTime << " seconds of that with naive matrix multiplication!" << endl;
//	cout << "That equals " << (float)(100*matMultTime/time) << " percent of the time!" << endl << endl;
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

inline void naiveMatmultTQ(double* matA, double* matBT, double* matC, int nn){
	//matA and matB have to be ALIGNMENT aligned!!!
	
//	siwir::Timer timer;

	double* temp;// = NULL;
	posix_memalign((void**)&temp, ALIGNMENT, 4*sizeof(double));

	for(int i = 0; i < nn; ++i){
		for(int j = 0; j < nn; ++j){
			__m256d sum = _mm256_setzero_pd();
			for(int k = 0; k < nn; k+=4){
				__m256d A = _mm256_load_pd(&matA [i*nn + k]);
				__m256d B = _mm256_load_pd(&matBT[j*nn + k]);
				__m256d C = _mm256_mul_pd(A, B);
				
				sum = _mm256_add_pd(sum, C);
			}

			_mm256_store_pd(temp, sum);
			matC[i*nn + j] = temp[0] + temp[1] + temp[2] + temp[3];
		}
	}

	free(temp);

//	matMultTime += timer.elapsed();
}

inline void naiveMatmultTQ_fast(double* matA, double* matBT, double* matC){
//	siwir::Timer timer;

	double* temp1;
	posix_memalign((void**)&temp1, ALIGNMENT, 8*4*sizeof(double));
	double* temp2 = temp1 + 4;
	double* temp3 = temp1 + 8;
	double* temp4 = temp1 + 12;
	double* temp5 = temp1 + 16;
	double* temp6 = temp1 + 20;
	double* temp7 = temp1 + 24;
	double* temp8 = temp1 + 28;

	
	for(int j = 0; j < THRESHOLD; j+=4){
		for(int i = 0; i < THRESHOLD; i+=2){
/*		
			__m256d sum1 = _mm256_setzero_pd();
			__m256d sum2 = _mm256_setzero_pd();
			__m256d sum3 = _mm256_setzero_pd();
			__m256d sum4 = _mm256_setzero_pd();
			__m256d sum5 = _mm256_setzero_pd();
			__m256d sum6 = _mm256_setzero_pd();
			__m256d sum7 = _mm256_setzero_pd();
			__m256d sum8 = _mm256_setzero_pd();
*/			
			
			__m256d B1 = _mm256_load_pd(&matBT[j*THRESHOLD]);	
			__m256d B2 = _mm256_load_pd(&matBT[(j+1)*THRESHOLD]);
			__m256d B3 = _mm256_load_pd(&matBT[(j+2)*THRESHOLD]);
			__m256d B4 = _mm256_load_pd(&matBT[(j+3)*THRESHOLD]);

			__m256d A1 = _mm256_load_pd(&matA[i*THRESHOLD]);
			__m256d A2 = _mm256_load_pd(&matA[(i+1)*THRESHOLD]);
			
			__m256d sum1 = _mm256_mul_pd(A1, B1);
			__m256d sum2 = _mm256_mul_pd(A2, B1);
			__m256d sum3 = _mm256_mul_pd(A1, B2);
			__m256d sum4 = _mm256_mul_pd(A2, B2);
			__m256d sum5 = _mm256_mul_pd(A1, B3);
			__m256d sum6 = _mm256_mul_pd(A2, B3);
			__m256d sum7 = _mm256_mul_pd(A1, B4);
			__m256d sum8 = _mm256_mul_pd(A2, B4);

			for(int k = 4; k < THRESHOLD; k+=4){
				B1 = _mm256_load_pd(&matBT[j*THRESHOLD + k]);	
				B2 = _mm256_load_pd(&matBT[(j+1)*THRESHOLD + k]);
				B3 = _mm256_load_pd(&matBT[(j+2)*THRESHOLD + k]);
				B4 = _mm256_load_pd(&matBT[(j+3)*THRESHOLD + k]);

				A1 = _mm256_load_pd(&matA[i*THRESHOLD + k]);
				A2 = _mm256_load_pd(&matA[(i+1)*THRESHOLD + k]);
				
				__m256d C1 = _mm256_mul_pd(A1, B1);
				__m256d C2 = _mm256_mul_pd(A2, B1);
				__m256d C3 = _mm256_mul_pd(A1, B2);
				__m256d C4 = _mm256_mul_pd(A2, B2);
				__m256d C5 = _mm256_mul_pd(A1, B3);
				__m256d C6 = _mm256_mul_pd(A2, B3);
				__m256d C7 = _mm256_mul_pd(A1, B4);
				__m256d C8 = _mm256_mul_pd(A2, B4);
			
				sum1 = _mm256_add_pd(sum1, C1);
				sum2 = _mm256_add_pd(sum2, C2);
				sum3 = _mm256_add_pd(sum3, C3);
				sum4 = _mm256_add_pd(sum4, C4);
				sum5 = _mm256_add_pd(sum5, C5);
				sum6 = _mm256_add_pd(sum6, C6);
				sum7 = _mm256_add_pd(sum7, C7);
				sum8 = _mm256_add_pd(sum8, C8);
			}

			_mm256_store_pd(temp1, sum1);
			matC[i*THRESHOLD + j] = temp1[0] + temp1[1] + temp1[2] + temp1[3];
			
			_mm256_store_pd(temp2, sum2);
			matC[(i+1)*THRESHOLD + j] = temp2[0] + temp2[1] + temp2[2] + temp2[3];
			
			_mm256_store_pd(temp3, sum3);
			matC[i*THRESHOLD + j + 1] = temp3[0] + temp3[1] + temp3[2] + temp3[3];
			
			_mm256_store_pd(temp4, sum4);
			matC[(i+1)*THRESHOLD + j + 1] = temp4[0] + temp4[1] + temp4[2] + temp4[3];

			_mm256_store_pd(temp5, sum5);
			matC[i*THRESHOLD + j + 2] = temp5[0] + temp5[1] + temp5[2] + temp5[3];
			
			_mm256_store_pd(temp6, sum6);
			matC[(i+1)*THRESHOLD + j + 2] = temp6[0] + temp6[1] + temp6[2] + temp6[3];
			
			_mm256_store_pd(temp7, sum7);
			matC[i*THRESHOLD + j + 3] = temp7[0] + temp7[1] + temp7[2] + temp7[3];
		
			_mm256_store_pd(temp8, sum8);
			matC[(i+1)*THRESHOLD + j + 3] = temp8[0] + temp8[1] + temp8[2] + temp8[3];
		}
	}

//	free(temp1);



//	matMultTime += timer.elapsed();	
}

void strassenMult(double* matA, double* matB, double* matC, int nn){
	if(nn == THRESHOLD){
		naiveMatmultTQ_fast(matA, matB, matC);
		return;
	} else if(nn < THRESHOLD){
		cerr << "Duerfte nicht passieren!" << endl;
		exit(1);
	}
	
	int nnh = nn/2;

	int nnhq = nnh*nnh; // = nn halb quadrat

	double* M;// = NULL;
	posix_memalign((void**)&M, ALIGNMENT, nnhq*7*sizeof(double));

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
	

	double* matT;// = NULL;
	posix_memalign((void**)&matT, ALIGNMENT, 7*2*nnhq*sizeof(double));

	double* matT0 = matT + nnhq*0;
	double* matT1 = matT + nnhq*2;
	double* matT2 = matT + nnhq*4;
	double* matT3 = matT + nnhq*6;
	double* matT4 = matT + nnhq*8;
	double* matT5 = matT + nnhq*10;
	double* matT6 = matT + nnhq*12;

	
/*
	double* matT0 = NULL;
	posix_memalign((void**)&matT0, ALIGNMENT, 2*nnhq*sizeof(double));
	double* matT1 = NULL;
	posix_memalign((void**)&matT1, ALIGNMENT, 2*nnhq*sizeof(double));
	double* matT2 = NULL;
	posix_memalign((void**)&matT2, ALIGNMENT, 2*nnhq*sizeof(double));
	double* matT3 = NULL;
	posix_memalign((void**)&matT3, ALIGNMENT, 2*nnhq*sizeof(double));
	double* matT4 = NULL;
	posix_memalign((void**)&matT4, ALIGNMENT, 2*nnhq*sizeof(double));
	double* matT5 = NULL;
	posix_memalign((void**)&matT5, ALIGNMENT, 2*nnhq*sizeof(double));
	double* matT6 = NULL;
	posix_memalign((void**)&matT6, ALIGNMENT, 2*nnhq*sizeof(double));
*/

	if(nnh%4 != 0){
		cerr << "nnh ist nicht durch 4 teilbar!!! Threshold zu klein?" << endl;
		exit(1);
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

	strassenMult(matT0, matT0 + nnhq, M + 0*nnhq, nnh);
	strassenMult(matT1, matT1 + nnhq, M + 1*nnhq, nnh);
	strassenMult(matT2, matT2 + nnhq, M + 2*nnhq, nnh);
	strassenMult(matT3, matT3 + nnhq, M + 3*nnhq, nnh);
	strassenMult(matT4, matT4 + nnhq, M + 4*nnhq, nnh);
	strassenMult(matT5, matT5 + nnhq, M + 5*nnhq, nnh);
	strassenMult(matT6, matT6 + nnhq, M + 6*nnhq, nnh);
	
	free(matT);

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
	
	free(M);
}
