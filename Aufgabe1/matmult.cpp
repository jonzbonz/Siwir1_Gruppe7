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
#define ALIGNMENT 32 //for AVX

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
	
	matA = new double[na*ma];
	//align(&matA);
	//out << matA << " after" << endl;

	//matA = NULL;
	//posix_memalign((void**)&matA, ALIGNMENT, na*ma*sizeof(double));

	//cout << "Geschafft!" << endl << "&matA = " << matA << endl;
	matB = new double[nb*mb];
	//matB = NULL;
	//posix_memalign((void**)&matB, ALIGNMENT, nb*mb*sizeof(double));


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
	//double* matC = NULL;
	//posix_memalign((void**)&matC, ALIGNMENT, nc*mc*sizeof(double));//
	double* matC = new double[nc*mc];
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

void strassenMult(double* matA, double* matB, double* matC, int nn){
	if(nn <= THRESHOLD){
		naiveMatmultTQ(matA, matB, matC, nn);
		return;
	}
	
	int nnh = nn/2;

	int nnhq = nnh*nnh; // = nn halb quadrat


	//double* M = NULL;
	//posix_memalign((void**)&M, ALIGNMENT, nnhq*7*sizeof(double));//
	
	double* M = new double[nnhq*7];

	//memset(M, 0, nnhq*7*sizeof(double));

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
	

	//double* matT = NULL;
	//posix_memalign((void**)&matT, ALIGNMENT, 7*2*nnhq);
	double* matT = new double[7*2*nnhq]; //nn = größe z.b. von matA, nnh = hälfte davon, nnhq = quadrat von nnh --> Größe von A11 z.b.
	
	//cout << "Matrix ist " << (7*2*nnhq*sizeof(double)/(1024*1024.0f)) << "MByte groß!" << endl;
	//memset(matT, 0, 7*2*nnhq*sizeof(double));


	if(nnh%4 != 0){
		cerr << "nnh ist nicht durch 4 teilbar!!! Threshold zu klein?" << endl;
		exit(-1);
	}
	
	/*cout << matA << endl;
	cout << matB << endl;
	cout << matC << endl;
	cout << matT << endl;
	cout << M << endl;*/

	//siwir::Timer timer;
	
	int pos = -1; //TODO -4 for AVX
	for(int i = 0; i < nnh; ++i){
		for(int j = 0; j < nnh; ++j){
			pos++;
			double A11 = matA[i*nn + j];
			double A12 = matA[i*nn + j + nnh];
			double A21 = matA[(i+nnh)*nn + j];
			double A22 = matA[(i+nnh)*nn + j + nnh];

			matT[0*nnhq + pos] = A11 + A22;
			matT[2*nnhq + pos] = A21 + A22;
			matT[4*nnhq + pos] = A11;
			matT[6*nnhq + pos] = A22;
			matT[8*nnhq + pos] = A11 + A12;
			matT[10*nnhq + pos] = A21 - A11;
			matT[12*nnhq + pos] = A12 - A22;
		}	
	}
	
	pos = -1;
	for(int i = 0; i < nnh; ++i){
		for(int j = 0; j < nnh; ++j){
			pos++;
			double B11 = matB[i*nn + j];
			double B21 = matB[i*nn + j + nnh];
			double B12 = matB[(i+nnh)*nn + j];
			double B22 = matB[(i+nnh)*nn + j + nnh];

			matT[1*nnhq + pos] = B11 + B22;
			matT[3*nnhq + pos] =  B11;
			matT[5*nnhq + pos] = B12 - B22;
			matT[7*nnhq + pos] = B21 - B11;
			matT[9*nnhq + pos] = B22;
			matT[11*nnhq + pos] = B11 + B12;
			matT[13*nnhq + pos] = B21 + B22;		
		}
	}
				

		/*for(int j = 0; j < nnh; j+=4){ //AVX!
			pos += 4;
			//TODO align to 32
			__m256d A11 = _mm256_loadu_pd(&matA[i*nn + j]);
			__m256d A12 = _mm256_loadu_pd(&matA[i*nn + j + nnh]);
			__m256d A21 = _mm256_loadu_pd(&matA[(i+nnh)*nn + j]);
			__m256d A22 = _mm256_loadu_pd(&matA[(i+nnh)*nn + j + nnh]);
			__m256d B11 = _mm256_loadu_pd(&matB[i*nn + j]);
			__m256d B12 = _mm256_loadu_pd(&matB[i*nn + j + nnh]);
			__m256d B21 = _mm256_loadu_pd(&matB[(i+nnh)*nn + j]);
			__m256d B22 = _mm256_loadu_pd(&matB[(i+nnh)*nn + j + nnh]);

			__m256d tmp0 = _mm256_add_pd(A11, A22);
			_mm256_storeu_pd(&matT[0*nnhq + pos], tmp0);
			__m256d tmp1 = _mm256_add_pd(B11, B22);
			_mm256_storeu_pd(&matT[1*nnhq + pos], tmp1);

			__m256d tmp2 = _mm256_add_pd(A21, A22);
			_mm256_storeu_pd(&matT[2*nnhq + pos], tmp2);
			_mm256_storeu_pd(&matT[3*nnhq + pos], B11);

			_mm256_storeu_pd(&matT[4*nnhq + pos], A11);
			__m256d tmp5 = _mm256_sub_pd(B12, B22);
			_mm256_storeu_pd(&matT[5*nnhq + pos], tmp5);

			_mm256_storeu_pd(&matT[6*nnhq + pos], A22);
			__m256d tmp7 = _mm256_sub_pd(B21, B11);
			_mm256_storeu_pd(&matT[7*nnhq + pos], tmp7);

			__m256d tmp8 = _mm256_add_pd(A11, A12);
			_mm256_storeu_pd(&matT[8*nnhq + pos], tmp8);
			_mm256_storeu_pd(&matT[9*nnhq + pos], B22);

			__m256d tmp10 = _mm256_sub_pd(A21, A11);
			_mm256_storeu_pd(&matT[10*nnhq + pos], tmp10);
			__m256d tmp11 = _mm256_add_pd(B11, B12);
			_mm256_storeu_pd(&matT[11*nnhq + pos], tmp11);

			__m256d tmp12 = _mm256_sub_pd(A12, A22);
			_mm256_storeu_pd(&matT[12*nnhq + pos], tmp12);
			__m256d tmp13 = _mm256_add_pd(B21, B22);
			_mm256_storeu_pd(&matT[13*nnhq + pos], tmp13);
		}
		}
		*/
	
   	//double time = timer.elapsed();
	//cout << "It took " << time << "seconds" << endl;


	//Not needed any more
	//delete[] matA;
	//delete[] matB;
	
	//This will be done recursive in the future
	/*
	naiveMatmultTQ(matT +  0*nnhq, matT +  1*nnhq, M + 0*nnhq, nnh);
	naiveMatmultTQ(matT +  2*nnhq, matT +  3*nnhq, M + 1*nnhq, nnh);
	naiveMatmultTQ(matT +  4*nnhq, matT +  5*nnhq, M + 2*nnhq, nnh);
	naiveMatmultTQ(matT +  6*nnhq, matT +  7*nnhq, M + 3*nnhq, nnh);
	naiveMatmultTQ(matT +  8*nnhq, matT +  9*nnhq, M + 4*nnhq, nnh);
	naiveMatmultTQ(matT + 10*nnhq, matT + 11*nnhq, M + 5*nnhq, nnh);
	naiveMatmultTQ(matT + 12*nnhq, matT + 13*nnhq, M + 6*nnhq, nnh);
	*/
	/*
	naiveMatmult(matT +  0*nnhq, matT +  1*nnhq, M + 0*nnhq, nnh, nnh, nnh);
	naiveMatmult(matT +  2*nnhq, matT +  3*nnhq, M + 1*nnhq, nnh, nnh, nnh);
	naiveMatmult(matT +  4*nnhq, matT +  5*nnhq, M + 2*nnhq, nnh, nnh, nnh);
	naiveMatmult(matT +  6*nnhq, matT +  7*nnhq, M + 3*nnhq, nnh, nnh, nnh);
	naiveMatmult(matT +  8*nnhq, matT +  9*nnhq, M + 4*nnhq, nnh, nnh, nnh);
	naiveMatmult(matT + 10*nnhq, matT + 11*nnhq, M + 5*nnhq, nnh, nnh, nnh);
	naiveMatmult(matT + 12*nnhq, matT + 13*nnhq, M + 6*nnhq, nnh, nnh, nnh);
	*/
	
	strassenMult(matT +  0*nnhq, matT +  1*nnhq, M + 0*nnhq, nnh);
	strassenMult(matT +  2*nnhq, matT +  3*nnhq, M + 1*nnhq, nnh);
	strassenMult(matT +  4*nnhq, matT +  5*nnhq, M + 2*nnhq, nnh);
	strassenMult(matT +  6*nnhq, matT +  7*nnhq, M + 3*nnhq, nnh);
	strassenMult(matT +  8*nnhq, matT +  9*nnhq, M + 4*nnhq, nnh);
	strassenMult(matT + 10*nnhq, matT + 11*nnhq, M + 5*nnhq, nnh);
	strassenMult(matT + 12*nnhq, matT + 13*nnhq, M + 6*nnhq, nnh);
	
	pos = -1;	
	for(int i = 0; i < nnh;  ++i){
		for(int j = 0; j < nnh; ++j){
			++pos;

			matC[i*nn + j] = M1[pos] + M4[pos] - M5[pos] + M7[pos];
			
			matC[i*nn + j + nnh] = M3[pos] + M5[pos];
			
			matC[(i+nnh)*nn + j] = M2[pos] + M4[pos];

			matC[(i+nnh)*nn + j + nnh] = M1[pos] - M2[pos] + M3[pos] + M6[pos];
		}
	}
}
