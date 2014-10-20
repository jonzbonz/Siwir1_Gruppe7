#include "Timer.h"

#include <iostream>
#include <vector>

extern "C" {
#include <cblas.h>
}

#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}
#endif



int main()
{
#ifdef USE_LIKWID
   likwid_markerInit();
   likwid_markerStartRegion( "matmult" );
#endif
   
   siwir::Timer timer;
  
	// bitte hier koten

   double time = timer.elapsed();
   
#ifdef USE_LIKWID
   likwid_markerStopRegion( "matmult" );
#endif   
   
    std::cout << "Calculation took " << time << " seconds\n";  
}
