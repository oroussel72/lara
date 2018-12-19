#include <time.h>
#include <omp.h>
#include <sys/times.h>
#include "Timer.h"

/*
______________________________________________________________________________________________

Constructor
______________________________________________________________________________________________

*/

Timer::Timer()
{
#ifdef _OPENMP
  _startOMPTime = omp_get_wtime();
#else
  _startCPUTime = clock();
#endif
  
  _sumCPUtime = 0.;
  _isTimerOn = true;

}

/*
______________________________________________________________________________________________

Distructor
______________________________________________________________________________________________

*/

Timer::~Timer()
{
}

/*
______________________________________________________________________________________________

Start CPU time
______________________________________________________________________________________________

*/

void Timer::start()
{
#ifdef _OPENMP
  _startOMPTime = omp_get_wtime();
#else
  _startCPUTime = clock();
#endif
  
  _isTimerOn = true;
  return;
}

/*
______________________________________________________________________________________________

Reset CPU time
______________________________________________________________________________________________

*/

void Timer::reset()
{
#ifdef _OPENMP
  _startOMPTime = omp_get_wtime();
#else
  _startCPUTime = clock();
#endif
  
  _isTimerOn = true;
  _sumCPUtime = 0.;
  return;
}
/*
______________________________________________________________________________________________

Stop CPU time
______________________________________________________________________________________________

*/

double Timer::stop()
{
  clock_t endCPUTime;
  double endOMPTime;
	
  // --- Execution ---

#ifdef _OPENMP
  endOMPTime = omp_get_wtime();
  _sumCPUtime += endOMPTime - _startOMPTime;
#else
  endCPUTime = clock();
  _sumCPUtime  += double((unsigned long int)endCPUTime-_startCPUTime)/ (unsigned long int)CLOCKS_PER_SEC;
#endif
  
  _isTimerOn = false;
  return _sumCPUtime;
}

/*
______________________________________________________________________________________________

Check CPU time
______________________________________________________________________________________________

*/
void Timer::check()
{
  clock_t endCPUTime;
  double endOMPTime;

  // --- Execution ---

#ifdef _OPENMP
  endOMPTime = omp_get_wtime();
  _sumCPUtime += endOMPTime - _startOMPTime;
  _startOMPTime = endOMPTime;
#else
  endCPUTime = clock();
  _sumCPUtime  += double((unsigned long int)endCPUTime-_startCPUTime)/ (unsigned long int)CLOCKS_PER_SEC;
  _startCPUTime =  endCPUTime;
#endif
}

/*
______________________________________________________________________________________________

	Return CPU time from previous start in seconds
______________________________________________________________________________________________

*/
double Timer::CPUTime()
{
  if (_isTimerOn)
    check();

  return _sumCPUtime;
}

