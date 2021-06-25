#include <time.h>
#include "timer.h"


double timer(clock_t startTime, clock_t endTime){
  
  double diffTicks  =  (double)(startTime - endTime);
  double diffms     =  (diffTicks) / CLOCKS_PER_SEC;

  return diffms;
}
