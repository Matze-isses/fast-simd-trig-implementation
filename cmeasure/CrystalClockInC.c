#include "./CrystalClockInC.h"
#include <limits.h>
#include <stdbool.h>

#include <stdio.h>
#include <x86intrin.h>

#include "./cmeasure.h"

double _freq = -1;

uint64_t
current() {
  uint64_t lCyc = 0;
  #ifdef __x86_64__
    lCyc =  __rdtsc();
  #elif (__ARM_ARCH == 8)
    asm volatile("mrs %0, cntvct_el0" : "=r" (lCyc));
  #else 
    #error CrystalClock: unsupported architecture
  #endif
  return lCyc;
}

//struct CrystalClock C ={_begin(0), _end(0)} ;

uint64_t
cycles1(struct CrystalClock C){
  return cycles(C._begin, C._end);
}

uint64_t
cycles(const uint64_t aBegin, const uint64_t aEnd) {
  uint64_t lRes = 0;
  if(aBegin > aEnd) {
    lRes = (ULLONG_MAX - aBegin) + aEnd;
  } else {
    lRes = aEnd - aBegin;
  }
  return lRes;
}

double
frequency() {
  if(0 > _freq) {
    _freq = measure_frequency();
  }
  return _freq;
}

double
cycle_time_ns() {
  return ((1 / frequency()) * ((double) 1000 * 1000 * 1000));
}

void
init() {
  _freq = measure_frequency();
}

double
measure_frequency() {
  const bool lTrace = false;
  struct cmeasure_t lMeas;
  cmeasure_start(&lMeas);
  const uint64_t lBegin = current();
  do {
    cmeasure_stop(&lMeas);
  } while(1.0 > cmeasure_total_s(&lMeas));
  const uint64_t lEnd = current();
  const uint64_t lNoCyc = cycles(lBegin, lEnd);
  const double lFreq = (((double) (lNoCyc)) / cmeasure_total_s(&lMeas));
  //if(lTrace) {
  //  std::cout << "measure_frequency:" << std::endl;
  //  std::cout << "   time " << cmeasure_total_s(&lMeas) << " [s]" << std::endl;
  //  if(lEnd < lBegin) {
  //    std::cout << "   wrap arround." << std::endl;
  //  }
  //  std::cout << "   cycles    :" << lNoCyc << std::endl;
  //  std::cout << "   frequency :" << lFreq   << std::endl;
  //}
  return lFreq;
}

double
duration_ns1(struct CrystalClock C) {
  return duration_s1(C) * ((double) 1000 * 1000 * 1000);
}

double
duration_us1(struct CrystalClock C) {
  return duration_s1(C) * ((double) 1000 * 1000);
}

double
duration_ms1(struct CrystalClock C) {
  return duration_s1(C) * 1000;
}

double
duration_s1(struct CrystalClock C) {
  return (((double) cycles1(C)) / frequency());
}

double
duration_ns(const uint64_t aNoCycles) {
    return duration_s(aNoCycles) * ((double) 1000 * 1000 * 1000);

}

double
duration_us(const uint64_t aNoCycles) {
    return duration_s(aNoCycles) * ((double) 1000 * 1000);
}

double
duration_ms(const uint64_t aNoCycles) {
    return duration_s(aNoCycles) * ((double) 1000);
}

double
duration_s(const uint64_t aNoCycles) {
    return (((double) aNoCycles) / frequency());
}

double
duration_min(const uint64_t aNoCycles) {
  return (duration_s(aNoCycles) / 60);
}

double
duration_h(const uint64_t aNoCycles) {
  return (duration_s(aNoCycles) / 3600);
}

uint64_t
duration_ns_to_cc(const double aDuration) {
  return (frequency() * aDuration / ((double) 1000000000));
}

uint64_t
duration_us_to_cc(const double aDuration) {
  return (frequency() * aDuration / ((double) 1000000));
}

uint64_t
duration_ms_to_cc(const double aDuration) {
  return (frequency() * aDuration / ((double) 1000));
}

uint64_t
duration_s_to_cc(const double aDuration) {
  return (frequency() * aDuration);
}

uint64_t
duration_min_to_cc(const double aDuration) {
  return (frequency() * aDuration * 60);
}

uint64_t
duration_h_to_cc(const double aDuration) {
  return (frequency() * aDuration * 3600);
}


