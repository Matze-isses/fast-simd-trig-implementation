#ifndef CRYSTAL_CLOCK_HH
#define CRYSTAL_CLOCK_HH

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __linux__
  #ifndef _GNU_SOURCE
    #define _GNU_SOURCE
  #endif
#endif

#include <stdio.h>
#include <time.h>
#include <inttypes.h>
#include <stdbool.h>

#ifdef __linux__
  #include <bits/types/struct_timespec.h>
#endif

typedef uint64_t (*cur)();

struct CrystalClock {
  uint64_t _begin;
  uint64_t _end;
  double _freq;
};

uint64_t cycles1(struct CrystalClock C) ; 
double duration_ns1(struct CrystalClock C) ; 
double duration_us1(struct CrystalClock C) ;
double duration_ms1(struct CrystalClock C) ; 
double duration_s1(struct CrystalClock C)  ; 
uint64_t current();
uint64_t cycles( uint64_t aBegin,  uint64_t aEnd); 
double    duration_ns( uint64_t aNoCycles);
double    duration_us( uint64_t aNoCycles);
double    duration_ms( uint64_t aNoCycles);
double    duration_s( uint64_t aNoCycles);
double    duration_min( uint64_t aNoCycles);
double    duration_h( uint64_t aNoCycles);
uint64_t  duration_ns_to_cc( double aDuration);
uint64_t  duration_us_to_cc( double aDuration);
uint64_t  duration_ms_to_cc( double aDuration);
uint64_t  duration_s_to_cc( double aDuration);
uint64_t  duration_min_to_cc( double aDuration);
uint64_t  duration_h_to_cc( double aDuration);
void init(); 
double frequency(); 
double cycle_time_ns();
double measure_frequency();

#ifdef __cplusplus
}
#endif

#endif
