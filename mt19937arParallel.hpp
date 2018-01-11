#ifndef MTPARALLEL
#define MTPARALLEL

#include <stdio.h>
#include <cmath>
#define N 624

typedef struct stateMT_st{
  unsigned long mt[N];
  int mti; 
} StateRN;

/* initializes mt[N] with a seed */
void init_sgenrand(unsigned long s, StateRN *state);

/* generates a random number on [0,0xffffffff]-interval */
unsigned long sgenrand_int32(StateRN *state);

/* generates a random number on [0,1]-real-interval */
double sgenrand_real1(StateRN *state);

/* generates a random number on [0,1)-real-interval */
double sgenrand_real2(StateRN *state);

/* generates a random number on (0,1)-real-interval */
double sgenrand_real3(StateRN *state);

double sgaussrand(StateRN *state);

#endif
