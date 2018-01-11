#include <stdio.h>
#include <math.h>
#include "mt19937arParallel.hpp"

/* Period parameters */  

#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

/* initializes mt[N] with a seed */
void init_sgenrand(unsigned long s, StateRN *state) {
  state->mti=N;
  state->mt[0]= s & 0xffffffffUL;
  for (state->mti=1; state->mti<N; state->mti++) {
    state->mt[state->mti] = 
      (1812433253UL * (state->mt[state->mti-1] ^ (state->mt[state->mti-1] >> 30)) + state->mti); 
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, MSBs of the seed affect   */
    /* only MSBs of the array mt[].                        */
    /* 2002/01/09 modified by Makoto Matsumoto             */
    state->mt[state->mti] &= 0xffffffffUL;
    /* for >32 bit machines */
  }
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long sgenrand_int32(StateRN *state)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (state->mti >= N) { /* generate N words at one time */
        int kk;

        if (state->mti == N+1)   /* if init_genrand() has not been called, */
	  init_sgenrand(5489UL, state); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (state->mt[kk]&UPPER_MASK)|(state->mt[kk+1]&LOWER_MASK);
            state->mt[kk] = state->mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (state->mt[kk]&UPPER_MASK)|(state->mt[kk+1]&LOWER_MASK);
            state->mt[kk] = state->mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (state->mt[N-1]&UPPER_MASK)|(state->mt[0]&LOWER_MASK);
        state->mt[N-1] = state->mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        state->mti = 0;
    }
  
    y = state->mt[state->mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,1]-real-interval */
double sgenrand_real1(StateRN *state)
{
    return sgenrand_int32(state)*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double sgenrand_real2(StateRN *state)
{
    return sgenrand_int32(state)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double sgenrand_real3(StateRN *state)
{
    return (((double)sgenrand_int32(state)) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a Gaussian number with mu = 0 and sigma 1 (it use sgenrand_real1) */
double sgaussrand(StateRN *state) {
  double v,u;
  u=sqrt(-2.0*log(sgenrand_int32(state)*(1.0/4294967295.0)));
  v=2.0*3.14159265359*(sgenrand_int32(state)*(1.0/4294967295.0));
  return u*cos(v);
}
