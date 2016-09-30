#ifndef MT
#include "MT.h"

#endif

#include <stdint.h>



/* initializes mt[MT_N] with a seed */
void init_genrand(uint32_t s)
{
	mt[0] = s & 0xffffffffUL;
	for (mti = 1; mti<MT_N; mti++) {
		mt[mti] =
			(1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array mt[].                        */
		/* 2002/01/09 modified by Makoto Matsumoto             */
		mt[mti] &= 0xffffffffUL;
		/* for >32 bit machines */
	}
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(uint32_t init_key[], int32_t key_length)
{
	int32_t i, j, k;
	init_genrand(19650218UL);
	i = 1; j = 0;
	k = (MT_N>key_length ? MT_N : key_length);
	for (; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL))
			+ init_key[j] + j; /* non linear */
		mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
		i++; j++;
		if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i = 1; }
		if (j >= key_length) j = 0;
	}
	for (k = MT_N - 1; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL))
			- i; /* non linear */
		mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
		i++;
		if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i = 1; }
	}

	mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
	uint32_t y;
	static uint32_t mag01[2] = { 0x0UL, MATRIX_A };
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if (mti >= MT_N) { /* generate N words at one time */
		int32_t kk;

		if (mti == MT_N + 1)   /* if init_genrand() has not been called, */
			init_genrand(5489UL); /* a default initial seed is used */

		for (kk = 0; kk<MT_N - MT_M; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (; kk<MT_N - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (MT_M - MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[MT_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}

	y = mt[mti++];

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
	return (long)(genrand_int32() >> 1);
}

static uint32_t next(int32_t bits)
{
	return genrand_int32() >> (32 - bits);    // hope that's right!
}

bool nextBoolean()
{
	return next(1) != 0;
}

double nextDouble()
{
	return (((uint64_t)next(27) << 26) + next(26))
		/ (double)((uint64_t)1 << 53);
}

/* generates a random number on [0,1) with 53-bit resolution*/
double nextDoubleIE()
{
	return nextDouble();
}

/* generates a random number on (0,1] with 53-bit resolution*/
double nextDoubleEI()
{
	return (((uint64_t)next(27) << 26) + next(26) + 1)
		/ (double)((uint64_t)1 << 53);
}

	/* generates a random number on (0,1) with 53-bit resolution*/
double nextDoubleEE()
{
	double retval;
	do {
		retval = nextDoubleIE();
	} while (retval == (double)0);
	return retval;
}

/* generates a random number on [0,1] with 53-bit resolution*/
double nextDoubleII()
{
	uint64_t num;
	const int64_t Q = 2047;
	const int64_t SHIFT32P1 = ((uint64_t)1 << 53) + (uint64_t)1;

	do{
		num = ((uint64_t)next(32) << 32) + ((uint64_t)next(32));
	} while (num >= ((SHIFT32P1)*Q));
	num %= SHIFT32P1;
	return ((double)num / ((uint64_t)1 << 53));
}



