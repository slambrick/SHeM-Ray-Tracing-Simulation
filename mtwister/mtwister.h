#ifndef __MTWISTER_H
#define __MTWISTER_H

#define STATE_VECTOR_LENGTH 624
#define STATE_VECTOR_M      397 /* changes to STATE_VECTOR_LENGTH also require changes to this */

typedef struct tagMTRand {
    unsigned long mt[STATE_VECTOR_LENGTH];
    int index;
} MTRand;

void seedRand(unsigned long seed, MTRand* rand);
void genRandLong(MTRand* rand, unsigned long* y);
void genRand(MTRand* rand, double* randNum);

#endif /* #ifndef __MTWISTER_H */
