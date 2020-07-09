#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "common_helpers.h"
#include "mtwister.h"
#include <sys/time.h>
#include <string.h>

int main (int argc, char **argv) {
    struct timeval tv;
    unsigned long t;
    MTRand myrng;
    double a;
    int n;
    double limit;
    FILE *fptr;
    
    limit = atof(argv[1]);
    
    printf("Limit is: %f\n", limit);
    
    /* Seed the random number generator with the current time */
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;
    
    /* Set up the MTwister random number generator */
    myrng = seedRand(t);
    
    printf("How many random numbers to generate?\n");
    scanf("%d", &n);
    
    fptr = fopen("gaussian_sampling.dat", "w");
    
    for (int i = 0; i < n; i++) {
        a = gaussian_random_tail(0.0, 1.0, limit, &myrng);
        //double Z[2];
        //gaussian_random(0.0, 1.0, Z, &myrng);
        //a = Z[0];
        fprintf(fptr,"%1.4f\n", a); 
        //assert(a > limit);
    }
    
    fclose(fptr);
    
    return 0;
}
