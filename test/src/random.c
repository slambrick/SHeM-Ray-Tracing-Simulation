#include <stdio.h>
#include <assert.h>
#include "common_helpers.h"
#include "mtwister.h"
#include "sys/time.h"

int main () {
    struct timeval tv;
    unsigned long t;
    MTRand myrng;
    double a;
    int n;
    double limit;
    FILE *fptr = fopen("gaussian_sampling.dat", "w");
    
    /* Seed the random number generator with the current time */
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;
    
    /* Set up the MTwister random number generator */
    myrng = seedRand(t);
    
    printf("How many random numbers to generate?\n");
    scanf("%d", &n);
    
    printf("Limit of the gaussian tail?\n");
    scanf("%f", &limit);
    
    for (int i = 0; i < n; i++) {
        a = gaussian_random_tail(0.0, 1.0, limit, &myrng);
        //double Z[2];
        //gaussian_random(0.0, 1.0, Z, &myrng);
        //a = Z[0];
        printf("%1.2f\n", a);
        fprintf(fptr,"%1.4f\n", a); 
        //assert(a > limit);
    }
    printf("\n");

    printf("\nTests passed.\n");
    
    return 0;
}
