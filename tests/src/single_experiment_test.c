/*
 * single_ray.c
 *
 *  Created on: 26 Oct 2020
 *      Author: Sam Lambrick
 */

#include "atom_ray_tracing3D.h"
#include "mtwister.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdint-gcc.h>

int main(int argc, char * argv []) {
	int n;
	if (argc == 1)
		printf("Defaulting to 1000 rays.\n");
	else if (argc != 2)
		printf("Unexpected number of arguments, ignoring after the first\n");

	n = argc > 1 ? atoi(argv[1]) : 1000;
    unsigned long t;
    struct timeval tv;
    MTRand myrng;

    Ray3D the_ray;
    int killed = 0;
    int32_t cntr_detected;
    int maxScatters = 20;   // Maximum allowed number of scatters
    Surface3D sample;
    NBackWall plate;
    AnalytSphere sphere;
    int sample_index = 0, plate_index = 1, sphere_index = 2;
    double c[3] = {0, 1, 0};
    SourceParam source;
    source.pinhole_r = 1e-3;
    for (int i = 0; i < 3; i++) {
    	source.pinhole_c[i] = c[i];
    }
    source.theta_max = 1/100;
    source.init_angle = 0;
    source.source_model = 1;
    source.sigma = 1/100;

    // A material to use
    Material standard_mat;
    standard_mat.name = "diffuse";
    standard_mat.func_name = "cosine";
    standard_mat.params = 0;
    standard_mat.n_params = 0;
    standard_mat.func = distribution_by_name("cosine");

    // Counter for the number of detected
    cntr_detected = 0;
    int32_t * numScattersRay;

    double aperture_c[3] = {0, 0};
    double aperture_axes[2] = {0.1, 0.1};
    plate.surf_index = plate_index;
    plate.n_detect = 1;
    plate.aperture_c = aperture_c;
    plate.aperture_axes = aperture_axes;
    plate.circle_plate_r = 2;
    plate.plate_represent = 1;
    plate.material = standard_mat;
    make_basic_sample(sample_index, 10, &sample);
    generate_empty_sphere(sphere_index, &sphere);
    sphere.sphere_c[0] = 0;
    sphere.sphere_c[1] = 0.1;
    sphere.sphere_c[2] = 0;
    sphere.sphere_r = 0.1;
    sphere.make_sphere = 0;

    // Seed the random number generator with the current time
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;

    // Set up the MTwister random number generator
    seedRand(t, &myrng);

    printf("The ray:\n");
    print_ray(&the_ray);
    printf("The sphere:\n");
    print_sphere(&sphere);
    printf("The plate:\n");
    print_nBackWall(&plate);
    printf("The sample:\n");
    print_surface(&sample);

    numScattersRay = (int32_t *)calloc(n, sizeof(int32_t));
    generating_rays_simple_pinhole(source, n, &killed, &cntr_detected, maxScatters, sample,
    		plate, sphere, &myrng, numScattersRay);

    printf("Number of detected rays is: %i\n", cntr_detected);
    printf("Sample is set-up to be specular, all of them should be detected.\n");

    free(numScattersRay);
    clean_up_surface_all_arrays(&sample);
}


