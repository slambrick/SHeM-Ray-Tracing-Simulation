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

void single_ray_test();

int main() {
    single_ray_test();
}


void single_ray_test() {
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
    int detector = 0;
    int detected = 0;
    int sample_index = 0, plate_index = 1, sphere_index = 2;

    // A material to use
    Material standard_mat;
    standard_mat.name = "diffuse";
    standard_mat.func_name = "cosine";
    standard_mat.params = 0;
    standard_mat.n_params = 0;
    standard_mat.func = distribution_by_name("cosine");

    // Counter for the number of detected
    cntr_detected = 0;

    double init_pos[3] = {0, 1, 0};
    double init_dir[3] = {0, -1, 0};

    // Lets create a ray
    new_Ray(&the_ray, init_pos, init_dir);

    double aperture_c[3] = {0, 0};
    double aperture_axes[2] = {0.1, 0.1};
    plate.surf_index = plate_index;
    plate.n_detect = 1;
    plate.aperture_c = aperture_c;
    plate.aperture_axes = aperture_axes;
    plate.circle_plate_r = 2;
    plate.plate_represent = 1;
    plate.material = standard_mat;
    make_basic_sample(sample_index, 2, &sample);
    generate_empty_sphere(sphere_index, &sphere);
    sphere.sphere_c[0] = 0;
    sphere.sphere_c[1] = 0.1;
    sphere.sphere_c[2] = 0;
    sphere.sphere_r = 0.1;
    sphere.make_sphere = 1;

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

    if (1) {
        int n = 1e7;
        for (int i = 0; i < n; i++) {
            Ray3D this_ray = the_ray;
            trace_ray_simple_multi(&this_ray, &killed, &cntr_detected, maxScatters,
                    &sample, &plate, &sphere, &detector, &myrng, &detected);
        }
    } else {
        int status;
        scatterOffSurface(&the_ray, &sample, &sphere, &myrng, &status);
        printf("\nstatus: %i\n", status);
        scatterOffSurface(&the_ray, &sample, &sphere, &myrng, &status);
        printf("\nstatus: %i\n", status);
    }
    // TODO: test results

    printf("Number of detected rays is: %i\n", cntr_detected);

    clean_up_surface_all_arrays(&sample);
}


