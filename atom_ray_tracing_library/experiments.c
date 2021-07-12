/*
 * Copyright (c) 2020, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 */

#include "trace_ray.h"
#include "ray_tracing_core3D.h"
#include "mex.h"

/*
 * Using C ray generation and a CAD model of the pinhole plate with a single
 * detector.
 *
 * TODO: const the objects passed around
 */
void generating_rays_cad_pinhole(SourceParam source, int nrays, int *killed,
		int * const cntr_detected, int maxScatters, Sample overall_sample, Surface3D plate,
		double const backWall[], MTRand * const myrng,
		int32_t * const numScattersRay) {
	int i;

	// TODO: this will be where memory is moved to the GPU

    for (i = 0; i < nrays; i++) {
        Ray3D the_ray;

        create_ray(&the_ray, &source, myrng);

        trace_ray_triag_plate(&the_ray, maxScatters, overall_sample, plate,
                backWall, myrng);

        /*
         * Add the number of scattering events the ray has undergone to the
         * histogram. But only if it is detected.
         */
        switch (the_ray.status) {
            case 2:
                numScattersRay[the_ray.nScatters - 1]++;
                *cntr_detected += 1;
                break;
            case 1:
                // The ray died naturally...
                break;
            case -1:
                *killed += 1;
                break;
            case 0:
                // This should not happen...
                printf("Warning, your ray didn't finish...\n");
                break;
        }
    }

    // TODO: this is where memory is extracted from the GPU
}

void generating_rays_simple_pinhole(SourceParam source, int n_rays, int * const killed,
        int * const cntr_detected, int maxScatters, Sample overall_sample, NBackWall plate,
        MTRand * const myrng, int32_t * const numScattersRay) {

    int i;
    // TODO: this will be where memory is moved to the GPU
    for (i = 0; i < n_rays; i++) {
        Ray3D the_ray;
        int ind;
        create_ray(&the_ray, &source, myrng);

        trace_ray_simple_multi(&the_ray, maxScatters, overall_sample, plate, myrng);

        /*
         * Add the number of scattering events the ray has undergone to the
         * histogram. But only if it is detected.
         */
        switch (the_ray.status) {
            case 2:
                ind = (the_ray.detector - 1)*maxScatters + (the_ray.nScatters - 1);
                numScattersRay[ind]++;
                cntr_detected[the_ray.detector - 1] += 1;
                break;
            case 1:
                // The ray died naturally...
                break;
            case -1:
                *killed += 1;
                break;
            case 0:
                // This should not happen...
                printf("Warning, your ray didn't finish...\n");
                break;
        }
    }

    // TODO: this is where memory is extracted from the GPU
}

void given_rays_simple_pinhole(Rays3D * const all_rays, int * killed,
        int * const cntr_detected, Sample overall_sample, NBackWall plate,
        int maxScatters, int32_t * const detected,
        int32_t * const which_detector, MTRand * const myrng) {
    int i;

    // TODO: this will be where memory is moved to the GPU

    for (i = 0; i < all_rays->nrays; i++) {
        //trace_ray_simple_multi(&all_rays->rays[i], maxScatters, overall_sample, plate, myrng);
        which_detector[i] = all_rays->rays[i].detector;
        switch (all_rays->rays[i].status) {
            case 2:
                detected[i] = 1;
                cntr_detected[all_rays->rays[i].detector - 1] += 1;
                break;
            case 1:
                // The ray died naturally...
                break;
            case -1:
                *killed += 1;
                break;
            case 0:
                // This should not happen...
                printf("Warning, your ray didn't finish...\n");
                break;
        }
    }

    // TODO: this is where memory is extracted from the GPU
}

void given_rays_cad_pinhole(Rays3D * const all_rays, int * const killed, int * const cntr_detected,
        Sample overall_sample, Surface3D plate, double const backWall[],
        int maxScatters, int32_t * const detected, MTRand * const myrng) {
    int i;

    // TODO: this will be where memory is moved to the GPU

    for (i = 0; i < all_rays->nrays; i++) {
        trace_ray_triag_plate(&all_rays->rays[i], maxScatters, overall_sample, plate,
                        backWall, myrng);

        switch (all_rays->rays[i].status) {
            case 2:
                detected[i] = 1;
                *cntr_detected += 1;
                break;
            case 1:
                // The ray died naturally...
                break;
            case -1:
                *killed += 1;
                break;
            case 0:
                // This should not happen...
                printf("Warning, your ray didn't finish...\n");
                break;
        }
    }

    // TODO: this is where memory is extracted from the GPU
}
