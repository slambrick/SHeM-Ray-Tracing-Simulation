/*
 * Copyright (c) 2020, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 */

#include "trace_ray.h"
#include "ray_tracing_core3D.h"

/*
 * Using C ray generation and a CAD model of the pinhole plate with a single
 * detector.
 *
 * TODO: const the objects passed around
 */
void generating_rays_cad_pinhole(SourceParam source, int nrays, int *killed,
		int * const cntr_detected, int maxScatters, Surface3D const * const sample, Surface3D const * const plate,
		AnalytSphere the_sphere, double const backWall[], MTRand * const myrng, int32_t * const numScattersRay) {
	int i;

	// TODO: this will be where memory is moved to the GPU

    for (i = 0; i < nrays; i++) {
        Ray3D the_ray;
        int detected;

        create_ray(&the_ray, &source, myrng);

        trace_ray_triag_plate(&the_ray, killed, cntr_detected, maxScatters,
        		sample, plate, &the_sphere, backWall, myrng, &detected);

        /*
         * Add the number of scattering events the ray has undergone to the
         * histogram. But only if it is detected.
         */
        if (detected)
            numScattersRay[the_ray.nScatters - 1]++;
    }

    // TODO: this is where memory is extracted from the GPU
}

void generating_rays_simple_pinhole(SourceParam source, int n_rays, int * const killed,
        int * const cntr_detected, int maxScatters, Surface3D const * const sample, NBackWall plate,
        AnalytSphere the_sphere, MTRand * const myrng, int32_t * const numScattersRay) {

    int i;
    // TODO: this will be where memory is moved to the GPU

    for (i = 0; i < n_rays; i++) {
        int detected = 0;
        int detector = 0;
        Ray3D the_ray;

        create_ray(&the_ray, &source, myrng);

        trace_ray_simple_multi(&the_ray, killed, cntr_detected, maxScatters,
                 sample, &plate, &the_sphere, &detector, myrng, &detected);

        /*
         * Add the number of scattering events the ray has undergone to the
         * histogram. But only if it is detected.
         */
        if (detected) {
            int ind = (detector - 1)*maxScatters + (the_ray.nScatters - 1);
            numScattersRay[ind]++;
        }
    }

    // TODO: this is where memory is extracted from the GPU
}
