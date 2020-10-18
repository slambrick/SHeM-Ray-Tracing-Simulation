/*
 * Copyright (c) 2020, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 */

/*
 * Using C ray generation and a CAD model of the pinhole plate with a single
 * detector.
 *
 * TODO: const the objects passed around
 */
void generating_rays_cad_pinhole(double *source_parameters, int nrays, int *killed,
		int *cntr_detected, int maxScatters, const Surface3D *sample, const Surface3D *plate,
		AnalytSphere* the_sphere, BackWall *backWall, MTRand *myrng, double *numScattersRay) {
	int i;

	// TODO: this will be where memory is moved to the GPU

    for (i = 0; i < nrays; i++) {
        Ray3D the_ray;
        int detected;

        the_ray = create_ray(source_parameters[0], &source_parameters[1],
            source_parameters[4], source_parameters[5], source_model,
            source_parameters[6], &myrng);

        trace_ray_triag_plate(&the_ray, killed, cntr_detected, maxScatters,
        		sample, plate, the_sphere, backWall, myrng, &detected);

        /*
         * Add the number of scattering events the ray has undergone to the
         * histogram. But only if it is detected.
         */
        if (detected)
            numScattersRay[the_ray.nScatters - 1]++;
    }

    // TODO: this is where memory is extracted from the GPU
}
