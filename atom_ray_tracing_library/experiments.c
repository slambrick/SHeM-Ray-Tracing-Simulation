/*
 * Copyright (c) 2020-21, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 * 
 * Each of the functions in this file performs a single "experiment" for a large
 * number of rays.
 * 
 * TODO: there is repeated code here, it can probably be simplified somewhat
 */

#include "trace_ray.h"
#include "ray_tracing_core3D.h"
#include "mex.h"

void record_ray_result(Ray3D the_ray, int * const cntr_detected, int * killed, 
        int32_t * const numScattersRay, int maxScatters);


/*
 * Records the appropriate results from the ray. Note that this is for the
 * outputs associated with C-based source models. The set of outputs are
 * different for the "given ray" experiments, they should not use this function.
 * 
 * INPUTS: 
 *  the_ray        - a ray that has finished its trajectory
 *  cntr_detected  - an array of the number of detections in each detector
 *  killed         - a counter of the number of rays that needed to be stopped
 *  numScattersRay - histogram of the number of scattering events detected rays
 *                   underwent during their trajectory
 */
void record_ray_result(Ray3D the_ray, int * const cntr_detected, int * killed, 
        int32_t * const numScattersRay, int maxScatters) {
    int ind;
    
    // Add the number of scattering events the ray has undergone to the
    // histogram. But only if it is detected.
    switch (the_ray.status) {
        case 2:
            // The ray is detected
            ind = (the_ray.detector - 1)*maxScatters + (the_ray.nScatters - 1);
            numScattersRay[ind]++;
            cntr_detected[the_ray.detector - 1] += 1;
            break;
        case 1:
            // The ray died naturally...
            break;
        case -1:
            // We had to stop the ray because it scattered too many times.
            *killed += 1;
            break;
        case 0:
            // This should not happen...
            printf("Warning, your ray didn't finish...\n");
            break;
    }
}


void single_start_ray_cad_pinhole(Ray3D start_ray, int nrays, int * killed,
        int * const cntr_detected, int maxScatters, Sample overall_sample,
        Surface3D plate, double const backWall[], MTRand * const myrng,
        int32_t * const numScattersRay) {
    int i;
    double nearest_n[3];
    double nearest_b[6];

    // Perform the first scattering event
    scatterOffSurfaceReturn(&start_ray, overall_sample, myrng, nearest_n, 
                            nearest_b);
    
    if (start_ray.status == 1) {
        // Initial ray misses the sample, 0 signal
        return;
    } else {
        start_ray.nScatters++;
    }
    
    // Copy the ray position, generate a direction from the sample material then scatter
    for (i = 0; i < nrays; i++) {
        Ray3D the_ray;
        Material const * composition;
        double new_direction[3];

        copy_ray(&the_ray, start_ray);

        // Get the correct distribution and generate a direction
        if (the_ray.on_surface == overall_sample.the_sphere->surf_index) {
            composition = &(overall_sample.the_sphere->material);
        } else if (the_ray.on_surface == overall_sample.the_circle->surf_index) {
            composition = &(overall_sample.the_sphere->material);
        } else {
            composition = overall_sample.triag_sample->compositions[the_ray.on_element];
        }
        
        /* Find the new direction and update position*/
        composition->func(nearest_n, nearest_b, the_ray.direction,
            new_direction, composition->params, myrng);
        update_ray_direction(&the_ray, new_direction);

        // Trace the rest of the ray trajectory
        trace_ray_triag_plate(&the_ray, maxScatters, overall_sample, plate,
                backWall, myrng);

        record_ray_result(the_ray, cntr_detected, killed, numScattersRay, maxScatters);
    }
}

void single_start_ray_simple_pinhole(Ray3D start_ray, int nrays, int * killed, 
        int * const cntr_detected, int maxScatters, Sample overall_sample, 
        NBackWall plate, MTRand * const myrng, int32_t * const numScattersRay) {
    int i;
    double nearest_n[3];
    double nearest_b[6];

    // Perform the first scattering event
    scatterOffSurfaceReturn(&start_ray, overall_sample, myrng, nearest_n, nearest_b);

    if (start_ray.status == 1) {
        // Initial ray misses the sample, 0 signal
        return;
    } else {
        start_ray.nScatters++;
    }

    // Copy the ray position, generate a direction from the sample material then scatter
    for (i = 0; i < nrays; i++) {
        Ray3D the_ray;
        Material const * composition;
        double new_direction[3];

        copy_ray(&the_ray, start_ray);

        // Get the correct distribution and generate a direction
        if (the_ray.on_surface == overall_sample.the_sphere->surf_index) {
            composition = &(overall_sample.the_sphere->material);
        } else if (the_ray.on_surface == overall_sample.the_circle->surf_index) {
            composition = &(overall_sample.the_sphere->material);
        } else {
            composition = overall_sample.triag_sample->compositions[the_ray.on_element];
        }
        
        /* Find the new direction and update position*/
        composition->func(nearest_n, nearest_b, the_ray.direction,
            new_direction, composition->params, myrng);
        update_ray_direction(&the_ray, new_direction);

        // Trace the rest of the ray trajectory
        trace_ray_simple_multi(&the_ray, maxScatters, overall_sample, plate, myrng);

        record_ray_result(the_ray, cntr_detected, killed, numScattersRay, maxScatters);
    }
}


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

    for (i = 0; i < nrays; i++) {
        Ray3D the_ray;

        create_ray(&the_ray, &source, myrng);

        trace_ray_triag_plate(&the_ray, maxScatters, overall_sample, plate,
                backWall, myrng);

        record_ray_result(the_ray, cntr_detected, killed, numScattersRay, maxScatters);
    }
}

void generating_rays_simple_pinhole(SourceParam source, int nrays, int * const killed,
        int * const cntr_detected, int maxScatters, Sample overall_sample, NBackWall plate,
        MTRand * const myrng, int32_t * const numScattersRay) {

    if (source.source_model == 3) {
        Ray3D start_ray;
        create_ray(&start_ray, &source, myrng);

        single_start_ray_simple_pinhole(start_ray, nrays, killed, cntr_detected, maxScatters, 
            overall_sample, plate, myrng, numScattersRay);
    } else {
        int i;
        
        for (i = 0; i < nrays; i++) {
            Ray3D the_ray;
            
            create_ray(&the_ray, &source, myrng);

            trace_ray_simple_multi(&the_ray, maxScatters, overall_sample, plate, myrng);

            record_ray_result(the_ray, cntr_detected, killed, numScattersRay, maxScatters);
        }
    }
}

void generating_rays_abstract_pinhole(SourceParam source, int nrays, int * const killed,
        int * const cntr_detected, int maxScatters, Sample overall_sample, AbstractHemi plate,
        MTRand * const myrng, int32_t * const numScattersRay) {

    if (source.source_model == 3) {
        Ray3D start_ray;
        create_ray(&start_ray, &source, myrng);

        // TODO: this
        //single_start_ray_abstract_pinhole(start_ray, nrays, killed, cntr_detected, maxScatters, 
        //    overall_sample, plate, myrng, numScattersRay);
    } else {
        int i;
        
        for (i = 0; i < nrays; i++) {
            Ray3D the_ray;
            
            create_ray(&the_ray, &source, myrng);

            trace_ray_abstract(&the_ray, maxScatters, overall_sample, plate, myrng);

            record_ray_result(the_ray, cntr_detected, killed, numScattersRay, maxScatters);
        }
    }
}

void given_rays_simple_pinhole(Rays3D * const all_rays, int * killed,
        int * const cntr_detected, Sample overall_sample, NBackWall plate,
        int maxScatters, int32_t * const detected,
        int32_t * const which_detector, MTRand * const myrng) {
    int i;

    for (i = 0; i < all_rays->nrays; i++) {
        trace_ray_simple_multi(&all_rays->rays[i], maxScatters, overall_sample, plate, myrng);
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
}

void given_rays_cad_pinhole(Rays3D * const all_rays, int * const killed, int * const cntr_detected,
        Sample overall_sample, Surface3D plate, double const backWall[],
        int maxScatters, int32_t * const detected, MTRand * const myrng) {
    int i;

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
}

