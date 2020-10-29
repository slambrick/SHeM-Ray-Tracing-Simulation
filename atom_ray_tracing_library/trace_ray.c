/*
 * Copyright (c) 2018-20, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 */

#include "trace_ray.h"
#include "tracing_functions.h"
#include "ray_tracing_core3D.h"
#include "common_helpers.h"
#include <math.h>
#include <stdlib.h>
#include "mtwister.h"

/*
 * For a simple model of the pinhole plate as a circle with multiple detectors.
 *
 * Traces a single ray
 *
 * NOTE: This function run by itself does cause seg faults
 * TODO: find the basterd pointer that causes this!
 */
void trace_ray_simple_multi(Ray3D *the_ray,
        int maxScatters, Surface3D sample, NBackWall plate,
		AnalytSphere the_sphere, MTRand * const myrng) {
    /*
     * The total number of scattering events undergone (sample and pinhole
     * plate) 1000 events are allowed in total. A separate limit is placed
     * on the number of scattering events off of the sample.
     */
    int n_allScatters = 0;
    int detector = 0;

    /*
     * Keep propagating the ray until it is deemed 'dead', by either not
     * intersecting either the sample or the pinhole plate.
     */
    while (!(the_ray->status)) {
        /* The ray is dead unless we hit something */
        the_ray->status = 1;

        /*
        * Try to scatter of sample. This only tries to scatter off of the
        * sample and not the pinhole plate.
        *
        * If the ray has not hit the sample then it is immediately dead.
        */
        if (the_ray->nScatters == 0) {
            scatterOffSurface(the_ray, sample, the_sphere, myrng);
            if (!the_ray->status) {
                /* Hit the sample */
                the_ray->nScatters += 1;
                n_allScatters++;
            } else {
                /* Move onto the next ray */
                continue;
            }
        }

        /* The number of scattering events is set to -1 if we exceed the overall
         * number of scattering events. */
        if ((the_ray->nScatters > maxScatters) || (n_allScatters > 50)) {
            /* Ray has exceeded the maximum number of scatters, kill it */
            the_ray->nScatters = -1;
            the_ray->status = -1;
            break;
        }


        /* Try to scatter of both surfaces. */
        // TODO: change of status inside escatterSimpleMulti
        scatterSimpleMulti(the_ray, sample, plate, the_sphere, &detector, myrng);

        if (the_ray->status == 2) {
            the_ray->detector = detector;
            break;
        }
        /******************************************************************/
        /* Update counters */

        if (the_ray->status == 0) {
            /* Hit a surface */
            n_allScatters++;

            /* Hit the sample */
            if ((the_ray->on_surface == sample.surf_index) ||
                    (the_ray->on_surface == the_sphere.surf_index)) {
                the_ray->nScatters += 1;
            }
        }
    }
}

/*
 * For representing the pinhole plate as a triangulated surface.
 *
 * Trace a single ray
 */
void trace_ray_triag_plate(Ray3D * the_ray, int maxScatters,
        Surface3D sample, Surface3D plate, AnalytSphere the_sphere,
        double const backWall[], MTRand * const myrng) {
    int n_allScatters;

    /*
     * The total number of scattering events undergone (sample and pinhole
     * plate) 1000 events are allowed in total. A separate limit is placed
     * on the number of scattering events off of the sample.
     */
    n_allScatters = 1000;

    // Keep propagating the ray until it doesn't hit something
    while (!(the_ray->status)) {
        /* The ray is dead unless we hit something */
        the_ray->status = 1;

        /*
        * Try to scatter of sample. This only tries to scatter off of the
        * sample and not the pinhole plate.
        *
        * If the ray has not hit the sample then it is immediately dead.
        */
        if (the_ray->nScatters == 0) {
            scatterOffSurface(the_ray, sample, the_sphere, myrng);

            if (!(the_ray->status)) {
                /* Hit the sample */
                the_ray->nScatters++;
                n_allScatters++;
            } else {
                /* Move onto the next ray */
                break;
            }
        }

        /* The number of scattering events is set to -1 if we exceed the overall
         * number of scattering events. */
        if ((the_ray->nScatters > maxScatters) || (n_allScatters > 1000)) {
            /* Ray has exceeded the maximum number of scatters, kill it */
            the_ray->nScatters = -1;
            the_ray->status = -1;
            break;
        }

        /* Try to scatter of both surfaces. */
        scatterSurfaces(the_ray, sample, plate, the_sphere, backWall, myrng);

        /******************************************************************/
        /* Update counters */

        if (the_ray->status == 0) {
            /* Hit a surface */
            n_allScatters++;

            /* Hit the sample */
            if ((the_ray->on_surface == sample.surf_index) ||
                    (the_ray->on_surface == the_sphere.surf_index)) {
                the_ray->nScatters++;
            }
        }
    }
}

/*
 * For scattering a ray off a surface only. So that the distribution may be
 * acquired from the scattering.
 *
 * Trace a single ray off only ths ample
 *
 * TODO: change to use the new ray status inside the struct
 */
void trace_ray_just_sample(Ray3D * the_ray, int * const killed, int maxScatters,
        Surface3D sample, AnalytSphere the_sphere, MTRand * const myrng) {

    while (!(the_ray->status)) {
        /* The ray is dead unless we hit something */
        the_ray->status = 1;

        /* Try to scatter off the sample */
        scatterOffSurface(the_ray, sample, the_sphere, myrng);
        if (the_ray->status == 0) {
            /* Hit the sample */
            the_ray->nScatters += 1;
        }

        /* The number of scattering events is set to -1 if we exceed the overall
         * number of scattering events. */
        if (the_ray->nScatters > maxScatters) {
            /* Ray has exceeded the maximum number of scatters, kill it */
            the_ray->nScatters = -1;
            *killed += 1;
            break;
        }
    }
}

/*
 * Trace the full trajectory of a ray scattering off a single surface. Records
 * all positions and directions of the ray.
 * TODO: this for de-bugging
 */
/*int trace_ray_trajectory() {

}*/
