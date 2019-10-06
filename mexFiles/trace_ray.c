/*
 * Copyright (c) 2018-19, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 */

#include "mex.h"
#include "tracing_functions.h"
#include "ray_tracing_structs3D.h"
#include "common_helpers.h"
#include "small_functions3D.h"
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

/* 
 * For the simple model of the pinhole plate as a circle.
 * 
 * Trace a single ray returns 1 or 0 depending on if the ray is detected. 
 */
int32_t trace_ray_simple(Ray3D *the_ray, int *killed, int *cntr_detected, 
        int maxScatters, Surface3D Sample, BackWall Plate, AnalytSphere the_sphere, 
        gsl_rng *my_rng) {
    int n_allScatters;
    int dead;
        
    /* 
     * The total number of scattering events undergone (sample and pinhole 
     * plate) 1000 events are allowed in total. A separate limit is placed
     * on the number of scattering events off of the sample.
     */
    n_allScatters = 0;
    
    /*
     * If the ray is 'dead' it is no longer of interest. It may be dead in two 
     * ways, either by not meeting any surface of by going into the detector. 
     *      dead = 1, it did not meet any surfaces 
     *      dead = 2, it went into the detector
     */
    dead = 0;
    
    /* 
     * Keep propagating the ray until it is deemed 'dead', by either not 
     * intersecting either the sample or the pinhole plate.
     */
    while (!dead) {        
        /* The ray is dead unless we hit something */
        dead = 1;
                
        /******************************************************************/
        /* 
        * Try to scatter of sample. This only tries to scatter off of the 
        * sample and not the pinhole plate, must be done first as the newly
        * generated rays will in about half of the cases hit the inside of
        * the pinhole plate, which will cause weird (and very wrong) 
        * results.
        * 
        * dead = 1 by scatterOffSurface -- did not hit the sample
        * dead = 0 by scatterOffSurface -- scattered off the sample
        * 
        * If the ray has not hit the sample then it is immediately dead.
        */
        if (the_ray->nScatters == 0) {
            dead = scatterOffSurface(the_ray, &Sample, the_sphere, my_rng);
            if (dead == 0) {
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
        if ((the_ray->nScatters > maxScatters) || (n_allScatters > 15)) {
            /* Ray has exceeded the maximum number of scatters, kill it */
            the_ray->nScatters = -1;
            *killed += 1;
            break;
        }
        
        /* Try to scatter of both surfaces. */
        dead = scatterSimpleSurfaces(the_ray, &Sample, Plate, the_sphere, my_rng);

        /******************************************************************/
        /* Update counters */
        
        switch (dead) {
            case 2:
                /* Detected */
                *cntr_detected += 1;
                break;
            case 1:
                /* Did not hit a surface or get detected */
                break;
            case 0:
                /* Hit a surface */
                n_allScatters++;
                
                /* Hit the sample */
                if ((the_ray->on_surface == Sample.surf_index) || 
                        (the_ray->on_surface == the_sphere.surf_index)) {
                    the_ray->nScatters += 1;
                }
                break;
        }
    }
    if (dead == 2)
        return(1);
    else
        return(0);
}

/* 
 * For a simple model of the pinhole plate as a circle with multiple detectors.
 * 
 * Traces a single ray
 */
int32_t trace_ray_simpleMulti(Ray3D *the_ray, int *killed, int *cntr_detected, 
        int maxScatters, Surface3D Sample, NBackWall Plate, AnalytSphere the_sphere,
        gsl_rng *my_rng, int *detector) {
    int n_allScatters;
    int dead;
        
    /* 
     * The total number of scattering events undergone (sample and pinhole 
     * plate) 1000 events are allowed in total. A separate limit is placed
     * on the number of scattering events off of the sample.
     */
    n_allScatters = 0;
    
    /*
     * If the ray is 'dead' it is no longer of interest. It may be dead in two 
     * ways, either by not meeting any surface of by going into the detector. 
     *      dead = 1, it did not meet any surfaces 
     *      dead = 2, it went into the detector
     */
    dead = 0;
    
    /* 
     * Keep propagating the ray until it is deemed 'dead', by either not 
     * intersecting either the sample or the pinhole plate.
     */
    while (!dead) {        
        /* The ray is dead unless we hit something */
        dead = 1;
                
        /******************************************************************/
        /* 
        * Try to scatter of sample. This only tries to scatter off of the 
        * sample and not the pinhole plate, must be done first as the newly
        * generated rays will in about half of the cases hit the inside of
        * the pinhole plate, which will cause weird (and very wrong) 
        * results.
        * 
        * dead = 1 by scatterOffSurface -- did not hit the sample
        * dead = 0 by scatterOffSurface -- scattered off the sample
        * 
        * If the ray has not hit the sample then it is immediately dead.
        */
        if (the_ray->nScatters == 0) {
            dead = scatterOffSurface(the_ray, &Sample, the_sphere, my_rng);
            if (dead == 0) {
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
        if ((the_ray->nScatters > maxScatters) || (n_allScatters > 15)) {
            /* Ray has exceeded the maximum number of scatters, kill it */
            the_ray->nScatters = -1;
            *killed += 1;
            break;
        }
        
        /* Try to scatter of both surfaces. */
        dead = scatterSimpleMulti(the_ray, &Sample, Plate, the_sphere, my_rng, 
            detector);

        /******************************************************************/
        /* Update counters */
        
        switch (dead) {
            case 2:
                /* Detected */
                cntr_detected[*detector - 1] += 1;
                break;
            case 1:
                /* Did not hit a surface or get detected */
                break;
            case 0:
                /* Hit a surface */
                n_allScatters++;
                
                /* Hit the sample */
                if ((the_ray->on_surface == Sample.surf_index) || 
                        (the_ray->on_surface == the_sphere.surf_index)) {
                    the_ray->nScatters += 1;
                }
                break;
        }
    }
    if (dead == 2)
        return(1);
    else
        return(0);
}

/* 
 * For representing the pinhole plate as a triangulated surface.
 * 
 * Trace a single ray 
 */
int32_t trace_ray_triagPlate(Ray3D *the_ray, int *killed, int *cntr_detected, int maxScatters,
        Surface3D Sample, Surface3D Plate, AnalytSphere the_sphere,
        double backWall[], gsl_rng *my_rng) {
    int n_allScatters;
    int dead;
    
    /* 
     * The total number of scattering events undergone (sample and pinhole 
     * plate) 1000 events are allowed in total. A separate limit is placed
     * on the number of scattering events off of the sample.
     */
    n_allScatters = 0;
    
    /*
     * If the ray is 'dead' it is no longer of interest. It may be dead in two 
     * ways, either by not meeting any surface of by going into the detector. 
     *      dead = 1, it did not meet any surfaces 
     *      dead = 2, it went into the detector
     */
    dead = 0;
    
    /* 
     * Keep propagating the ray until it is deemed 'dead', by either not 
     * intersecting either the sample or the pinhole plate.
     */
    while (!dead) {        
        /* The ray is dead unless we hit something */
        dead = 1;
        
        /******************************************************************/
        /* 
        * Try to scatter of sample. This only tries to scatter off of the 
        * sample and not the pinhole plate, must be done first as the newly
        * generated rays will in about half of the cases hit the inside of
        * the pinhole plate, which will cause weird (and very wrong) 
        * results.
        * 
        * dead = 1 by scatterOffSurface -- did not hit the sample
        * dead = 0 by scatterOffSurface -- scattered off the sample
        * 
        * If the ray has not hit the sample then it is immediately dead.
        */
        if (the_ray->nScatters == 0) {
            dead = scatterOffSurface(the_ray, &Sample, the_sphere, my_rng);
            
            if (dead == 0) {
                /* Hit the sample */
                the_ray->nScatters++;
                n_allScatters++;
            } else {
                /* Move onto the next ray */
                continue;
            }
        }
        
        /* The number of scattering events is set to -1 if we exceed the overall
         * number of scattering events. */
        if ((the_ray->nScatters > maxScatters) || (n_allScatters > 1000)) {
            /* Ray has exceeded the maximum number of scatters, kill it */
            the_ray->nScatters = -1;
            *killed += 1;
            break;
        }
        
        /* Try to scatter of both surfaces. */
        dead = scatterSurfaces(the_ray, &Sample, &Plate, the_sphere, 
            backWall, my_rng);
        
        /******************************************************************/
        /* Update counters */
        
        switch (dead) {
            case 2:
                /* Detected */
                *cntr_detected += 1;
                break;
            case 1:
                /* Did not hit a surface or get detected */
                break;
            case 0:
                /* Hit a surface */
                n_allScatters++;
                
                /* Hit the sample */
                if ((the_ray->on_surface == Sample.surf_index) || 
                        (the_ray->on_surface == the_sphere.surf_index)) {
                    the_ray->nScatters++;
                }
                break;
        }
    }
    if (dead == 2)
        return(1);
    else
        return(0);
}

/* 
 * For scattering a ray off a surface only so that the distribution may be
 * aquired from the scattering.
 * 
 * Trace a single ray
 */
void trace_ray_justSample(Ray3D *the_ray, int *killed, int maxScatters, Surface3D Sample,
        AnalytSphere the_sphere, gsl_rng *my_rng) {
    int dead;
    
    /* 
     * If the ray is 'dead' it is no longer of interest. Here it may obly be
     * dead by not hitting the single sample surface or being killed.
     * 
     * dead = 1 Did not hit anything
     * dead = 3 Was killed
     */
    dead = 0;
    
    while (!dead) {
        /* The ray is dead unless we hit something */
        dead = 1;
        
        /* Try to scatter off the sample */
        dead = scatterOffSurface(the_ray, &Sample, the_sphere, my_rng);
        if (dead == 0) {
            /* Hit the sample */
            the_ray->nScatters += 1;
        }
        
        /* The number of scattering events is set to -1 if we exceed the overall
         * number of scattering events. */
        if (the_ray->nScatters > maxScatters) {
            /* Ray has exceeded the maximum number of scatters, kill it */
            the_ray->nScatters = -1;
            *killed += 1;
            dead = 3;
            break;
        }
    }
}

/* 
 * Trace the full trajectory of a ray scattering off a single surface. Records
 * all positions and directions of the ray.
 * 
 */
/*int32_t trace_ray_trajectory() {
    
}*/
