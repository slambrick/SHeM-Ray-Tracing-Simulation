/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the  
 * GNU/GPL-3.0-or-later.
 * 
 * Contains the algorithm for scattering rays off a surface in 2D.
 */
#include "small_functions2D.h"
#include "ray_tracing_structs2D.h"
#include "intersect_detection2D.h"
#include "distributions2D.h"
#include <stdlib.h>
#include "mtwister.h"

/* 
 * Function that performs the simulation, scatters the rays stored in the Rays
 * struct off of the sample surface stored in the Surface struct. Loops through 
 * all the rays and traces them one at a time.
 * 
 * INPUTS:
 *  Sample   - a Surface2D struct containing the 1D surface to scatter the rays
 *             off of
 *  all_rays - A Rays2D struct containing all the rays to scatter off the sample
 */
void scatterRays2D(Surface2D Sample, Rays2D all_rays) {
    int i;
    
    /* Loop through and trace each ray */
    for (i = 0; i < all_rays.nrays; i++) {
        trace_ray2D(&all_rays.rays[i], Sample);
    }
    
    return;
}

/* 
 * Trace a single ray by scattering it off the sample surface.
 * 
 * INPUTS:
 *  the_ray - a Ray2D struct that contains the position, direction, and number of
 *            scattering events of a single ray
 *  Sample  - a Surface2D struct that contains the sample 1D surface to scatter
 *            the ray off
 */
void trace_ray2D(Ray2D *the_ray, Surface2D Sample) {
    int dead;
    int cnt;
        
    /* Trace the ray until it doesn't hit the sample again */
    cnt = 0;
    dead = 0;
    while (!dead) {
        double normal[2];
        double intersection[2];
        int nearest_element;
        
        /* Artifical limit on the number of scatters */
        if (cnt > 10000) {
            break;
        }
        
        cnt++;
        
        /* Try and hit the sample, if unable to then ray is dead */
        dead = !intersect2D(the_ray, Sample, intersection, normal, &nearest_element);
        
        /* 
         * Generate a new direction for the ray, update its position and 
         * direction, add one to the number of scattering events it has
         * undergone, and specify which element of surface it is now on.
         */
        if (!dead) {
            the_ray->position[0] = intersection[0];
            the_ray->position[1] = intersection[1];
            /* Aquire a new direction for the ray */
            new_direction2D(the_ray, normal, Sample.scattering, 
                Sample.scattering_parameters, Sample.my_rng);
            
            /* Increment the number of scattering events */
            the_ray->nscatters++;
            
            /* Update which surface element the ray is on */
            the_ray->on_element = nearest_element;
        }
    }
    
    return;
}

/* 
 * Tries to intersect a single ray with the sample surface
 * 
 * INPUTS:
 *  the_ray     - a Ray2D struct containing the position and direction of the ray
 *  Sample      - a Surface2D struct containing the smaple surface that is being 
 *                intersected with
 *  intersction - a two element double array to contain the position of any
 *                intersection
 *  normal      - a two element double array to contain the unit normal to the 
 *                surface at the point of intersection
 *  nearest_element - pointer to an int specifying the index of the surface 
 *                    element that is intersected with
 * 
 * OUTPUT:
 *  meets - int, 1 if the ray meets the surface, 0 if it does not
 */
int intersect2D(Ray2D *the_ray, Surface2D Sample, double intersection[2], 
              double normal[2], int *nearest_element) {
    int i;
    double distance;
    int meets;
    
    /* By default we do not meet anything */
    meets = 0;
    
    /* The initial distance is set to be negative */
    distance = -2;
    
    /* The initial index of the surface that the ray intersects is set to be none */
    *nearest_element = -1;
    
    /* Loop through each element of the 1D surface */
    for (i = 0; i < Sample.n_elements; i++) {
        double v1[2];
        double v2[2];
        double nn[2];
        int j;
        double m, c, t;
        double intersect_tmp[2];
        
        /* If the ray is on the surface element of this iteration then skip 
         * this iteration. */
        if (i == the_ray->on_element)
            continue;
        
        /* Extract the information about the ith element of surface */
        get_element2D(Sample, i, v1, v2, nn);
        
        /* Calculate the intersection point between the ray and the line defined
         * by the surface element. */
        
        /* Gradient and y intercept of the line of the surface element */
        m = (v2[1] - v1[1])/(v2[0] - v1[0]);
        c = v1[1] - m*v1[0];
        
        /* Distance the ray has to travel to intersect */
        t = (the_ray->position[1] - m*the_ray->position[0] - c) / 
            (m*the_ray->direction[0] - the_ray->direction[1]);
        
        /* Ray does not travel backwards */
        if (t < 0)
            continue;
        
        /* The intersection position */
        intersect_tmp[0] = the_ray->position[0] + the_ray->direction[0]*t;
        intersect_tmp[1] = the_ray->position[1] + the_ray->direction[1]*t;

        /* Check if the intersection is within the range of the element */
        if ((intersect_tmp[0] < v1[0]) || (intersect_tmp[0] > v2[0]))
            continue;
        
        /* 
         * Is the distance to the intersection smaller than any previous 
         * intersection, if so this is the best intersection so far, if the 
         * distance so far is negative then this intersection is the first 
         * successful one.
         */
        if ((t < distance) || (distance < 0)) {
            intersection[0] = intersect_tmp[0];
            intersection[1] = intersect_tmp[1];
            meets = 1;
            for (j = 0; j < 2; j++) 
                normal[j] = nn[j];
            distance = t;
            *nearest_element = i;
        }
    }
    /* Return if the ray has met a surface */
    return(meets);
}

