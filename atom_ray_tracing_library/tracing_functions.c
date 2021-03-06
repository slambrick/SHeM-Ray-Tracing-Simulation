/*
 * Copyright (c) 2018-2020, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 *
 * Functions for performing a single iteration of the ray tracing algorithm. One
 * ray is scattered once. There are different ways to do that depending on the
 * simulation being run.
 */
#include "ray_tracing_core3D.h"
#include "intersect_detection3D.h"
#include "distributions3D.h"
#include <math.h>
#include "mtwister.h"
#include <stdbool.h>


/*
 * Scatters the given ray off a single triangulated surface, the sample, and an
 * analytic sphere, if that is desired. Returns true if the ray did not hit the
 * surface, i.e. it is 'dead' and returns false if the ray does, i.e. the ray
 * is 'alive'. This function has undergone some low level optimisation, so it
 * may not be written in the most intuitive and simple manner.
 *
 * NOTE: this function run by itself does not cause seg faults
 *
 * INPUTS:
 *  the_ray - pointer to a ray struct
 *  sample - pointer to a const Surface3D struct of the sample
 *  the_sphere - pointer to a const AnalytSphere struct of the sphere
 *  myrng - pointer to a random number generator object
 *  status - int 1, 0 declaring whether the ray is 'dead', 1 is dead (has not
 *         met), 0 is alive (has met)
 */
void scatterOffSurface(Ray3D * the_ray, Surface3D sample, AnalytSphere the_sphere,
        MTRand * const myrng) {
    double min_dist;
    int meets = 0;
    int tri_hit = -1;
    double nearest_n[3];
    double nearest_inter[3];
    double new_direction[3];
    bool meets_sphere = 0;
    int which_surface = -1;

    /* Much further than any of the triangles */
    min_dist = 10.0e10;

    /* Try to scatter of the sample */
    scatterTriag(the_ray, sample, &min_dist, nearest_inter, nearest_n, &meets,
        &tri_hit, &which_surface);

    /* Should the sphere be represented */
    if (the_sphere.make_sphere) {
        /* Check the sphere if we are not on it */
        if (the_ray->on_surface != the_sphere.surf_index) {
            scatterSphere(the_ray, the_sphere, &min_dist, nearest_inter,
            	nearest_n, &tri_hit, &which_surface, &meets_sphere);
        }
    }

    /* If we have met a triangle/sphere we must scatter off of it */
    if (meets || meets_sphere) {
        Material const * composition;

        if (meets_sphere) {
            /* sphere is defined to be uniform */
            composition = &(the_sphere.material);
        } else {
            composition = sample.compositions[tri_hit];
        }

        /* Find the new direction and update position*/
        composition->func(nearest_n, the_ray->direction,
            new_direction, composition->params, myrng);
        update_ray_direction(the_ray, new_direction);
        update_ray_position(the_ray, nearest_inter);

        /* Updates the current triangle and surface the ray is on */
        the_ray->on_element = tri_hit;
        the_ray->on_surface = which_surface;
    }

    the_ray->status = !(meets || meets_sphere);
}


/*
 * Scatters the ray off only the pinhole plate.
 *
 * INPUTS:
 *
 * OUTPUTS:
 *  dead - int 2, 1, 0, declaring whether the ray is dead. 1 is dead (has not
 *         met), 0 is alive (has met), 2 is detected (has hit the detector
 *         surface)
 */
void scatterPinholeSurface(Ray3D * the_ray, Surface3D plate, double const backWall[],
        MTRand * const myrng) {

    double min_dist;
    int meets;
    int tri_hit;
    double nearest_n[3];
    double nearest_inter[3];
    double new_direction[3];
    int which_surface;

    /* tri_hit stores which triangle has been hit */
    tri_hit = -1;

    /* which_surface stores which surface has been hit */
    which_surface = -1;

    /* meets is 0/1, have we met a triangle */
    meets = 0;

    /* Much further than any of the triangles */
    min_dist = 10.0e10;

    /* Try to scatter off the pinhole plate */
    scatterTriag(the_ray, plate, &min_dist, nearest_inter, nearest_n, &meets,
        &tri_hit, &which_surface);

    /* Update position/direction etc. */
    if (meets) {
        Material * composition;

        /* Get the scattering process */
        composition = plate.compositions[tri_hit];

        /* Update the direction and position of the ray */
        composition->func(nearest_n, the_ray->direction,
            new_direction, composition->params, myrng);
        update_ray_direction(the_ray, new_direction);
        update_ray_position(the_ray, nearest_inter);

        /* Updates the current triangle the ray is on */
        the_ray->on_element = tri_hit;

        /* Update the current surface the ray is on */
        the_ray->on_surface = which_surface;
    } else {
        /*
         * We must consider if the ray has been detected if it hasn't hit
         * either surface
         */

        /* First consider if the ray is propagating in the +ve y direction */
        if (the_ray->direction[1] > 0) {
            double alpha;
            double wall_hit[3];

            /* Find where, just behind the pinhole plate, the ray will hit,
             * backWall[0] is the y coordinate of the back of the pinhole plate
             */
            alpha = (backWall[0] - the_ray->position[1])/the_ray->direction[1];

            propagate(the_ray->position, the_ray->direction, alpha, wall_hit);

            /*
             * Now find if this point is covered by the plate, if it is then the
             * ray must be detected, if it is not, then the ray must be dead.
             * backWall[1] and backWall[2] are the depth in x and z of the
             * pinhole plate
             */
            if ((fabs(wall_hit[0]) < (backWall[1]/2)) && (fabs(wall_hit[2]) <
                    (backWall[2]/2))) {
                /* We update position only and keep the direction the same */
                update_ray_position(the_ray, wall_hit);

                // Detection is status = 2
                the_ray->status = 2;
                return;
            }
        } else {
            /* Ray be dead */
            meets = 0;
        }
    }

    the_ray->status = !meets;
}


/*
 * Scatters the ray off of two surfaces, one of the sample and one of the
 * pinhole plate. The pinhole plate surface includes a detection surface.
 *
 * INPUTS:
 *
 * OUTPUTS:
 *  dead - int 2, 1, 0, declaring whether the ray is dead. 1 is dead (has not
 *         met), 0 is alive (has met), 2 is detected (has hit the detector
 *         surface)
 */
void scatterSurfaces(Ray3D * the_ray, Surface3D sample, Surface3D plate,
		AnalytSphere the_sphere, double const backWall[], MTRand * const myrng) {

    double min_dist;
    int meets;
    int tri_hit;
    double nearest_n[3];
    double nearest_inter[3];
    double new_direction[3];
    bool meets_sphere;
    int which_surface;

    /* tri_hit stores which triangle has been hit */
    tri_hit = -1;

    /* which_surface stores which surface has been hit */
    which_surface = -1;

    /* meets is 0/1 have we met a triangle */
    meets = 0;
    meets_sphere = 0;

    /* Much further than any of the triangles */
    min_dist = 10.0e10;

    /* Try to scatter off the sample */
    scatterTriag(the_ray, sample, &min_dist, nearest_inter, nearest_n, &meets,
        &tri_hit, &which_surface);

    /* Try to scatter off the pinhole plate */
    scatterTriag(the_ray, plate, &min_dist, nearest_inter, nearest_n, &meets,
        &tri_hit, &which_surface);

    /* Should the sphere be represented */
    if (the_sphere.make_sphere) {
        /* Check the sphere if we are not on it */
        if (the_ray->on_surface != the_sphere.surf_index) {
            scatterSphere(the_ray, the_sphere, &min_dist, nearest_inter,
            	nearest_n, &tri_hit, &which_surface, &meets_sphere);
        }
    }

    /* Update position/direction etc. */
    if (meets || meets_sphere) {
        Material const * composition;

        if (meets_sphere) {
            /* sphere is defined to be uniform */
            composition = &(the_sphere.material);
        } else {
            if (which_surface == plate.surf_index) {
                composition = plate.compositions[tri_hit];
            } else {
                composition = sample.compositions[tri_hit];
            }
        }

        /* Find the new direction and update position*/
        composition->func(nearest_n, the_ray->direction,
            new_direction, composition->params, myrng);
        update_ray_direction(the_ray, new_direction);
        update_ray_position(the_ray, nearest_inter);

        /* Updates the current triangle and surface the ray is on */
        the_ray->on_element = tri_hit;
        the_ray->on_surface = which_surface;
    } else {
        /*
         * We must consider if the ray has been detected if it hasn't hit
         * either surface
         */
        /* First consider if the ray is propagating in the +ve y direction */
        if (the_ray->direction[1] > 0) {
            double alpha;
            double wall_hit[3];

            /* Find where, just behind the pinhole plate, the ray will hit,
             * backWall[0] is the y coordinate of the back of the pinhole plate
             */
            alpha = (backWall[0] - the_ray->position[1])/the_ray->direction[1];

            propagate(the_ray->position, the_ray->direction, alpha, wall_hit);

            /*
             * Now find if this point is covered by the plate, if it is then the
             * ray must be detected, if it is not, then the ray must be dead.
             * backWall[1] and backWall[2] are the depth in x and z of the
             * pinhole plate
             */
            if ((fabs(wall_hit[0]) < (backWall[1]/2)) && (fabs(wall_hit[2]) <
                    (backWall[2]/2))) {
                /* We update position only and keep the direction the same */
                update_ray_position(the_ray, wall_hit);

                the_ray->status = 2;
                return;
            }
        } else {
            /* Ray be dead */
            meets = 0;
        }
    }

    the_ray->status = !(meets || meets_sphere);
}

/*
 * Scatters the ray off a sample triangulated surface and a simple flat model of
 * the pinhole plate.
 *
 * INPUTS:

 *
 * OUTPUTS:
 *  dead - int 2, 1, 0, declaring whether the ray is dead. 1 is dead (has not
 *         met), 0 is alive (has met), n is detected from the detector (n-1)
 */
void scatterSimpleMulti(Ray3D * the_ray, Surface3D sample, NBackWall plate,
		AnalytSphere the_sphere, int * detector, MTRand * const myrng) {

    double nearest_n[3];
    double nearest_inter[3];
    double new_direction[3];

    /* tri_hit stores which triangle has been hit */
    int tri_hit = -1;

    /* By default no detection */
    int detected = 0;

    /* which_surface stores which surface has been hit */
    int which_surface = -1;

    /* By default don't hit the sphere */
    bool meets_sphere = false;

    /* meets is false/true have we met a triangle */
    int meets = 0;

    /* Much further than any of the triangles */
    double min_dist = 10.0e10;

    /* Try to scatter off the sample */
    int meets_sample = 0;
    scatterTriag(the_ray, sample, &min_dist, nearest_inter, nearest_n, &meets_sample,
        &tri_hit, &which_surface);
    meets = meets_sample;

    /* Should the sphere be represented */
    if (the_sphere.make_sphere) {
        /* Check the sphere if we are not on it */
        if (the_ray->on_surface != the_sphere.surf_index) {
            scatterSphere(the_ray, the_sphere, &min_dist, nearest_inter,
            	nearest_n, &tri_hit, &which_surface, &meets_sphere);
        }
    }

    /* Try to scatter off the simple pinhole plate */
    if (the_ray->on_surface != plate.surf_index) {
        int meets_wall = 0;
        multiBackWall(the_ray, plate, &min_dist, nearest_inter, nearest_n,
        	&meets_wall, &tri_hit, &which_surface, &detected);
        meets = meets || meets_wall;
    }

    /* If we are detected */
    if (detected) {
        /* Update the ray position but not the direction */
        update_ray_position(the_ray, nearest_inter);
        *detector = detected;

        /* 2 = detected ray */
        the_ray->status = 2;
        return;
    }

    /* Update position/direction etc. */
    if (meets || meets_sphere) {
        Material const * composition;

        if (meets_sphere) {
            /* sphere is defined to be uniform */
            composition = &(the_sphere.material);
        } else {
            if (which_surface == plate.surf_index)
                composition = &(plate.material);
            else
                composition = sample.compositions[tri_hit];
        }

        /* Find the new direction and update position*/
        composition->func(nearest_n, the_ray->direction,
            new_direction, composition->params, myrng);
        /* Updates the current triangle and surface the ray is on */
        the_ray->on_element = tri_hit;
        the_ray->on_surface = which_surface;
        update_ray_direction(the_ray, new_direction);
        update_ray_position(the_ray, nearest_inter);
    } else {
        /* Ray be dead */
        meets = 0;
    }

    the_ray->status = !(meets || meets_sphere);
}


/*
 * Detects on a hemisphere centered on the sample plane, the aperture is
 * specified by two angles and a half cone angle.
 */
//void scatterAbstractSurfaces(Ray3D *the_ray, Surface3D const * sample, AbstractHemi const * plate,
//		AnalytSphere const * the_sphere, MTRand * const myrng, int * status) {
//
//    int meets = 0;
//
//    *status =!meets;
//}
