/* 
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 *
 * 
 * Functions for performing a single iteration of the ray tracing algorithm. One
 * ray is scattered once. There are different ways to do that depending on the
 * simulation being run. 
 */
#include "small_functions3D.h"
#include "ray_tracing_structs3D.h"
#include "scattering3D.h"
#include "scattering_processes3D.h"
#include <math.h>


/*
 * Scatters the given ray off a single triangulated surface, the sample, and an
 * analytic sphere, if that is desired. Returns true if the ray did not hit the 
 * surface, i.e. it is 'dead' and returns false if the ray does, i.e. the ray
 * is 'alive'. This function has undergone some low level optimisation, so it 
 * may not be written in the most intuitive and simple manner.
 * 
 * INPUTS:

 * 
 * OUTPUTS:
 *  dead - int 1, 0 declaring whether the ray is 'dead', 1 is dead (has not 
 *         met), 0 is alive (has met)
 */
int scatterOffSurface(Ray3D *the_ray, Surface3D *Sample, AnalytSphere 
        the_sphere, MTRand *myrng) {
    double min_dist;
    int meets;
    int tri_hit;
    double nearest_n[3];
    double nearest_inter[3];
    int meets_sphere;
    int which_surface;
    
    tri_hit = -1;
    meets = 0;
    meets_sphere = 0;
    
    /* Much further than any of the triangles */
    min_dist = 10.0e10;
    
    /* Try to scatter of the sample */
    scatterTriag(the_ray, Sample, &min_dist, nearest_inter, nearest_n, &meets,
        &tri_hit, &which_surface);
    
    /* Should the sphere be represented */
    if (the_sphere.make_sphere) {
        /* Check the sphere if we are not on it */
        if (the_ray->on_surface != the_sphere.surf_index) {
            meets_sphere = scatterSphere(the_ray, the_sphere, &min_dist, 
                nearest_inter, nearest_n, &tri_hit, &which_surface);
        }
    }
    
    /* If we have met a triangle/sphere we must scatter off of it */
    if (meets || meets_sphere) {
        double composition;
        double scattering_parameters;
        
        if (meets_sphere == 1) {
            /* sphere is defined to be uniform */
            composition = the_sphere.composition;
            scattering_parameters = the_sphere.scattering_parameters;
        } else {
            composition = Sample->composition[tri_hit];
            scattering_parameters = Sample->scattering_parameters[tri_hit];
        }
        
        /* Find the new direction and update position*/
        new_direction3D(the_ray, nearest_n, composition, scattering_parameters, myrng);
        update_ray_position(the_ray, nearest_inter);
        
        /* Updates the current triangle and surface the ray is on */
        the_ray->on_element = tri_hit;
        the_ray->on_surface = which_surface;
    }
    
    return(!(meets || meets_sphere));
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
int scatterPinholeSurface(Ray3D *the_ray, Surface3D *Plate, double backWall[],
        MTRand *myrng) {
    
    double min_dist;
    int meets;
    int tri_hit;
    double nearest_n[3];
    double nearest_inter[3];
    int which_surface;
    
    /* tri_hit stores which triangle has been hit */
    tri_hit = -1;
    
    /* which_surface stores which surface has been hit */
    which_surface = -1;
    
    /* meets is 0/1, have we met a triangle */
    meets = 0;
    
    /* Much further than any of the triangles */
    min_dist = 10.0e10;
    
    /* Try to scater off the pinhole plate */
    scatterTriag(the_ray, Plate, &min_dist, nearest_inter, nearest_n, &meets, 
        &tri_hit, &which_surface);
    
    /* Update position/direction etc. */
    if (meets) {
        double scattering;
        double scattering_parameters;

        /* Get the scattering proccess */
        scattering = Plate->composition[tri_hit];
        scattering_parameters = Plate->scattering_parameters[tri_hit];
        
        /* Update the direction and position of the ray */
        new_direction3D(the_ray, nearest_n, scattering, scattering_parameters, myrng);
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
        
        /* First consider if the ray is propogating in the +ve y direction */
        if (the_ray->direction[1] > 0) {
            double alpha;
            double wall_hit[3];
            
            /* Find where, just behind the pinhole plate, the ray will hit, 
             * backWall[0] is the y coordinate of the back of the pinhole plate
             */
            alpha = (backWall[0] - the_ray->position[1])/the_ray->direction[1];
            
            propogate(the_ray->position, the_ray->direction, alpha, wall_hit);
            
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
                
                return(2);
            }    
        } else {
            /* Ray be dead */
            meets = 0;
        }
    }
    
    return(!meets);
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
int scatterSurfaces(Ray3D *the_ray, Surface3D *Sample, Surface3D *Plate, 
        AnalytSphere the_sphere, double backWall[], MTRand *myrng) {
    
    double min_dist;
    int meets;
    int tri_hit;
    double nearest_n[3];
    double nearest_inter[3];
    int meets_sphere;
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
    scatterTriag(the_ray, Sample, &min_dist, nearest_inter, nearest_n, &meets,
        &tri_hit, &which_surface);
    
    /* Try to scater off the pinhole plate */
    scatterTriag(the_ray, Plate, &min_dist, nearest_inter, nearest_n, &meets,
        &tri_hit, &which_surface);
    
    /* Should the sphere be represented */
    if (the_sphere.make_sphere) {
        /* Check the sphere if we are not on it */
        if (the_ray->on_surface != the_sphere.surf_index) {
            meets_sphere = scatterSphere(the_ray, the_sphere, &min_dist, 
                nearest_inter, nearest_n, &tri_hit, &which_surface);
        }
    }
    
    /* Update position/direction etc. */
    if (meets || meets_sphere) {
        double composition;
        double scattering_parameters;
        
        if (meets_sphere) {
            /* sphere is defined to be uniform */
            composition = the_sphere.composition;
            scattering_parameters = the_sphere.scattering_parameters;
        } else {
            if (which_surface == Plate->surf_index) {
                composition = Plate->composition[tri_hit];
                scattering_parameters = Plate->scattering_parameters[tri_hit];
            } else {
                composition = Sample->composition[tri_hit];
                scattering_parameters = Sample->scattering_parameters[tri_hit];
            }
        }
        
        /* Find the new direction and update position*/
        new_direction3D(the_ray, nearest_n, composition, scattering_parameters, myrng);
        update_ray_position(the_ray, nearest_inter);
        
        /* Updates the current triangle and surface the ray is on */
        the_ray->on_element = tri_hit;
        the_ray->on_surface = which_surface;
    } else {
        /* 
         * We must consider if the ray has been detected if it hasn't hit 
         * either surface 
         */
        /* First consider if the ray is propogating in the +ve y direction */
        if (the_ray->direction[1] > 0) {
            double alpha;
            double wall_hit[3];
            
            /* Find where, just behind the pinhole plate, the ray will hit, 
             * backWall[0] is the y coordinate of the back of the pinhole plate
             */
            alpha = (backWall[0] - the_ray->position[1])/the_ray->direction[1];
            
            propogate(the_ray->position, the_ray->direction, alpha, wall_hit);
            
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
                
                return(2);
            }    
        } else {
            /* Ray be dead */
            meets = 0;
        }
    }
    
    return(!(meets || meets_sphere));
}


/*
 * Scatters the ray off a sample triangulated surface and a simple flat model of
 * the pinhole plate.
 * 
 * INPUTS:

 * 
 * OUTPUTS:
 *  dead - int 2, 1, 0, declaring whether the ray is dead. 1 is dead (has not 
 *         met), 0 is alive (has met), 2 is detected (has hit the detector
 *         surface)
 */
int scatterSimpleSurfaces(Ray3D *the_ray, Surface3D *Sample, BackWall Plate, 
        AnalytSphere the_sphere, MTRand *myrng) {
    
    double min_dist;
    int meets;
    int tri_hit;
    double nearest_n[3];
    double nearest_inter[3];
    int meets_sphere;
    int which_surface;
    int detected;
    
    /* tri_hit stores which triangle has been hit */
    tri_hit = -1;
    
    /* By default no detection */
    detected = 0;
    
    /* which_surface stores which surface has been hit */
    which_surface = -1;
    
    /* By default don't hit the sphere */
    meets_sphere = 0;
    
    /* meets is 0/1 have we met a triangle */
    meets = 0;
    
    /* Much further than any of the triangles */
    min_dist = 10.0e10; 

    /* Try to scatter off the sample */
    scatterTriag(the_ray, Sample, &min_dist, nearest_inter, nearest_n, &meets,
        &tri_hit, &which_surface);

    /* Should the sphere be represented */
    if (the_sphere.make_sphere) {
        /* Check the sphere if we are not on it */
        if (the_ray->on_surface != the_sphere.surf_index) {
            meets_sphere = scatterSphere(the_ray, the_sphere, &min_dist, 
                nearest_inter, nearest_n, &tri_hit, &which_surface);
        }
    }

    /* Try to scatter off the simple pinhole plate */
    if (the_ray->on_surface != Plate.surf_index) {
        detected =  backWallScatter(the_ray, Plate, &min_dist, nearest_inter,
            nearest_n, &meets, &tri_hit, &which_surface) ;
    }

    /* If we are detected */
    if (detected) {
        /* Update the ray position but not the direction */
        update_ray_position(the_ray, nearest_inter);
        
        /* 2 = detected ray */
        return 2;
    }
    
    /* Update position/direction etc. */
    if (meets || meets_sphere) {
        double composition;
        double scattering_parameters;
        
        if (meets_sphere) {
            /* sphere is defined to be uniform */
            composition = the_sphere.composition;
            scattering_parameters = the_sphere.scattering_parameters;
        } else {
            if (which_surface == Plate.surf_index) {
                composition = Plate.composition;
                scattering_parameters = Plate.scattering_parameters;
            } else {
                composition = Sample->composition[tri_hit];
                scattering_parameters = Sample->scattering_parameters[tri_hit];
            }
        }
        
        /* Find the new direction and update position*/
        new_direction3D(the_ray, nearest_n, composition, scattering_parameters, myrng);
        update_ray_position(the_ray, nearest_inter);
        
        /* Updates the current triangle and surface the ray is on */
        the_ray->on_element = tri_hit;
        the_ray->on_surface = which_surface;
    } else {
        /* Ray be dead */
        meets = 0;
    }
    
    return(!(meets || meets_sphere));
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
int scatterSimpleMulti(Ray3D *the_ray, Surface3D *Sample, NBackWall Plate, 
        AnalytSphere the_sphere, int *detector, MTRand *myrng) {
    
    double min_dist;
    int meets;
    int tri_hit;
    double nearest_n[3];
    double nearest_inter[3];
    int meets_sphere;
    int which_surface;
    int detected;
    
    /* tri_hit stores which triangle has been hit */
    tri_hit = -1;
    
    /* By default no detection */
    detected = 0;
    
    /* which_surface stores which surface has been hit */
    which_surface = -1;
    
    /* By default don't hit the sphere */
    meets_sphere = 0;
    
    /* meets is 0/1 have we met a triangle */
    meets = 0;
    
    /* Much further than any of the triangles */
    min_dist = 10.0e10; 

    /* Try to scatter off the sample */
    scatterTriag(the_ray, Sample, &min_dist, nearest_inter, nearest_n, &meets,
        &tri_hit, &which_surface);
    
    /* Should the sphere be represented */
    if (the_sphere.make_sphere) {
        /* Check the sphere if we are not on it */
        if (the_ray->on_surface != the_sphere.surf_index) {
            meets_sphere = scatterSphere(the_ray, the_sphere, &min_dist, 
                nearest_inter, nearest_n, &tri_hit, &which_surface);
        }
    }
    
    /* Try to scatter off the simple pinhole plate */
    if (the_ray->on_surface != Plate.surf_index) {
        detected =  multiBackWall(the_ray, Plate, &min_dist, nearest_inter,
            nearest_n, &meets, &tri_hit, &which_surface);
    }

    /* If we are detected */
    if (detected) {
        /* Update the ray position but not the direction */
        update_ray_position(the_ray, nearest_inter);
        
        *detector = detected;
        
        /* 2 = detected ray */
        return 2;
    }
    
    /* Update position/direction etc. */
    if (meets || meets_sphere) {
        double composition;
        double scattering_parameters;
        
        if (meets_sphere) {
            /* sphere is defined to be uniform */
            composition = the_sphere.composition;
            scattering_parameters = the_sphere.scattering_parameters;
        } else {
            if (which_surface == Plate.surf_index) {
                composition = Plate.composition;
                scattering_parameters = Plate.scattering_parameters;
            } else {
                composition = Sample->composition[tri_hit];
                scattering_parameters = Sample->scattering_parameters[tri_hit];
            }
        }
        
        /* Find the new direction and update position*/
        new_direction3D(the_ray, nearest_n, composition, scattering_parameters, myrng);
        update_ray_position(the_ray, nearest_inter);
        
        /* Updates the current triangle and surface the ray is on */
        the_ray->on_element = tri_hit;
        the_ray->on_surface = which_surface;
    } else {
        /* Ray be dead */
        meets = 0;
    }
    
    return(!(meets || meets_sphere));
}


/* 
 * Detects on a hemisphere centered on the sample plane, the aperture is 
 * specified by two angles and a half cone angle.
 */
int scatterAbstractSurfaces(Ray3D *the_ray, Surface3D *Sample, AbstractHemi Plate,
        AnalytSphere the_sphere, MTRand *myrng) {
    
    int meets;
    meets = 0;
    
    return (!meets);
}









