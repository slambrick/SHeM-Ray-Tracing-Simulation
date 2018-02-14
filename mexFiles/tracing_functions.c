/* 
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 *
 * 
 * Functions for tracing a single ray off of either one or two triangulated 
 * surfaces. The functions also allow for a single analytically defined sphere
 * to be included as well. Such a feature should be used with care, it is 
 * designed only to be used with a flat sample, though it will work with any 
 * other triangulated sample topology.
 * 
 * These functions make use of the library small_functions and are called in the
 * mex file tracingMex.c.
 */
#include <gsl/gsl_rng.h>
#include "small_functions.h"
#include "tracing_functions.h"
#include <math.h>

/* 
 * Finds the distance to, the normal to, and the position of a rays intersetion
 * with an analytically defined sphere. Returns 0 if the ray does not intersect
 * the sphere. The inclusion of the analytic sphere is not very general and can
 * only really be used for the specific probable it was written for without 
 * alteration--a single sphere places just touching a flat surface in the centre
 * of the scan region.
 * 
 * INPUTS:
 *  e              - double array, the initial position of the ray
 *  d              - double array, the initial direction of the ray
 *  normal         - double array, a variable to store the normal to the point 
 *                   of intersection in
 *  t              - double, a variable to store the distance to the 
 *                   intersection point in
 *  scan_pos_x     - double, the x position of the scan this pixel is
 *  scan_pos_z     - double, the z position of the scan this pixel is
 *  dist_to_sphere - double, the distance between the pinhole plate and the flat
 *                   surface below the sphere
 *  sphere_r       - double, the radius of the sphere
 * 
 * OUTPUTS: 
 *  intersect - int, 1 or 0 depending on if the ray intersects the sphere.
 *  
 */
static int scatterSphere(double e[3], double d[3], double normal[3], double *t, 
        double scan_pos_x, double scan_pos_z, double dist_to_sphere, 
        double sphere_r) {
    double a,b,c;
    double alpha, beta, gamma;
    int meets;
    double distance;
    double offset;
    
    /* Offset of the flat region from specular */
    offset = dist_to_sphere - 2.121;
    
    /* Centre of the sphere */
    a = scan_pos_x + offset;
    b = -2.121 - offset;
    c = scan_pos_z;
    b = b + sphere_r;
    
    /* Coefficients of the quadratic equation */
    alpha = 1; /* Is always 1 because of normailsed direction*/
    beta = 2*(d[0]*(e[0] - a) + d[1]*(e[1] - b) + d[2]*(e[2] - c));
    gamma = -sphere_r*sphere_r - 2*e[0]*a - 2*e[1]*b - 2*e[2]*c + e[0]*e[0] + 
        e[1]*e[1] + e[2]*e[2] + a*a + b*b + c*c;
    
    /* Do we hit the sphere */
    if (beta*beta - 4*gamma < 0) {
        return(0);
    }
    
    /* Solve the quadratic equation. Take the smaller root */
    distance = (-beta - sqrt(beta*beta - 4*gamma))/2;
    
    /* Normal to the sphere at that point */
    normal[0] = (e[0] + distance*d[0] - a)/sphere_r;
    normal[1] = (e[1] + distance*d[1] - b)/sphere_r;
    normal[2] = (e[2] + distance*d[2] - c)/sphere_r;
    
    *t = distance;
    
    normalise(normal);
    
    return(1);
}

/*
 * Finds the distance to, the normal to, and the position of a ray's intersetion
 * with an triangulated surface. 
 * 
 * INPUTS:
 *  e               - double array, the initial position of the ray
 *  d               - double array, the initial direction of the ray
 *  ntriag          - int, the number of triangles in the surface
 *  V               - double array, a 2d array of the coordinates of the 
 *                    vertices in the surface
 *  N               - double array, a 2d array of the normals to all the faces
 *                    in the surface
 *  F               - double array, a 2d array of the indices of which vertices
 *                    make up each triangle
 *  current_tri     - int, the index of the triangle the current ray is on, -1 
 *                    indicates that the ray is not on any triangle
 *  current_surface - int, an index stating which surface the ray is on, 
 *                    0=sample, 1=pinhole plate, -1=none, -2=sphere
 *  min_dist        - double pointer, to store the minium distance to a surface
 *                    in
 *  nearest_inter   - double array, to store the location of the nearest 
 *                    intersection in
 *  nearest_n       - double array, to store the normal to the surface at the 
 *                    nearest instersection
 *  meets           - int, 1 or 0 to store if the ray hits any surfaces
 *  tri_hit         - int pointer, to store the index of the triangle that the 
 *                    ray hits (if it hits any)
 *  this_surface    - int, the index of this triangulated surface
 *  which_surface   - int pointer, which surface does the ray intersect (if any)
 * 
 * NOTE: this function is messy as attempts (mostly successful) have been made 
 *       to improve the speed of the simulation as this is the section of code
 *       called the highest number of times, hence the rather low level looking
 *       code.
 */
static void scatterTriag(double e[3], double d[3], int ntriag, double V[], double N[], 
        double F[], int current_tri, int current_surface, double *min_dist, 
        double nearest_inter[3], double nearest_n[3], int *meets, int *tri_hit,
        int this_surface, int *which_surface) {
    int j;
    double normal[3];
    
    /* Loop through all triangles in the surface */
    for (j = 0; j < ntriag; j++) {
            int tri[3];
        double a[3];
        double b[3];
        double c[3];
        double AA[3][3];
        double v[3];
        double u[3];
        double epsilon;
        int k;
        
        /* Skip this tiangle if the ray is already on it */
        if ((current_tri == j) && (current_surface == this_surface)) {
            continue;
        }
        
        /* 
         * Specify which triangle and get its normal.
         * F is imported as double, we cast to int in tri to use as indecies
         * when retriving the locations of the vertices of the triangles.
         */
        k = j*3 + 0;
        tri[0] = ((int)F[k] - 1)*3;
        normal[0] = N[k];
        k += 1;
        tri[1] = ((int)F[k] - 1)*3;
        normal[1] = N[k];
        k += 1;
        tri[2] = ((int)F[k] - 1)*3;
        normal[2] = N[k];
        
        /* If the triangle is 'back-facing' then the ray cannot hit it */
        if (dot(normal, d) > 0) { 
            continue;
        }
        
        /* Vertices of the triangle */
        a[0] = V[tri[0] + 0];
        a[1] = V[tri[0] + 1];
        a[2] = V[tri[0] + 2];
        b[0] = V[tri[1] + 0];
        b[1] = V[tri[1] + 1];
        b[2] = V[tri[1] + 2];
        c[0] = V[tri[2] + 0];
        c[1] = V[tri[2] + 1];
        c[2] = V[tri[2] + 2];
        
        /* 
         * If the triangle is behind the current ray position then the ray 
         * cannot hit it. To do this we have to test each of the three vertices 
         * to find if they are `behind' the ray. If any one of the vertices is 
         * infornt of the ray we have to consider it. Re-use variable v.
         */
        v[0] = a[0] - e[0];
        v[1] = a[1] - e[1];
        v[2] = a[2] - e[2];
        if (v[0]*d[0] + v[1]*d[1] + v[2]*d[2] < 0) {
            v[0] = b[0] - e[0];
            v[1] = b[1] - e[1];
            v[2] = b[2] - e[2];
            if (v[0]*d[0] + v[1]*d[1] + v[2]*d[2] < 0) {
                v[0] = c[0] - e[0];
                v[1] = c[1] - e[1];
                v[2] = c[2] - e[2];
                if (v[0]*d[0] + v[1]*d[1] + v[2]*d[2] < 0) {
                    continue;
                }
            }
        }
        
        /* 
         * Construct the linear equation
         * AA u = v, where u contains (alpha, beta, t) for the propogation 
         * equation:
         * e + td = a + beta(b - a) + gamma(c - a)
         */
        /*propogate(a, e, -1, v);*/
        v[0] = a[0] - e[0];
        v[1] = a[1] - e[1];
        v[2] = a[2] - e[2];
        
        /* This could be pre-calculated and stored, however it would involve an 
         * array of matrices
         */
        AA[0][0] = a[0] - b[0];
        AA[0][1] = a[0] - c[0];
        AA[0][2] = d[0];
        AA[1][0] = a[1] - b[1];
        AA[1][1] = a[1] - c[1];
        AA[1][2] = d[1];
        AA[2][0] = a[2] - b[2];
        AA[2][1] = a[2] - c[2];
        AA[2][2] = d[2];
        
        /* 
         * Tests to see if this triangle is parrallel to the ray, if it is the
         * determinant of the matrix AA will be zero, we must set a tolerance for 
         * size of determinant we will allow.
         */
        epsilon = 0.00000000000001;
        
        if (!solve3x3(AA, u, v, epsilon)) {
            continue;
        }
        
        /* Find if the point of intersection is inside the triangle */
        /* Must also find if the ray is propogating forwards */
        if ((u[0] >= 0) && (u[1] >= 0) && ((u[0] + u[1]) <= 1) && (u[2] > 0)) {
            double new_loc[3];
            double movment[3];
            double dist;
            
            /* We have hit a triangle */
            *meets = 1;
            
            /* Store the location and normal of the nearest intersection */
            /*propogate(e, d, u[2], new_loc);*/
            new_loc[0] = e[0] + (u[2]*d[0]);
            new_loc[1] = e[1] + (u[2]*d[1]);
            new_loc[2] = e[2] + (u[2]*d[2]);
            
            /* Movment is the vector from the current location to the possible 
             * new location */
            movment[0] = new_loc[0] - e[0];
            movment[1] = new_loc[1] - e[1];
            movment[2] = new_loc[2] - e[2];
            
            /* NOTE: we are comparing the square of the distance */
            dist = movment[0]*movment[0] + movment[1]*movment[1] + 
                movment[2]*movment[2];
            
            if (dist < *min_dist) {
                *min_dist = dist;
                
                *tri_hit = j;
                nearest_n[0] = normal[0];
                nearest_n[1] = normal[1];
                nearest_n[2] = normal[2];
                nearest_inter[0] = new_loc[0];
                nearest_inter[1] = new_loc[1];
                nearest_inter[2] = new_loc[2];
                
                *which_surface = this_surface;
            }
        }
    }
}

/*
 * Scatters the given ray off a single triangulated surface, the sample, and an
 * analytic sphere, if that is desired. Returns true if the ray did not hit the 
 * surface, i.e. it is 'dead' and returns false if the ray does, i.e. the ray
 * is 'alive'. This function has undergone some low level optimisation, so it 
 * may not be written in the most intuitive and simple manner.
 * 
 * INPUTS:
 *  e               - double array, position of the ray, if the ray is scattered
 *                    this is updated
 *  d               - double array, direction of the ray, if the ray is 
 *                    scattered this is updated
 *  ntriag          - int, the number of triangles in the sample surface
 *  V               - double array, list of the coordinates of the vertices in 
 *                    the triangulation of the surface
 *  N               - double array, list of normals to the triangles in the 
 *                    surface
 *  F               - double array, list containing the reference to the faces 
 *                    of the triangles
 *  C               - double, array, list of the type of scattering off of each
 *                    triangle in the surface
 *  myrng           - random number generator object to be used
 *  current_tri     - int, the index of the triangle the current ray is on, -1 
 *                    indicates that the ray is not on any triangle
 *  current_surface - int, an index stating which surface the ray is on, 
 *                    0=sample, 1=pinhole plate, -1=none, -2=sphere
 *  scan_pos_x      - double, the x position of the scan this pixel is
 *  scan_pos_z      - double, the z position of the scan this pixel is
 *  make_sphere     - int, 1 or 0, shoule an analytic sphere be included
 *  dist_to_sample  - double, the distance between the pinhole plate and the
 *                    flat surface below the sphere
 *  sphere_r        - double, radius of the analytic sphere
 *  sphere_diffuse  - double, the type of scattering off of the analytic sphere
 * 
 * OUTPUTS:
 *  dead - int 1, 0 declaring whether the ray is 'dead', 1 is dead (has not 
 *         met), 0 is alive (has met)
 */
int scatterOffSurface(double e[3], double d[3], int ntriag, double V[], 
        double N[], double F[], double C[], gsl_rng *myrng, int *current_tri, 
        int *current_surface, double scan_pos_x, double scan_pos_z, 
        int make_sphere, double dist_to_sample, double sphere_r, 
        double sphere_diffuse) {
    
    double min_dist;
    int meets;
    int j, tri_hit;
    double nearest_n[3];
    double nearest_inter[3];
    double new_dir[3];
    double t;
    double normal[3];
    int meets_sphere;
    int which_surface;
    
    tri_hit = -1;
    meets = 0;
    
    /* Much further than any of the triangles */
    min_dist = 10.0e10;
    
    /* Try to scatter of the sample */
    scatterTriag(e, d, ntriag, V, N, F, *current_tri, *current_surface, 
        &min_dist, nearest_inter, nearest_n, &meets, &tri_hit, 0, &which_surface);
    
    /* Should the sphere be represented */
    if (make_sphere == 1) {
        /* Check the sphere if we are not on it */
        if (*current_surface != -2) {
            meets_sphere = scatterSphere(e, d, normal, &t, scan_pos_x, scan_pos_z, dist_to_sample, sphere_r);
        }
        
        if (meets_sphere && (t > 0)) {
            /* We have hit a surface, the sphere */
            meets = 1;
            
            /* Is the sphere the closest surface */
            if (t < min_dist) {
                min_dist = t*t;
                
                /* Update the position */
                nearest_inter[0] = e[0] + d[0]*t;
                nearest_inter[1] = e[1] + d[1]*t;
                nearest_inter[2] = e[2] + d[2]*t;
                
                /* Update the surfac normal */
                nearest_n[0] = normal[0];
                nearest_n[1] = normal[1];
                nearest_n[2] = normal[2];
                
                which_surface = -2;
                tri_hit = -1;
            } else {
                meets_sphere = 0;
            }
        }
    }
    
    /* If we have met a triangle/sphere we must scatter off of it */
    if (meets) {
        double tester;
        double diffuseLvl;
        
        /* Update the position */
        for (j = 0; j < 3; j++) {
            e[j] = nearest_inter[j];
        }
        
        if (meets_sphere == 1) {
            /* sphere is defined to be uniform */
            diffuseLvl = sphere_diffuse;
            tri_hit = -2;
        } else {
            diffuseLvl = C[tri_hit];
        }
        
        /* Find the new direction */
        /* If diffuseLvl is 2 then use a uniform distribution, if it is [0,1] 
         * then use either a cosine or specular distribution, or a combination 
         * of both
         */
        if (fabs(diffuseLvl - 2) < 0.000001) {
            uniformScatter(nearest_n, new_dir, myrng);
        } else {
            tester = gsl_rng_uniform(myrng);
            if (tester < diffuseLvl) {
                cosineScatter(nearest_n, new_dir, myrng);
            } else {
                reflect(nearest_n, d, new_dir);
            }
        }
        
        /* Update the direction */
        for (j = 0; j < 3; j++) {
            d[j] = new_dir[j];
        }
        
        /* Updates the current triangle and surface the ray is on */
        *current_tri = tri_hit;
        *current_surface = which_surface;
    }
    
    return (!meets);
}

/*
 * INPUTS:
 *  e               - double array, position of the ray, if the ray is scattered
 *                    this is updated
 *  d               - double array, direction of the ray, if the ray is 
 *                    scattered this is updated
 *  ntriag_sample   - int, the number of triangles in the sample surface
 *  ntriag_plate    - int, the number of triangles in the pinhole plate
 *  V               - double array, list of the coordinates of the vertices in 
 *                    the triangulation of the sample surface
 *  N               - double array, list of normals to the triangles in the 
 *                    sample surface
 *  F               - double array, list containing the reference to the faces 
 *                    of the triangles in the sample surface
 *  C               - double, array, list of the type of scattering off of each
 *                    triangle in the sample surface
 *  VS              - double array, list of the coordinates of the vertices in 
 *                    the triangulation of the pinhole plate
 *  NS              - double array, list of normals to the triangles in the 
 *                    pinhole plate
 *  FS              - double array, list containing the reference to the faces 
 *                    of the triangles in the pinhole plate
 *  CS              - double, array, list of the type of scattering off of each
 *                    triangle in the pinhole plate
 *  myrng           - random number generator object to be used
 *  current_tri     - int, the index of the triangle the current ray is on, -1 
 *                    indicates that the ray is not on any triangle
 *  current_surface - int, an index stating which surface the ray is on, 
 *                    0=sample, 1=pinhole plate, -1=none, -2=sphere
 *  backWall        - double array, contains information about the location of 
 *                    the detector surface, first element gives the y 
 *                    coordinate of the back of the pinhole plate, second 
 *                    element gives the depth of the pinhole plate in x, and the
 *                    third element gives the dpeth of the pinhole plate in z
 *                    gives the 
 *  scan_pos_x      - double, the x position of the scan this pixel is
 *  scan_pos_z      - double, the z position of the scan this pixel is
 *  make_sphere     - int, 1 or 0, shoule an analytic sphere be included
 *  dist_to_sample  - double, the distance between the pinhole plate and the
 *                    flat surface below the sphere
 *  sphere_r        - double, radius of the analytic sphere
 *  sphere_diffuse  - double, the type of scattering off of the analytic sphere
 * 
 * OUTPUTS:
 *  dead - int 2, 1, 0, declaring whether the ray is dead. 1 is dead (has not 
 *         met), 0 is alive (has met), 2 is detected (has hit the detector
 *         surface)
 */
int scatterSurfaces(double e[3], double d[3], int ntriag_sample, int ntriag_plate, 
        double V[], double N[], double F[], double C[], double VS[], double NS[], 
        double FS[], double CS[], gsl_rng *myrng, int *current_tri, int *current_surface, 
        double backWall[], double scan_pos_x, double scan_pos_z, int make_sphere, 
        double dist_to_sample, double sphere_r, double sphere_diffuse) {
    
    double min_dist;
    int meets;
    int j, tri_hit;
    double nearest_n[3];
    double nearest_inter[3];
    double new_dir[3];
    double normal[3];
    double t;
    int meets_sphere;
    int which_surface;
    
    /* tri_hit stores which triangle has been hit */
    tri_hit = -1;
    
    /* which_surface stores which surface has been hit */
    which_surface = -1;
    
    /* meets is 0/1 have we met a triangle */
    meets = 0;
    
    /* Much further than any of the triangles */
    min_dist = 10.0e10; 
    
    /* Try to scatter off the sample */
    scatterTriag(e, d, ntriag_sample, V, N, F, *current_tri, *current_surface, 
        &min_dist, nearest_inter, nearest_n, &meets, &tri_hit, 0, &which_surface);
    
    /* Try to scater off the pinhole plate */
    scatterTriag(e, d, ntriag_plate, VS, NS, FS, *current_tri, *current_surface, 
        &min_dist, nearest_inter, nearest_n, &meets, &tri_hit, 1, &which_surface);
    
    /* Represent the sphere? */
    if (make_sphere == 1) {
        /* Check the sphere if we are not on it */
        if (*current_surface != -2) {
            meets_sphere = scatterSphere(e, d, normal, &t, scan_pos_x, 
                                         scan_pos_z, dist_to_sample, sphere_r);
        }
        
        if (meets_sphere == 1) {
            /* We have hit a surface, the sphere */
            meets = 1;
            
            /* Is the sphere the closest surface */
            if (t < min_dist) {
                min_dist = t;
                
                /* Update the position */
                nearest_inter[0] = e[0] + d[0]*t;
                nearest_inter[1] = e[1] + d[1]*t;
                nearest_inter[2] = e[2] + d[2]*t;
                
                /* Update the surface normal */
                nearest_n[0] = normal[0];
                nearest_n[1] = normal[1];
                nearest_n[2] = normal[2];
                
                
                which_surface = -2;
                tri_hit = -2;
            } else {
                meets_sphere = 0;
            }
        }
    }
    
    /* Update position/direction etc. */
    if (meets) {
        double tester;
        double diffuseLvl;
        
        /* Update the position */
        for (j = 0; j < 3; j++) {
            e[j] = nearest_inter[j];
        }
        
        /* Find the new direction */
        tester = gsl_rng_uniform(myrng);
        
        if (meets_sphere == 1) {
            /* sphere is defined to be uniform */
            diffuseLvl = sphere_diffuse;
        } else {
            if (which_surface == 1) {
                /* We hit the pinhole plate */
                diffuseLvl = CS[tri_hit];
            } else {
                /* We hit the sample */
                diffuseLvl = C[tri_hit];
            }
        }
        
        /* If diffuseLvl is 2 then use a uniform distribution, if it is [0,1] 
         * then use either a cosine or specular distribution, or a combination 
         * of both
         */
        if (fabs(diffuseLvl - 2) < 0.000001) {
            uniformScatter(nearest_n, new_dir, myrng);
        } else {
            if (tester < diffuseLvl) {
                cosineScatter(nearest_n, new_dir, myrng);
            } else {
                reflect(nearest_n, d, new_dir);
            }
        }
        
        /* Update the direction */
        for (j = 0; j < 3; j++) {
            d[j] = new_dir[j];
        }
        
        /* Updates the current triangle the ray is on */
        *current_tri = tri_hit;
        
        /* Update the current surface the ray is on */
        *current_surface = which_surface;
    } else {
        /* 
         * We must consider if the ray has been detected if it hasn't hit 
         * either surface 
         */
        
        /* First consider if the ray is propogating in the +ve y direction */
        if (d[1] > 0) {
            double alpha;
            double wall_hit[3];
            
            /* Find where, just behind the pinhole plate, the ray will hit, 
             * backWall[0] is the y coordinate of the back of the pinhole plate
             */
            alpha = (backWall[0] - e[1])/d[1];
            
            propogate(e, d, alpha, wall_hit);
            
            /* Now find if this point is covered by the plate, if it is then the
             * ray must be detected, if it is not, then the ray must be dead.
             * backWall[1] and backWall[2] are the depth in x and z of the 
             * pinhole plate
             */
            if ((fabs(wall_hit[0]) < (backWall[1]/2)) && (fabs(wall_hit[2]) < 
                    (backWall[2]/2))) {
                /* We have detection! */
                
                for (j = 0; j < 3; j++) {
                    /* We update position only and keep the direction the same */
                    e[j] = wall_hit[j];
                }
                
                return 2;
            }
            
        } else {
            /* Ray be dead */
            meets = 0;
        }
    }
    
    return (!meets);
}
