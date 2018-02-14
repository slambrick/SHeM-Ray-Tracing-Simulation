/* 
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 *
 * 
 * Function declerations for the larger ray tracing function for the SHeM
 * Simulation.
 */

#ifndef HEADER_FILE
#define HEADER_FILE

/* 
 * Finds the intersection, normal at the point of intersection and distance to 
 * the intersection between the ray and an analytic sphere. returns 0 if the 
 * ray does not intersect the sphere. 
 */
int scatterSphere(double e[3], double d[3], double normal[3], double *t, 
    double scan_pos_x, scan_pos_z, double dist_to_sample, double sphere_r);

/* 
 * Finds the intersection, normal at the point of intersection and distance to 
 * the intersection between the ray and a triangulated surface. 
 */
void scatterTriag(double e[3], double d[3], int ntriag, double V[], double N[], 
    double F[], int current_tri, int current_surface, double *min_dist, 
    double nearest_inter[3], double nearest_n[3], int *meets, int *tri_hit,
    int this_surface, int *which_surface)

/* 
 * Scatters a ray off a single triangulated surface, and an analytic sphere if
 * desired. 
 */
int scatterOffSurface(double e[3], double d[3], int ntriag, double V[], 
    double N[], double F[], double C[], gsl_rng *myrng, int *current_tri, 
    int *current_surface, double scan_pos_x, double scan_pos_z, int make_sphere,
    double dist_to_sample, double sphere_r, double sphere_difuse);

/* 
 * Scatters a ray off two triangulared surfaces, and an analytic sphere if 
 * desired. 
 */
int scatterSurfaces(double e[3], double d[3], int ntriag_sample, int ntriag_plate, 
    double V[], double N[], double F[], double C[], double VS[], double NS[], 
    double FS[], double CS[], gsl_rng *myrng, int *current_tri, int *current_surface, 
    double backWall[], double scan_pos_x, double scan_pos_z, int make_sphere, 
    double dist_to_sample, double sphere_r, double sphere_diffuse) 

#endif
