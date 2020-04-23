/*
 * Copyright (c) 2018-19, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 *
 * 
 * Functions for setting up the structures for the SHeM Ray Tracing Simulation,
 * they take arguments directly aquired from MATLAB and put them into structures
 * for use in the C ray tracing code, also contains functions for cleaning up
 * at the end of the simulation.
 */
#include "mex.h"
#include "common_helpers.h"
#include "small_functions3D.h"
#include "ray_tracing_structs3D.h"
#include "scattering_processes3D.h"
#include <stdlib.h>
//#include <gsl/gsl_rng.h>
#include <math.h>

/* 
 * Set up a surface containing the information on a triangulated surface.
 * 
 * INPUTS: 
 *  V      - 3xm array of vertex positions
 *  N      - 3xn array of which vertices form faces of triangles
 *  N      - 3xn array of the unit normals to the surface elements
 *  C      - 1xn array of the composition of the surface elements
 *  ntriag - The number of triangles in the surface
 * 
 * OUTPUT:
 *  Sample - a Surface struct that contains information of the surface.
 */
Surface3D set_up_surface(double V[], double N[], double F[], double C[], 
        double P[], int ntraig, int surf_index) {
    Surface3D Surf;
    
    /* Allocate the components of the surface. */
    Surf.surf_index = surf_index;
    Surf.n_elements = ntraig;
    Surf.vertices = V;
    Surf.normals = N;
    Surf.faces = F;
    Surf.composition = C;
    Surf.scattering_parameters = P;
    
    /* Return the struct itslef */
    return(Surf);
}

/* Set up a Sphere struct */
AnalytSphere set_up_sphere(int make_sphere, double *sphere_c, double sphere_r,
        double sphere_scattering, double sphere_parameters, int surf_index) {
    AnalytSphere the_sphere;
    
    the_sphere.surf_index = surf_index;
    the_sphere.sphere_c = sphere_c;
    the_sphere.sphere_r = sphere_r;
    the_sphere.composition = sphere_scattering;
    the_sphere.scattering_parameters = sphere_parameters;
    the_sphere.make_sphere = make_sphere;
    
    return(the_sphere);
}

/* 
 * Set up a struct of rays using the input vectors from MATLAB.
 * 
 */
Rays3D compose_rays3D(double ray_pos[], double ray_dir[], int nrays) {
    int i;
    Rays3D all_rays;
    Ray3D *rays;
    
    /* Allocated the memory to the  */
    rays = (Ray3D*)malloc(nrays * sizeof(*rays));
    
    /* Put the rays defined by the arrays into an array of structs rays */
    for (i = 0; i < nrays; i++) {
        int n;
        
        n = i*3;
        rays[i].position[0] = ray_pos[n];
        rays[i].position[1] = ray_pos[n+1];
        rays[i].position[2] = ray_pos[n+2];
        rays[i].direction[0] = ray_dir[n];
        rays[i].direction[1] = ray_dir[n+1];
        rays[i].direction[2] = ray_dir[n+2];
        
        /* Initialise ray not to be on any surface element or any surface */
        rays[i].on_element = -1;
        rays[i].on_surface = -1;
        rays[i].nScatters = 0;
    }
    
    /* Put the data into the struct */
    all_rays.rays = rays;
    all_rays.nrays = nrays;
    
    return(all_rays);
}

/* Updates the ray postion */
void update_ray_position(Ray3D *the_ray, double new_pos[3]) {
    int j;
    
    for (j = 0; j < 3; j++)
        the_ray->position[j] = new_pos[j];
}

/* Updates the ray direction */
void update_ray_direction(Ray3D *the_ray, double new_dir[3]) {
    int j;
    
    for (j = 0; j < 3; j++)
        the_ray->direction[j] = new_dir[j];
}

/* 
 * Gets the nth element of a Surface struct and puts the information in the 
 * provided arrays.
 * 
 * INPUTS:
 *  Sample      - a pointer to a Surface struct with all the information on the
 *               surface in it
 *  n           - integer, the index of the surface element
 *  v1          - three element array for putting the first vertex in
 *  v2          - three element array for putting the second vertex in
 *  v3          - three element array for putting the third vertex in
 *  nn          - three element array for putting the unit normal in
 *  composition - pointer to a variable for puttin the element composition in 
 */
void get_element3D(Surface3D *Sample, int n, double v1[3], double v2[3], 
        double v3[3], double nn[3]) {
    int j;
    int tri[3];
    
    /* Get the indices to the three vertices and get the surface normal */
    j = n*3 + 0;
    tri[0] = ((int)Sample->faces[j] - 1)*3;
    nn[0] = Sample->normals[j];
    j += 1;
    tri[1] = ((int)Sample->faces[j] - 1)*3;
    nn[1] = Sample->normals[j];
    j += 1;
    tri[2] = ((int)Sample->faces[j] - 1)*3;
    nn[2] = Sample->normals[j];
    
    /* Vertices of the triangle */
    v1[0] = Sample->vertices[tri[0] + 0];
    v1[1] = Sample->vertices[tri[0] + 1];
    v1[2] = Sample->vertices[tri[0] + 2];
    v2[0] = Sample->vertices[tri[1] + 0];
    v2[1] = Sample->vertices[tri[1] + 1];
    v2[2] = Sample->vertices[tri[1] + 2];
    v3[0] = Sample->vertices[tri[2] + 0];
    v3[1] = Sample->vertices[tri[2] + 1];
    v3[2] = Sample->vertices[tri[2] + 2];
}


/* Cleanup a struct of rays */
void clean_up_rays(Rays3D all_rays) {
    free(all_rays.rays);
}

/* Get the ray positions and put them in an output array */
void get_positions(Rays3D *all_rays, double *final_pos) {
    int i;
    
    /* Loop through all the rays */
    for (i = 0; i < all_rays->nrays; i++) {
        int k;
        Ray3D *current_ray;
        
        current_ray = &all_rays->rays[i];
        for (k = 0; k < 3; k++) {
            int n;
            n = k + 3*i;
            final_pos[n] = current_ray->position[k];
        }
    }
}

/* Get the ray directions and put them in an output array */
void get_directions(Rays3D *all_rays, double *final_dir) {
    int i;
    
    /* Loop through all the rays */
    for (i = 0; i < all_rays->nrays; i++) {
        int k;
        Ray3D *current_ray;
        
        current_ray = &all_rays->rays[i];
        for (k = 0; k < 3; k++) {
            int n;
            n = k + 3*i;
            final_dir[n] = current_ray->direction[k];
        }
    }
}

/* Get the number of scattering events per ray */
void get_scatters(Rays3D *all_rays, int *nScatters) {
    int i;
    
    /* Loop through all the rays */
    for (i = 0; i < all_rays->nrays; i++) {
        Ray3D *current_ray;
        
        current_ray = &all_rays->rays[i];
        nScatters[i] = (int)current_ray->nScatters;
    }
}

/* Prints all the information about the ray to the terminal */
void print_ray(Ray3D *the_ray) {
    mexPrintf("\non_element = %i\n", the_ray->on_element);
    mexPrintf("on_surface = %i\n", the_ray->on_surface);
    mexPrintf("nScatters = %i\n", the_ray->nScatters);
    mexPrintf("Position: ");
    print1D_double(the_ray->position, 3);
    mexPrintf("Direction: ");
    print1D_double(the_ray->direction, 3);
}

/* Prints all the information about a BackWall struct */
void print_BackWall(BackWall *wall) {
    mexPrintf("\nSurface index = %i\n", wall->surf_index);
    mexPrintf("Plate represent = %i\n", wall->plate_represent);
    mexPrintf("Aperture centre: ");
    print1D_double(wall->aperture_c, 2);
    mexPrintf("Aperture axes: "); 
    print1D_double(wall->aperture_axes, 2);
    mexPrintf("Radius of the plate = %f\n", wall->circle_plate_r);
    mexPrintf("Composition = %f\n", wall->composition);
    mexPrintf("Scattering parameters = %f\n", wall->scattering_parameters);
}

/* Prints all the information on all the apertues in the NBackWall struct */
void print_nBackWall(NBackWall *all_apertures) {
    int i;
    
    mexPrintf("\nNumber of apertures = %i\n", all_apertures->n_detect);
    for (i = 0; i < all_apertures->n_detect; i++) {
        mexPrintf("Aperture %i:\n", i);
        mexPrintf("Centre: ");
        print1D_double(&all_apertures->aperture_c[2*i], 2);
        mexPrintf("Axes: ");
        print1D_double(&all_apertures->aperture_axes[2*i], 2);
    }
}

/* 
 * Gets the centre and the axes of the nth detector in the series and returns
 * a struct to that aperture.
 */
void get_nth_aperture(int n, NBackWall *allApertures, BackWall *this_wall) {
    /* The x and z coordinate of the apertures */
    this_wall->aperture_c = &allApertures->aperture_c[2*n];
        
    /* The axes of the aperture */
    this_wall->aperture_axes = &allApertures->aperture_axes[2*n];
    
    /* Other parameters */
    this_wall->surf_index = allApertures->surf_index;
    this_wall->circle_plate_r = allApertures->circle_plate_r;
    this_wall->composition = allApertures->composition;
    this_wall->scattering_parameters = allApertures->scattering_parameters;
    
    /* We do not represent the plate to be scattered off, just the apertures */
    this_wall->plate_represent = 0;
}

/* 
 * Creates a ray according to a simple model of the source, can use either a uniform
 * or Gaussian virtual source.
 * 
 * INPUTS:
 *  pinhole_r    - 
 *  pinhole_c    - 
 *  theta_max    -  
 *  init_angle   - 
 *  source_model -
 *  my_rng       - 
 *  sigma        -  
 * 
 * OUTPUT:
 *  gen_ray - ray3D struct with information on a ray in it
 */
Ray3D create_ray_source(double pinhole_r, double *pinhole_c, double theta_max, 
        double init_angle, int source_model, double sigma) {
    Ray3D gen_ray;
    double r, theta, phi;
    double rot_angle;
    double B;
    double normal[3];
    double dir[3];
    double rand1, rand2, rand3;
    
    /* Enuse theta is initialized */
    theta = 0;
    
    /* Generate the position of the ray */
    rand1 = (double)rand() / (double)RAND_MAX;
    rand2 = (double)rand() / (double)RAND_MAX;
    phi = 2*M_PI*rand1;
    r = pinhole_r*sqrt(rand2);
    gen_ray.position[0] = r*cos(phi);
    gen_ray.position[1] = 0;
    gen_ray.position[2] = r*sin(phi);
    
    /* Generate the direction of the ray */
    rand3 = (double)rand() / (double)RAND_MAX;
    phi = 2*M_PI*rand3;
    switch (source_model) {
        case 0:
            /* Uniform virtual source model */
            rand1 = (double)rand() / (double)RAND_MAX;
            theta = theta_max*sqrt(rand1);
            break;
        case 1:
            /* Gaussian virtual source model */
            B = 1/(1 - exp(-M_PI*M_PI/(2*sigma*sigma)));
            rand1 = (double)rand() / (double)RAND_MAX;
            theta = sigma*sqrt(-2*log((B - rand1/B)));
            break;
        case 2:
            /* Diffuse cosine model */
            normal[0] = 0;
            normal[1] = -1;
            normal[2] = 0;
            cosineScatter3D(normal, gen_ray.direction);
            break;
    }
    
    if (source_model != 2) {
        dir[0] = cos(theta);
        dir[1] = sin(theta)*cos(phi);
        dir[2] = sin(theta)*sin(phi);
        
        /* Need to rotate the rays direction according to the incidence angle */
        rot_angle = 0.5*M_PI - init_angle;
        gen_ray.direction[0] = cos(-rot_angle)*dir[0] - 
            sin(-rot_angle)*dir[1];
        gen_ray.direction[1] = sin(-rot_angle)*dir[0] + 
            cos(-rot_angle)*dir[1];
        gen_ray.direction[2] = dir[2];
    }
    
    /* Need to move the rays into the pinhole */
    gen_ray.position[0] += pinhole_c[0];
    gen_ray.position[2] += pinhole_c[2];
    
    /* Initialise other elements of the ray struct */
    gen_ray.on_surface = -1;
    gen_ray.on_element = -1;
    gen_ray.nScatters = 0;
    
    /* Return the ray struct all ready to use */
    return(gen_ray);
}

