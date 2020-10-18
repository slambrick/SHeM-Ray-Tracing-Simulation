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
#include "ray_tracing_structs3D.h"
#include "common_helpers.h"
#include "small_functions3D.h"
#include "distributions3D.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mex.h>
#include "mtwister.h"

/*
 * Set up a surface containing the information on a triangulated surface.
 * At the end of the program, MUST call clean_up_surface to free allocated memory.
 *
 * INPUTS:
 *  V      - 3xm array of vertex positions
 *  N      - 3xn array of which vertices form faces of triangles
 *  N      - 3xn array of the unit normals to the surface elements
 *  C      - 1xn array of strings pointing to material names
 *  M      - array of Material objects holding material info
 *  ntriag - The number of triangles in the surface
 *  nmaterials - number of materials in M
 *  surf_index - surface index for scattering algorithm
 *
 * OUTPUT:
 *  surf - a Surface struct that contains information of the surface.
 */
void set_up_surface(double V[], double N[], int F[], char ** C, Material M[],
		int nmaterials, int ntriag, int nvert, int surf_index, Surface3D * surf) {

    /* Allocate the components of the surface. */
    surf->surf_index = surf_index;
    surf->n_faces = ntriag;
    surf->n_vertices = nvert;
    surf->vertices = V;
    surf->normals = N;
    surf->faces = F;

    // assign references to the correct material
    // loop through faces and look for the material that fits the name
    surf->compositions = mxCalloc(ntriag, sizeof(Material*));
    
    for (int iface = 0; iface < ntriag; iface++) {
        bool found = false;
        int imat = 0;
        do {
            if(strcmp(C[iface], M[imat].name) == 0) {
                found = true;
                surf->compositions[iface] = &M[imat];
            }
            imat++;
        } while(!found && imat < nmaterials);
        if(!found)
            mexErrMsgIdAndTxt("MyToolbox:tracingMex:compositions",
                              "Composition of face %d not resolved.", iface);
    }
}

void clean_up_surface(Surface3D * surface) {
    mxFree(surface->compositions);
}

/* Set up a Sphere struct */
void set_up_sphere(int make_sphere, double *sphere_c, double sphere_r,
                           Material material, int surf_index, AnalytSphere * sph) {
    sph->surf_index = surf_index;
    sph->sphere_c = sphere_c;
    sph->sphere_r = sphere_r;
    sph->make_sphere = make_sphere;
    sph->material = material;
}


/*
 * Initialise a Material with given props. The function name will also
 * be resolved, i.e. the scattering distributions will be searched by name
 * and the distribution function will be assigned to the func field.
 */
void set_up_material(char * name, char * function, double * params,
		int n_params, Material *mat) {
    mat->name = name;
    mat->func_name = function;
    mat->params = params;
    mat->n_params = n_params;

    mat->func = distribution_by_name(mat->func_name);
    if(mat->func == NULL)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:material",
                          "Distribution name %s could not be resolved.", mat->func_name);
}


/*
 * Set up a struct of rays using the input vectors from MATLAB.
 *
 */
void compose_rays3D(double ray_pos[], double ray_dir[], int nrays, Rays3D *all_rays) {
    int i;
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
    all_rays->rays = rays;
    all_rays->nrays = nrays;
}

/* Updates the ray postion */
void update_ray_position(Ray3D *the_ray, const double new_pos[3]) {
    for (int j = 0; j < 3; j++)
        the_ray->position[j] = new_pos[j];
}

/* Updates the ray direction */
void update_ray_direction(Ray3D *the_ray, const double new_dir[3]) {
    for (int j = 0; j < 3; j++)
        the_ray->direction[j] = new_dir[j];
}

/*
 * Gets the nth element of a Surface struct and puts the information in the
 * provided arrays.
 *
 * INPUTS:
 *  sample      - a pointer to a Surface struct with all the information on the
 *               surface in it
 *  idx         - integer, the index of the surface element
 *  v1          - three element array for putting the first vertex in
 *  v2          - three element array for putting the second vertex in
 *  v3          - three element array for putting the third vertex in
 *  normal      - three element array for putting the unit normal in
 *  composition - pointer to a variable for puttin the element composition in
 */
void get_element3D(const Surface3D *sample, int idx, double v1[3], double v2[3],
        double v3[3], double normal[3]) {
    int j;
    int vertices[3];

    /* Get the indices to the three vertices and get the surface normal */
    j = idx*3 + 0;
    vertices[0] = (sample->faces[j] - 1)*3;
    normal[0] = sample->normals[j];
    j += 1;
    vertices[1] = (sample->faces[j] - 1)*3;
    normal[1] = sample->normals[j];
    j += 1;
    vertices[2] = (sample->faces[j] - 1)*3;
    normal[2] = sample->normals[j];

    /* Vertices of the triangle */
    v1[0] = sample->vertices[vertices[0]];
    v1[1] = sample->vertices[vertices[0] + 1];
    v1[2] = sample->vertices[vertices[0] + 2];
    v2[0] = sample->vertices[vertices[1]];
    v2[1] = sample->vertices[vertices[1] + 1];
    v2[2] = sample->vertices[vertices[1] + 2];
    v3[0] = sample->vertices[vertices[2]];
    v3[1] = sample->vertices[vertices[2] + 1];
    v3[2] = sample->vertices[vertices[2] + 2];
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

/* print details of Material struct */
void print_material(Material const * mat) {
    mexPrintf("\tMAT %-10s func %-12s", mat->name, mat->func_name);
    for (int i = 0; i < mat->n_params; i++)
        mexPrintf(" %.2f ", mat->params[i]);
}

/* Print details of whole sample to console */
void print_surface(Surface3D const * s) {
    for (int ivert = 0; ivert < s->n_vertices; ivert++) {
        int ind0, ind1, ind2;
        lin(ivert, 0, &ind0);
        lin(ivert, 1, &ind1);
        lin(ivert, 2, &ind2);
        mexPrintf("\n VERT % .2f % .2f % .2f", s->vertices[ind0],
            s->vertices[ind1], s->vertices[ind2]);
    }
    for (int iface = 0; iface < s->n_faces; iface++) {
        mexPrintf("\n FACE %2d", iface);
        int ind0, ind1, ind2;
        lin(iface, 0, &ind0);
        lin(iface, 1, &ind1);
        lin(iface, 2, &ind2);
        mexPrintf("\tV %3d %3d %3d", s->faces[ind1],
                    s->faces[ind1], s->faces[ind2]);
        mexPrintf("\tN % .2f % .2f % .2f", s->normals[ind0],
                    s->normals[ind1], s->normals[ind2]);
        print_material(s->compositions[iface]);
    }
    mexPrintf("\n");
}

/* Prints all the information about the ray to the terminal */
void print_ray(Ray3D const * the_ray) {
    mexPrintf("\non_element = %i\n", the_ray->on_element);
    mexPrintf("on_surface = %i\n", the_ray->on_surface);
    mexPrintf("nScatters = %i\n", the_ray->nScatters);
    mexPrintf("Position: ");
    print1D_double(the_ray->position, 3);
    mexPrintf("Direction: ");
    print1D_double(the_ray->direction, 3);
}

/* Prints all the information about a BackWall struct */
void print_BackWall(BackWall const * wall) {
    mexPrintf("\nSurface index = %i\n", wall->surf_index);
    mexPrintf("Plate represent = %i\n", wall->plate_represent);
    mexPrintf("Aperture centre: ");
    print1D_double(wall->aperture_c, 2);
    mexPrintf("Aperture axes: ");
    print1D_double(wall->aperture_axes, 2);
    mexPrintf("Radius of the plate = %f\n", wall->circle_plate_r);
    print_material(&(wall->material));
    mexPrintf("\n");
}

/* Prints all the information on all the apertues in the NBackWall struct */
void print_nBackWall(NBackWall const * all_apertures) {
    mexPrintf("\nNumber of apertures = %i\n", all_apertures->n_detect);
    for (int i = 0; i < all_apertures->n_detect; i++) {
        mexPrintf("Aperture %i:\n", i);
        mexPrintf("Centre: ");
        print1D_double(&all_apertures->aperture_c[2*i], 2);
        mexPrintf("Axes: ");
        print1D_double(&all_apertures->aperture_axes[2*i], 2);
    }
    // print_material(&(all_apertures->material));
    mexPrintf("\n");
}


/* Print the position, radius, material etc of a sphere */
void print_sphere(AnalytSphere const * sphere){
    mexPrintf("\n\t Sphere make %d \t R %3.3f\t C %3.3f %3.3f %3.3f", sphere->make_sphere,
              sphere->sphere_r, sphere->sphere_c[0],
              sphere->sphere_c[1], sphere->sphere_c[2]);
    print_material(&(sphere->material));
}

/*
 * Gets the centre and the axes of the nth detector in the series and returns
 * a struct to that aperture.
 */
void get_nth_aperture(int n, const NBackWall *allApertures, BackWall *this_wall) {
    /* The x and z coordinate of the apertures */
    this_wall->aperture_c = &(allApertures->aperture_c[2*n]);

    /* The axes of the aperture */
    this_wall->aperture_axes = &(allApertures->aperture_axes[2*n]);

    /* Other parameters */
    this_wall->surf_index = allApertures->surf_index;
    this_wall->circle_plate_r = allApertures->circle_plate_r;
    this_wall->material = allApertures->material;

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
void create_ray(Ray3D * gen_ray, SourceParam const * source, MTRand *myrng) {
    double r, theta=0, phi;
    double rot_angle;
    double normal[3];
    double dir[3];
    
    /* Generate the position of the ray */
    double tmp;
    genRand(myrng, &tmp);
    phi = 2*M_PI*tmp;
    genRand(myrng, &tmp);
    r = source->pinhole_r*sqrt(tmp);
    gen_ray->position[0] = source->pinhole_c[0] + r*cos(phi);
    gen_ray->position[1] = source->pinhole_c[1];
    gen_ray->position[2] = source->pinhole_c[2] + r*sin(phi);

    /* Generate the direction of the ray */
    genRand(myrng, &tmp);
    phi = 2*M_PI*tmp;
    switch (source->source_model) {
        case 0:
            /* Uniform virtual source model */
            genRand(myrng, &tmp);
            theta = source->theta_max*sqrt(tmp);
            break;
        case 1:
            /* Gaussian virtual source model */
            genRand(myrng, &tmp);
            theta = source->sigma*sqrt(-2*log((1 - tmp/1)));
            break;
        case 2:
            /* Diffuse cosine model */
            normal[0] = 0;
            normal[1] = -1;
            normal[2] = 0;
            cosine_scatter(normal, NULL, gen_ray->direction, NULL, myrng);
            break;
    }

    if (source->source_model != 2) {
        dir[0] = cos(theta);
        dir[1] = sin(theta)*cos(phi);
        dir[2] = sin(theta)*sin(phi);

        /* Need to rotate the rays direction according to the incidence angle */
        rot_angle = M_PI_2 - source->init_angle;
        gen_ray->direction[0] = cos(-rot_angle)*dir[0] -
            sin(-rot_angle)*dir[1];
        gen_ray->direction[1] = sin(-rot_angle)*dir[0] +
            cos(-rot_angle)*dir[1];
        gen_ray->direction[2] = dir[2];
    }

    /* Initialise other elements of the ray struct */
    gen_ray->on_surface = -1;
    gen_ray->on_element = -1;
    gen_ray->nScatters = 0;
}
