/*
 * Copyright (c) 2018-20, Sam Lambrick.
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
#include "ray_tracing_core3D.h"
#include "common_helpers.h"
#include "distributions3D.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mtwister.h"
#include <stdbool.h>

static void get_normal_ptr(Surface3D const * const s, int ind, double ** n);
static void get_lattice_ptr(Surface3D const * const s, int ind, double ** b);
static void get_face_ptr(Surface3D const * const s, int ind, int32_t ** f);
static void get_vertex_ptr(Surface3D const * const s, int ind, double ** v);

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
void set_up_surface(double V[], double N[], double B[], int32_t F[], char * C[], Material M[],
        int nmaterials, int ntriag, int nvert, int surf_index, Surface3D * const surf) {

    /* Allocate the components of the surface. */
    surf->surf_index = surf_index;
    surf->n_faces = ntriag;
    surf->n_vertices = nvert;
    surf->vertices = V;
    surf->normals = N;
    surf->lattice = B;
    surf->faces = F;

    // assign references to the correct material
    // loop through faces and look for the material that fits the name
    surf->compositions = calloc(ntriag, sizeof(Material*));
    
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
        /*if(!found) {
            mexErrMsgIdAndTxt("MyToolbox:tracingMex:compositions",
                              "Composition of face %d not resolved.", iface);
        }*/
    }
}

void clean_up_surface(Surface3D * const surface) {
    free(surface->compositions);
}

void clean_up_surface_all_arrays(Surface3D * const surface) {
    clean_up_surface(surface);
    free(surface->vertices);
    free(surface->normals);
    free(surface->faces);
}

/* Set up a Sphere struct */
void set_up_sphere(int make_sphere, double * const sphere_c, double sphere_r,
        Material M, int surf_index, AnalytSphere * const sph) {
    sph->surf_index = surf_index;
    sph->sphere_c = sphere_c;
    sph->sphere_r = sphere_r;
    sph->make_sphere = make_sphere;
    sph->material = M;
}

void set_up_circle(int make_circle, double * const circle_c, double circle_r,
		double * const circle_n, Material M, int surf_index, Circle * const circ) {
	circ->surf_index = surf_index;
	circ->centre = circle_c;
	circ->r = circle_r;
	circ->normal = circle_n;
	circ->make_circle = make_circle;
	circ->material = M;
}

void generate_empty_sphere(int surf_index, AnalytSphere * const sph) {
    double c[3] = {0, 0, 0};
    Material standard_mat;
    standard_mat.name = "diffuse";
    standard_mat.func_name = "cosine";
    standard_mat.params = 0;
    standard_mat.n_params = 0;
    standard_mat.func = distribution_by_name("cosine");

    sph->surf_index = surf_index;
    sph->sphere_c = c;
    sph->sphere_r = 0;
    sph->make_sphere = 0;
    sph->material = standard_mat;
}


/*
 * Initialise a Material with given props. The function name will also
 * be resolved, i.e. the scattering distributions will be searched by name
 * and the distribution function will be assigned to the func field.
 */
void set_up_material(char * const name, char * const function, double * const params, int n_params,
        Material * const mat) {
    mat->name = name;
    mat->func_name = function;
    mat->params = params;
    mat->n_params = n_params;

    mat->func = distribution_by_name(mat->func_name);
    /*if(mat->func == NULL) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:material",
                          "Distribution name %s could not be resolved.", mat->func_name);
    }*/
}


/*
 * Set up a struct of rays using the input vectors from MATLAB.
 *
 */
void compose_rays3D(double const ray_pos[], double const ray_dir[], int nrays, Rays3D * const all_rays) {
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
        rays[i].status = 0;
        rays[i].detector = 0;
    }

    /* Put the data into the struct */
    all_rays->rays = rays;
    all_rays->nrays = nrays;
}

/* Updates the ray postion */
void update_ray_position(Ray3D * const the_ray, const double new_pos[3]) {
    for (int j = 0; j < 3; j++)
        the_ray->position[j] = new_pos[j];
}

/* Updates the ray direction */
void update_ray_direction(Ray3D * const the_ray, const double new_dir[3]) {
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
void get_element3D(Surface3D const * const sample, int idx, Triangle * element) {
    int j;
    int vertices[3];

    /* Get the indices to the three vertices and get the surface normal */
    element->tri_index = idx;
    j = idx*3 + 0;
    vertices[0] = (sample->faces[j] - 1)*3;
    element->normal[0] = sample->normals[j];
    j += 1;
    vertices[1] = (sample->faces[j] - 1)*3;
    element->normal[1] = sample->normals[j];
    j += 1;
    vertices[2] = (sample->faces[j] - 1)*3;
    element->normal[2] = sample->normals[j];
    j = idx*6 + 0;
    for (int i = 0; i < 6; i++)
        element->lattice[i] = sample->lattice[j + i];
    
    /* Vertices of the triangle */
    element->v1[0] = sample->vertices[vertices[0]];
    element->v1[1] = sample->vertices[vertices[0] + 1];
    element->v1[2] = sample->vertices[vertices[0] + 2];
    element->v2[0] = sample->vertices[vertices[1]];
    element->v2[1] = sample->vertices[vertices[1] + 1];
    element->v2[2] = sample->vertices[vertices[1] + 2];
    element->v3[0] = sample->vertices[vertices[2]];
    element->v3[1] = sample->vertices[vertices[2] + 1];
    element->v3[2] = sample->vertices[vertices[2] + 2];
}

/* Cleanup a struct of rays */
void clean_up_rays(Rays3D all_rays) {
    free(all_rays.rays);
}

/* Get the ray positions and put them in an output array */
void get_positions(Rays3D const * const all_rays, double * const final_pos) {
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

void get_positions_indexed(Rays3D const * const all_rays, bool const * const index,
		double * const final_pos) {
    int i;
    int j = 0;

    /* Loop through all the rays */
    for (i = 0; i < all_rays->nrays; i++) {
        int k;
        Ray3D *current_ray;

        if (index[i]) {
        	current_ray = &all_rays->rays[i];
        	for (k = 0; k < 3; k++) {
        		int n;
        		n = k + 3*j;
        		final_pos[n] = current_ray->position[k];
        	}
        	j++;
        }
    }
}

/* Get the ray directions and put them in an output array */
void get_directions(Rays3D const * const all_rays, double * const final_dir) {
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

void get_directions_indexed(Rays3D const * const all_rays, bool * const index,
		double * const final_dir) {
    int i;
    int j = 0;

    /* Loop through all the rays */
    for (i = 0; i < all_rays->nrays; i++) {
        int k;
        Ray3D *current_ray;

        current_ray = &all_rays->rays[i];
        for (k = 0; k < 3; k++) {
            int n;
            n = k + 3*j;
            final_dir[n] = current_ray->direction[k];
        }
        j++;
    }
}

/* Get the number of scattering events per ray */
void get_scatters(Rays3D const * const all_rays, int * const nScatters) {
    int i;

    for (i = 0; i < all_rays->nrays; i++) {
        Ray3D *current_ray;

        current_ray = &all_rays->rays[i];
        nScatters[i] = (int)current_ray->nScatters;
    }
}

/* print details of Material struct */
void print_material(Material const * const mat) {
    printf("\tMAT %-10s func %-12s", mat->name, mat->func_name);
    for (int i = 0; i < mat->n_params; i++)
        printf(" %.2f ", mat->params[i]);
}

/* Print details of whole sample to console */
void print_surface(Surface3D const * const s) {
    printf("Surface index = %i\n", s->surf_index);
    for (int ivert = 0; ivert < s->n_vertices; ivert++) {
    	double v[3];
    	get_vertex(s, ivert, v);
        //int ind0, ind1, ind2;
        //lin(ivert, 0, &ind0);
        //lin(ivert, 1, &ind1);
        //lin(ivert, 2, &ind2);
        //printf("\n VERT % .2f % .2f % .2f", s->vertices[ind0],
        //    s->vertices[ind1], s->vertices[ind2]);
    	printf("\n VERT % .2f % .2f % .2f", v[0], v[1], v[2]);
    }
    for (int iface = 0; iface < s->n_faces; iface++) {
        int32_t f[3];
        double n[3];
        double b[6];
    	printf("\n FACE %2d", iface);
        get_face(s, iface, f);
        get_normal(s, iface, n);
        get_lattice(s, iface, b);
        //int ind0, ind1, ind2;
        //lin(iface, 0, &ind0);
        //lin(iface, 1, &ind1);
        //lin(iface, 2, &ind2);
        printf("\tV %3d %3d %3d\n", f[0], f[1], f[2]);
        //s->faces[ind0], s->faces[ind1], s->faces[ind2]);
        printf("\tN % .2f % .2f % .2f\n", n[0], n[1], n[2]);
        printf("\tB % .2f % .2f % .2f % .2f % .2f % .2f\n", b[0], b[1], b[2],
               b[3], b[4], b[5]);
        //s->normals[ind0],s->normals[ind1], s->normals[ind2]);
        print_material(s->compositions[iface]);
    }
    printf("\n");
}

/* Prints all the information about the ray to the terminal */
void print_ray(Ray3D const * const the_ray) {
    printf("\non_element = %i\n", the_ray->on_element);
    printf("on_surface = %i\n", the_ray->on_surface);
    printf("nScatters = %i\n", the_ray->nScatters);
    printf("Position: ");
    print1D_double(the_ray->position, 3);
    printf("Direction: ");
    print1D_double(the_ray->direction, 3);
}

/* Prints all the information about a BackWall struct */
void print_BackWall(BackWall const * const wall) {
    printf("\nSurface index = %i\n", wall->surf_index);
    printf("Plate represent = %i\n", wall->plate_represent);
    printf("Aperture centre: ");
    print1D_double(wall->aperture_c, 2);
    printf("Aperture axes: ");
    print1D_double(wall->aperture_axes, 2);
    printf("Radius of the plate = %f\n", wall->circle_plate_r);
    print_material(&(wall->material));
    printf("\n");
}

/* Prints all the information on all the apertues in the NBackWall struct */
void print_nBackWall(NBackWall const * const all_apertures) {
    printf("\nNumber of apertures = %i\n", all_apertures->n_detect);
    for (int i = 0; i < all_apertures->n_detect; i++) {
        printf("Aperture %i:\n", i);
        printf("Centre: ");
        print1D_double(&all_apertures->aperture_c[2*i], 2);
        printf("Axes: ");
        print1D_double(&all_apertures->aperture_axes[2*i], 2);
    }
    // print_material(&(all_apertures->material));
    printf("\n");
}


/* Print the position, radius, material etc of a sphere */
void print_sphere(AnalytSphere const * const sphere){
    printf("\n\t Sphere make %d \n R %3.3f\n C %3.3f %3.3f %3.3f\n", sphere->make_sphere,
              sphere->sphere_r, sphere->sphere_c[0],
              sphere->sphere_c[1], sphere->sphere_c[2]);
    print_material(&(sphere->material));
}

void print_circle(Circle const * const circle) {
	printf("\n\t Circle make %d \n R %3.3f\n C %3.3f %3.3f %3.3f\n N %3.3f %3.3f %3.3f\n",
			circle->make_circle, circle->r, circle->centre[0], circle->centre[1],
			circle->centre[2], circle->normal[0], circle->normal[1], circle->normal[2]);
}

void print_triangle(Triangle const * const tri) {
	printf("Triangle element:");
	printf("\n VERT1 % .2f % .2f % .2f", tri->v1[0], tri->v1[1], tri->v1[2]);
	printf("\n VERT2 % .2f % .2f % .2f", tri->v2[0], tri->v2[1], tri->v2[2]);
	printf("\n VERT3 % .2f % .2f % .2f", tri->v3[0], tri->v3[1], tri->v3[2]);
	printf("\n NORM  % .2f % .2f % .2f\n", tri->normal[0], tri->normal[1], tri->normal[2]);
}

/*
 * Gets the centre and the axes of the nth detector in the series and returns
 * a struct to that aperture.
 */
void get_nth_aperture(int n, NBackWall const * const allApertures, BackWall * const this_wall) {
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
void create_ray(Ray3D * const gen_ray, SourceParam const * const source, MTRand * const myrng) {
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
            cosine_scatter(normal, NULL, NULL, gen_ray->direction, NULL, myrng);
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

    // TODO: use the new_Ray function

    /* Initialise other elements of the ray struct */
    gen_ray->on_surface = -1;
    gen_ray->on_element = -1;
    gen_ray->nScatters = 0;
    gen_ray->status = 0;
    gen_ray->detector = 0;
}

void new_Ray(Ray3D * const gen_Ray, double const pos[3], double const dir[3]) {
	int i;

	for (i = 0; i < 3; i++) {
		gen_Ray->position[i] = pos[i];
		gen_Ray->direction[i] = dir[i];
	}
	gen_Ray->on_element = -1;
	gen_Ray->on_surface = -1;
	gen_Ray->nScatters = 0;
	gen_Ray->status = 0;
	gen_Ray->detector = 0;
}

/*
 * Creates a basic flat sample with 3 triangles, useful for debugging. If this
 * function is used rather than the conversion from mxArrays then
 * clean_up_surface_all_arrays() must be used to free up memory.
 *
 * INPUTS:
 *  sample_index - the surface index of the sample in use
 *  size         - length of the sides of the surface
 *  sample       - pointer to Surface3D struct to be initialised
 */
void make_basic_sample(int sample_index, double size, Surface3D * const sample) {

    int num_materials = 1;
    int i;
    Material standard_mat;
    standard_mat.name = "shiny";
    standard_mat.func_name = "pure_specular";
    standard_mat.params = 0;
    standard_mat.n_params = 0;
    standard_mat.func = distribution_by_name("pure_specular");

    int ntriag_sample = 3;
    int nvert = 5;
    double * V;
    double * N;
    double * B;
    int32_t * F;
    V = (double *)malloc(sizeof(double)*nvert*3);
    double VV[15] = {-size/2, -1, -size/2,
                     -size/2, -1,  size/2,
                      0,      -1,  size/2,
                      size/2, -1, -size/2,
                      size/2, -1,  size/2};
    for (i = 0; i < nvert*3; i++) {
        V[i] = VV[i];
    }
    double NN[9] = {0,1,0, 0,1,0, 0,1,0};
    double BB[18] = {1,0,0, 0,0,1,
                     1,0,0, 0,0,1,
                     1,0,0, 0,0,1};
    int32_t FF[9] = {1,2,3, 1,3,4, 3,5,4};
    N = (double *)malloc(sizeof(double)*ntriag_sample*3);
    B = (double *)malloc(sizeof(double)*ntriag_sample*6);
    F = (int32_t *)malloc(sizeof(int32_t)*ntriag_sample*3);
    for (i = 0; i < ntriag_sample*3; i++) {
        N[i] = NN[i];
        F[i] = FF[i];
        B[i] = BB[i];
    }
    char *C[3] = {"shiny", "shiny", "shiny"};
    Material * M;
    M = (Material *)malloc(sizeof(Material)*ntriag_sample);
    Material MM[3] = {standard_mat, standard_mat, standard_mat};
    for (i = 0; i < ntriag_sample; i++) {
        M[i] = MM[i];
    }
    set_up_surface(V, N, B, F, C, M, num_materials, ntriag_sample,
            nvert, sample_index, sample);
}


/*
 * Adds a 3 element array to another 3 element array multiplied by a scalar
 * (propagates an array).
 *
 * EQUATION:
 *  result = init + direc*a
 *
 * INPUTS:
 *  init   - double array, vector 1 (initial position)
 *  direc  - double array, vector 2 (direction)
 *  a      - double, amount to multiply vector 2 by (move in the direction of
 *           vector direc by this amount.)
 *  result - double array, array to store the result (final position) in
 */
void propagate(const double init[3], const double direc[3], double a, double result[3]) {
    int i;
    for (i = 0; i < 3; i++) {
        result[i] = init[i] + direc[i]*a;
    }
}

/*
 * Calculates the dot product of two three element vectors passed as arrays.
 *
 * INPUTS:
 *  a - double array, first vector
 *  b - double array, second vector
 *  result - double pointer, the value of the dot product
 */
void dot(const double a[3], const double b[3], double * const result) {
    int i;
    *result = 0;
    for (i = 0; i < 3; i++) {
        *result += a[i]*b[i];
    }
}

/* Calculates the cross product of two 3-vectors and writes it to c */
void cross(const double a[3], const double b[3], double c[3]) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

/*
 * Returns the square of the norm of a three vector.
 *
 * INPUTS:
 *  vect - double array, three element array to calculate the norm^2 of.
 *
 * OUTPUTS:
 *  result - the value of the norm^2
 */
void norm2(const double vect[3], double * const result) {
    *result = 0;
    int i;
    for (i = 0; i < 3; i++) {
        *result += vect[i]*vect[i];
    }
}

/*
 * Normalizes a 3 element vector. Stores the normalized vector in the same
 * array as the original.
 *
 * INPUTS:
 *  vect - double array, three element array to normalize
 */
void normalise(double vect[3]) {
    double magnitude;
    int i;

    norm2(vect, &magnitude);
    for (i = 0; i < 3; i++) {
        vect[i] = vect[i]/sqrt(magnitude);
    }
}


/*
 * Reflects a vector through the normal provided storing the result in the
 * array provided.
 *
 * INPUTS:
 *  normal   - double array, the normal to the surface at the point of
 *             reflection
 *  init_dir - double array, the initial direction of the ray
 *  new_dir  - double array, an array to store the final direction of the ray
 */
void reflect3D(const double normal[3], const double init_dir[3], double new_dir[3]) {
    double tmp;
    dot(normal, init_dir, &tmp);
    propagate(init_dir, normal, -2*tmp, new_dir);
    normalise(new_dir);
}


/*
 * Find two (unit) directions perpendicular to a given unit vector n,
 * and write them to v1 and v2.
 *
 * The idea is that (x, y, z) . (-y, x, 0) = 0
 * but have to check that x and y are not both 0 to avoid the zero vector
 * credit https://math.stackexchange.com/q/137362/558299
 */
void perpendicular_plane(const double n[3], double v1[3], double v2[3]) {
    double epsilon = 1e-3;

    // if both nx and ny are zero, pick nx and nz to switch
    if(fabs(n[0]) < epsilon && fabs(n[1]) < epsilon) {
        v1[0] = n[2];
        v1[1] = 0;
        v1[2] = - n[0];
    }
    else { // if at least one of nx and ny is nonzero, switch them around
        v1[0] = n[1];
        v1[1] = - n[0];
        v1[2] = 0;
    }

    normalise(v1);

    // cross product for the second vector
    // no need to normalise this as the lengths of n and v1 are 1
    v2[0] = n[1]*v1[2] - n[2]*v1[1];
    v2[1] = n[2]*v1[0] - n[0]*v1[2];
    v2[2] = n[0]*v1[1] - n[1]*v1[0];
}

/*
 * Performs a matrix mutiplication:
 *  result = M*vec
 */
void matrix_mult(const double vec[3], double M[3][3], double result[3]) {
	result[0] = M[0][0]*vec[0] + M[0][1]*vec[1] + M[0][2]*vec[2];
	result[1] = M[1][0]*vec[0] + M[1][1]*vec[1] + M[1][2]*vec[2];
	result[2] = M[2][0]*vec[0] + M[2][1]*vec[1] + M[2][2]*vec[2];
}

/*
 * Performs a general roation of the provided vector about a specified axis
 * through the origin. Note that it takes the cosine of the angle as an argument
 */
void general_rotation(const double vec[3], const double axis[3],
		double result[3], double c_theta) {
	double R[3][3];
	double s_theta, c_th_1;
	double xy, xz, yz;

	s_theta = sqrt(1 - c_theta*c_theta);
	c_th_1 = 1 - c_theta;
	xy = axis[0]*axis[1];
	xz = axis[0]*axis[2];
	yz = axis[1]*axis[2];

	// First row
	R[0][0] = c_theta + axis[0]*axis[0]*c_th_1;
	R[0][1] = xy*c_th_1 - axis[3]*s_theta;
	R[0][2] = xz*c_th_1 + axis[1]*s_theta;

	// Second row
	R[1][0] = xy*c_th_1 + axis[2]*s_theta;
	R[1][1] = c_theta + axis[1]*axis[1]*c_th_1;
	R[1][2] = yz*c_th_1 - axis[0]*s_theta;

	// Third row
	R[2][0] = xz*c_th_1 - axis[1]*s_theta;
	R[2][1] = yz*c_th_1 + axis[0]*s_theta;
	R[2][2] = c_theta + axis[2]*axis[2]*c_th_1;

	matrix_mult(vec, R, result);
}


/*
 * Solves a 3D matrix equation Au=v for u, we write
 *      (a b c)   (u[0])   (j)
 *      (d e f) . (u[1]) = (k)
 *      (g h i)   (u[2])   (l)
 * Uses Cramer's rule and explicit formula for the determinants to solve the
 * equation. Checks to see if the determinant of the matrix is less than the
 * supplied tolerance epsilon, returns 0 if it is less and 1 if it is more than
 * the tolerance.
 *
 * INPUTS:
 *  A       - double array, the 3x3 matrix, as an array, for the equation
 *  u       - double array, a 3 element array for the answer to be put in
 *  v       - double array, the 3-vector on the other side of the equation
 *  epsilon - double, the tolerance for the size of the determinant of the
 *            matrix A
 *  success - int, sets to 0 or 1 depending on the size of the determinant of A
 */
void solve3x3(double A[3][3], double u[3], double v[3], double epsilon, int * const success) {
    double M, Dx, Dy, Dz, X1, X2, X3;
    /*double a,b,c,d,e,f,g,h,i,j,k,l;

    a = A[0][0];
    b = A[0][1];
    c = A[0][2];

    d = A[1][0];
    e = A[1][1];
    f = A[1][2];

    g = A[2][0];
    h = A[2][1];
    i = A[2][2];

    j = v[0];
    k = v[1];
    l = v[2];*/

    /* These are used twice - only compute them once */
    /*X1 = e*i - h*f;
    X2 = c*h - b*i;
    X3 = b*f - c*e;
    M = a*X1 + d*X2 + g*X3;*/

    X1 = A[1][1]*A[2][2] - A[2][1]*A[1][2];
    X2 = A[0][2]*A[2][1] - A[0][1]*A[2][2];
    X3 = A[0][1]*A[1][2] - A[0][2]*A[1][1];
    M = A[0][0]*X1 + A[1][0]*X2 + A[2][0]*X3;

    /* fabs() is the math.h abs function for floats */
    if (fabs(M) < epsilon) {
        *success = 0;
        return;
    }
    /*Dx = j*X1 + k*X2 + l*X3;

    X1 = a*k - j*d;
    X2 = j*g - a*l;
    X3 = d*l - k*g;

    Dy = i*X1 + f*X2 + c*X3;
    Dz = - h*X1 - e*X2 - b*X3;*/

    Dx = v[0]*X1 + v[1]*X2 + v[2]*X3;

    X1 = A[0][0]*v[1] - v[0]*A[1][0];
    X2 = v[0]*A[2][0] - A[0][0]*v[2];
    X3 = A[1][0]*v[2] - v[1]*A[2][0];

    Dy = A[2][2]*X1 + A[1][2]*X2 + A[0][2]*X3;
    Dz = - A[2][1]*X1 - A[1][1]*X2 - A[0][1]*X3;

    u[0] = Dx/M;
    u[1] = Dy/M;
    u[2] = Dz/M;

    *success = 1;
}

void get_elementPlane(Plane const * const plane, Triangle * const element) {
	int i;

	element->tri_index = -1;
	for (i = 0; i < 3; i++) {
		element->v1[i] = plane->example_point1[i];
		element->v2[i] = plane->example_point2[i];
		element->v3[i] = plane->example_point3[i];
		element->normal[i] = plane->normal[i];
	}
}

void get_elementCircle(Circle the_circle, Triangle * const element) {
    int i;
    double v1[3], v2[3];

    element->tri_index = -1;
	perpendicular_plane(the_circle.normal, v1, v2);
	for (i = 0; i < 3; i++) {
    	element->normal[i] = the_circle.normal[i];
    	element->v1[i] = the_circle.centre[i];
    	element->v2[i] = the_circle.centre[i] + v1[i];
    	element->v3[i] = the_circle.centre[i] + v2[i];
    }
}

void constructPlate(double * point, double * normal, Plane * plane, int index, Material material) {
	double point2[3], point3[3];


	plane->surf_index = index;
	plane->example_point1 = point;
	plane->example_point2 = point2;
	plane->example_point3 = point3;
	plane->normal = normal;
	plane->material = material;
}

static void get_normal_ptr(Surface3D const * const s, int ind, double ** n){
	int ind0;
    lin(ind, 0, &ind0);

    *n = &s->normals[ind0];
}

static void get_lattice_ptr(Surface3D const * const s, int ind, double ** b) {
    int ind0;
    ind0 = 6*ind;
    
    *b = &s->lattice[ind0];
}

void get_lattice(Surface3D const * const s, int ind, double b[6]) {
    double * b_ptr;
    int j;
    
    get_lattice_ptr(s, ind, &b_ptr);
    
    for (j = 0; j < 6; j++)
        b[j] = b_ptr[j];
}

void get_normal(Surface3D const * const s, int ind, double n[3]) {
	double * n_ptr;
	int j;

	get_normal_ptr(s, ind, &n_ptr);

	for (j = 0; j < 3; j++)
		n[j] = n_ptr[j];
}

void set_normal(Surface3D const * const s, int ind, double new_n[3]) {
	double * n_ptr;
	int j;

	get_normal_ptr(s, ind, &n_ptr);

	for (j = 0; j < 3; j++)
		n_ptr[j] = new_n[j];
}

static void get_face_ptr(Surface3D const * const s, int ind, int32_t ** f){
	int ind0;
    lin(ind, 0, &ind0);

    *f = &s->faces[ind0];
}

void get_face(Surface3D const * const s, int ind, int32_t f[3]) {
	int32_t * f_ptr;
	int j;

	get_face_ptr(s, ind, &f_ptr);

	for (j = 0; j < 3; j++)
		f[j] = f_ptr[j];
}

void set_face(Surface3D const * const s, int ind, int32_t new_f[3]) {
	int32_t * f_ptr;
	int j;

	get_face_ptr(s, ind, &f_ptr);

	for (j = 0; j < 3; j++)
		f_ptr[j] = new_f[j];
}

/*
 * Gets a pointer to the desired vertex of a triangulated surface.This will allow editing of
 * the vertex.
 */
static void get_vertex_ptr(Surface3D const * const s, int ind, double ** v){
	int ind0;
    lin(ind, 0, &ind0);

    *v = &s->vertices[ind0];
}

/*
 * Gets the values of a vertex in a new variable.
 */
void get_vertex(Surface3D const * const s, int ind, double v[3]) {
	double * v_ptr;
	int j;

	get_vertex_ptr(s, ind, &v_ptr);

	for (j = 0; j < 3; j++)
		v[j] = v_ptr[j];
}

void set_vertex(Surface3D const * const s, int ind, double new_v[3]) {
	double * v_ptr;
	int j;

	get_vertex_ptr(s, ind, &v_ptr);

	for (j = 0; j < 3; j++)
		v_ptr[j] = new_v[j];
}

void moveSurface(Surface3D * const s, double displace[3]) {
	int ivert;

	for (ivert = 0; ivert < s->n_vertices; ivert++) {
		double * v;
		int j;
		get_vertex_ptr(s, ivert, &v);

		for (j = 0; j < 3; j++) {
			v[j] += displace[j];
		}
	}
}

