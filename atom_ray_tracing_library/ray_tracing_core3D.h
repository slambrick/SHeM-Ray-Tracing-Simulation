/*
 * Copyright (c) 2018-20, Sam Lambrick. 2020 Dan Serment.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 *
 * Structues for use with the SHeM Ray Tracing Simulation.
 */

#ifndef _ray_tracing_structures3D_h
#define _ray_tracing_structures3D_h

#include "mtwister.h"
#include <stdint.h>
#include <stdbool.h>
#include "distributions3D.h"

/******************************************************************************/
/*                          Structure declarations                            */
/******************************************************************************/


/* Structure holding information of a material's scattering properties */
typedef struct _material {
    char * name;        // the name (key) of the material
    char * func_name;   // the name of scattering distribution function
    double * params;    // extra parameters to pass to the function
    int n_params;       // number of parameters -- must be length of params[]
    distribution_func func; // the actual scattering probability distribution
} Material;

/* General surface type for the others to inherit from? */
typedef struct _generalSurface {
	int surf_index;
	Material material;
} GeneralSurf;

/*
 * A structure for holding information on a 3D sample surface constructed of
 * planar triangles.
 */
typedef struct _surface3d {
    int surf_index;       /* Index of the surface */
    int n_faces;          /* Number of faces in the surface */
    int n_vertices;       /* Number of vertices defined */
    double * vertices;    /* Vertices of the surface */
    int * faces;          /* Faces of the surface. */
    double * normals;     /* Normals to the elements of the surface */
    double * lattice;     /* Reciprocal lattice parameters of the surface elements */
    Material ** compositions; /* The type of scattering off the elements of this surface */
} Surface3D;

/* Contains information on a single triangle */
typedef struct _triangle {
	double v1[3];
	double v2[3];
	double v3[3];
	double normal[3];
    double lattice[6];
	int tri_index;
} Triangle;

/* Information on the flat plate model of detection */
typedef struct _backWall {
    int surf_index;         /* Index of this surface */
    double *aperture_c;     /* The centre of the detector aperture */
    double *aperture_axes;  /* The two axes of the elliptical aperture */
    double circle_plate_r;  /* The radius of the plate that the aperture sits on */
    Material material;
    int plate_represent;    /* Should the plate be scattered off, 0 or 1 */
} BackWall;

/* Contains information on a whole series of back wall apertures */
typedef struct _nBackWall{
    int surf_index;
    int n_detect;
    double *aperture_c;
    double *aperture_axes;
    double circle_plate_r;
    Material material;
    int plate_represent;
} NBackWall;

/* Information on the abstract hemisphere model of detection */
typedef struct _abstractHemi {
    int surf_index;
    double aperture_theta;
    double aperture_phi;
    double half_cone_angle;
    double dist_to_sample;
} AbstractHemi;

/* A structure to hold information on the analytic sphere */
typedef struct _sphere {
    int make_sphere;               /* Does the sphere actually exist? */
    int surf_index;
    double sphere_r;               /* The radius of the sphere */
    double  *sphere_c;             /* x y z of the centre of the sphere */
    Material material;           /* The material composition of the sphere */
} AnalytSphere;

/* Information on a circular surface */
typedef struct _circle {
	int make_circle;
	int surf_index;
	double r;
	double *centre;
	double *normal;
	Material material;
} Circle;

/* The total information on a potential sample */
typedef struct _sample {
    Circle *the_circle;
    AnalytSphere *the_sphere;
    Surface3D *triag_sample;
} Sample;

/* Information on a plane */
typedef struct _plane {
	int surf_index;
	double *normal;
	double *example_point1;
	double *example_point2;
	double *example_point3;
	Material material;
} Plane;

/* A structure for holding a single ray */
typedef struct _ray3d {
    double position[3];   /* First vertex of the surface element */
    double direction[3];  /* Second vertex of the surface element */
    int nScatters;        /* The number of scattering events the ray has undergone */
    int on_element;       /* The index of the surface element that the ray is on */
    int on_surface;       /* The index of the surface that the ray is on */
    int status;           /* Is the ray alive (0), dead (1), or detected (2) */
    int detector;         /* If the ray is detected, which one? none (0) */
} Ray3D;

/* A structure to hold an array of Ray3D structs */
typedef struct _rays3d {
    Ray3D *rays;  /* Pointer to an array of Ray structs */
    int nrays;  /* The number of rays that are stored in the struct */
} Rays3D;

// Contains the source information
typedef struct _sourceParam {
	double pinhole_r;
	double pinhole_c[3];
	double theta_max;
	double init_angle;
	int source_model;
	double sigma;
} SourceParam;

// Basic information on a scattering position
typedef struct _scatterInfo {
    double nearest_inter[3]; // Position of the nearest intersection
    double nearest_n[3];     // Normal to the surface at the intersection
    double nearest_b[6];     // Lattice parameters at the intersection
    int meets;               // Has a surface been met?
    int tri_hit;             // Index of the triangle that has been hit
    int which_surface;       // Index of the surface that has been hit
} ScatterInfo;

/******************************************************************************/
/*                           Function declarations                            */
/******************************************************************************/

/*  Set up a surface containing the information on a triangulated surface. */
void set_up_surface(double V[], double N[], double B[], int32_t F[], char * C[], Material M[],
		int nmaterials, int ntriag, int nvert, int surf_index, Surface3D * const surf);

void clean_up_surface(Surface3D * const surface);

void clean_up_surface_all_arrays(Surface3D * const surface);

/* Set up a Sphere struct */
void set_up_sphere(int make_sphere, double * const sphere_c, double sphere_r,
        Material M, int surf_index, AnalytSphere * const sph);

void set_up_circle(int make_circle, double * const circle_c, double circle_r,
		double * const circle_n, Material M, int surf_index, Circle * const circ);

void generate_empty_sphere(int surf_index, AnalytSphere * const sph);

/* Initialise a Material with given properties */
void set_up_material(char * const name, char * const function, double * const params, int n_params,
		Material * const mat);

/* Sets up a struct containing information on all the rays to be simulated. */
void compose_rays3D(double const ray_pos[], double const ray_dir[], int nrays, Rays3D * const all_rays);

/* Updates the ray position */
void update_ray_position(Ray3D * const the_ray, double const new_pos[3]);

/* Updates the ray direction */
void update_ray_direction(Ray3D * const the_ray, double const new_dir[3]);

/* Gets an element of a triangulated surface. Vertices are written to
 * v1, v2, v3, and the normal to `normal`.
 */
void get_element3D(Surface3D const * const sample, int index, Triangle * const element);

/* Cleans up the allocated memory inside the ray struct */
void clean_up_rays(Rays3D all_rays);

/* Get the final possitions */
void get_positions(Rays3D const * const all_rays, double * const final_pos);

void get_positions_indexed(Rays3D const * const all_rays, bool const * const index,
		double * const final_pos);

/* Get the final directions */
void get_directions(Rays3D const * const all_rays, double * const final_dir);

void get_directions_indexed(Rays3D const * const all_rays, bool * const index,
		double * const final_dir);

/* Get the number of scattering events for the rays */
void get_scatters(Rays3D const * const all_rays, int * const nScatters);

/* Print details of material */
void print_material(Material const * const mat);

/* Print details of whole sample to console */
void print_surface(Surface3D const * const s);

/* Prints information on a ray to the terminal */
void print_ray(Ray3D const * const the_ray);

void print_BackWall(BackWall const * const wall);

/* Prints the centres and axes of all the apertures in the NBackWall struct */
void print_nBackWall(NBackWall const * const all_apertures);

/* Print the position, radius, material etc of a sphere */
void print_sphere(AnalytSphere const * const sphere);

void print_circle(Circle const * const circle);

void print_triangle(Triangle const * const tri);

/* Gets information on one aperure out of a series of apertures */
void get_nth_aperture(int n, NBackWall const * const allApertures, BackWall * const this_wall);

/* Creates a ray in the pinhole */
void create_ray(Ray3D * const gen_ray, SourceParam const * const source, MTRand * const myrng);

void new_Ray(Ray3D * const gen_Ray, double const pos[3], double const dir[3]);

// Creates a flat sample with 3 triangles
void make_basic_sample(int sample_index, double size, Surface3D * const sample);

/* Adds a 3 element array to another 3 element array multiplied by a scalar
 * (propagates an array). */
void propagate(const double init[3], const double direc[3], double a, double result[3]);

/* Calculates the dot product of two 3 element double vectors. */
void dot(const double a[3], const double b[3], double * const result);

/* Calculates the cross product of two 3-vectors and writes it to c */
void cross(const double a[3], const double b[3], double c[3]);

/* Returns the square of the norm of a three vector. */
void norm2(const double vect[3], double * const result);

/* Normalises a three vector */
void normalise(double vect[3]);

/* reflect direction through normal */
void reflect3D(const double normal[3], const double init_dir[3], double new_dir[3]);

/* find two (unit) directions perpendicular to a given unit vector,
 * and write them to v1 and v2. */
void perpendicular_plane(const double n[3], double v1[3], double v2[3]);

void matrix_mult(const double vec[3], double M[3][3], double result[3]);

void general_rotation(const double vec[3], const double axis[3],
		double result[3], double c_theta);

/* Solves a 3D matrix equation Au=v. */
void solve3x3(double A[3][3], double u[], double v[], double epsilon, int * const success);

void get_elementPlane(Plane const * const plane, Triangle * const element);

void get_elementCircle(Circle * the_circle, Triangle * const element);

void constructPlate(double * point, double * normal, Plane * plane, int index, Material material);

void get_normal(Surface3D const * const s, int ind, double n[3]);

void get_lattice(Surface3D const * const s, int ind, double b[6]);

void set_normal(Surface3D const * const s, int ind, double new_n[3]);

void get_face(Surface3D const * const s, int ind, int32_t f[3]);

void set_face(Surface3D const * const s, int ind, int32_t new_f[3]);

void get_vertex(Surface3D const * const s, int ind, double v[3]);

void set_vertex(Surface3D const * const s, int ind, double new_v[3]);

void moveSurface(Surface3D * const s, double displace[3]);

void set_up_plate(int plate_represent, int n_detect, double circle_plate_r, 
        double * aperture_axes, double * aperture_c, Material M, int surf_index, 
        NBackWall * const plate);

#endif
