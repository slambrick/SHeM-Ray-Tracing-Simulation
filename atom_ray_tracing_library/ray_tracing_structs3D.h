/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 *
 * Structues for use with the SHeM Ray Tracing Simulation.
 */

#ifndef _ray_tracing_structures3D_h
#define _ray_tracing_structures3D_h

#include <mex.h>
#include "mtwister.h"
#include <stdint.h>
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


/*
 * A structure for holding information on a 3D sample surface constructed of
 * planar triangles.
 */
typedef struct _surface3d {
    int surf_index;       /* Index of the surface */
    int n_faces;          /* Number of faces in the surface */
    int n_vertices;       /* Number of vertices defined */
    double *vertices;     /* Vertices of the surface */
    int *faces;           /* Faces of the surface. */
    double *normals;      /* Normals to the elements of the surface */
    Material ** compositions; /* The type of scattering off the elements of this surface */
} Surface3D;

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


/* A structure for holding a single ray */
typedef struct _ray3d {
    double position[3];   /* First vertex of the surface element */
    double direction[3];  /* Second vertex of the surface element */
    int nScatters;        /* The number of scattering events the ray has undergone */
    int on_element;       /* The index of the surface element that the ray is on */
    int on_surface;       /* The index of the surface that the ray is on */
} Ray3D;

/* A structure to hold an array of Ray3D structs */
typedef struct _rays3d {
    Ray3D *rays;  /* Pointer to an array of Ray structs */
    int nrays;  /* The number of rays that are stored in the struct */
} Rays3D;

/******************************************************************************/
/*                           Function declarations                            */
/******************************************************************************/

/*  Set up a surface containing the information on a triangulated surface. */
Surface3D set_up_surface(double V[], double N[], int F[], char * C[],
        Material * materials, int nmaterials, int ntriag, int nvert, int surf_index);

void clean_up_surface(Surface3D * surface);

/* Set up a Sphere struct */
AnalytSphere set_up_sphere(int make_sphere, double *sphere_c, double sphere_r,
        Material M, int surf_index);

/* Initialise a Material with given properties */
Material set_up_material(char * name, char * function, double * params, int n_params);

/* Sets up a struct containing information on all the rays to be simulated. */
Rays3D compose_rays3D(double ray_pos[], double ray_dir[], int nrays);

/* Updates the ray position */
void update_ray_position(Ray3D *the_ray, const double new_pos[3]);

/* Updates the ray direction */
void update_ray_direction(Ray3D *the_ray, const double new_dir[3]);

/* Gets an element of a triangulated surface. Vertices are written to
 * v1, v2, v3, and the normal to `normal`.
 */
void get_element3D(Surface3D *sample, int index, double v1[3], double v2[3],
        double v3[3], double normal[3]);

/* Cleans up the allocated memory inside the ray struct */
void clean_up_rays(Rays3D all_rays);

/* Get the final possitions */
void get_positions(Rays3D *all_rays, double *final_pos);

/* Get the final directions */
void get_directions(Rays3D *all_rays, double *final_dir);

/* Get the number of scattering events for the rays */
void get_scatters(Rays3D *all_rays, int *nScatters);

/* Print details of material */
void print_material(const Material * mat);

/* Print details of whole sample to console */
void print_surface(const Surface3D * sample);

/* Prints information on a ray to the terminal */
void print_ray(Ray3D *the_ray);

/* Prints information on a BackWall to the terminal */
void print_BackWall(BackWall *wall);

/* Prints the centres and axes of all the apertures in the NBackWall struct */
void print_nBackWall(NBackWall *all_apertures);

/* Print the position, radius, material etc of a sphere */
void print_sphere(const AnalytSphere * sphere);

/* Gets information on one aperure out of a series of apertures */
void get_nth_aperture(int n, NBackWall *allApertures, BackWall *this_wall);

/* Creates a ray in the pinhole */
void create_ray(Ray3D * gen_ray, double pinhole_r, const double *pinhole_c, double theta_max,
        double init_angle, int source_model, double sigma, MTRand *myrng);

#endif
