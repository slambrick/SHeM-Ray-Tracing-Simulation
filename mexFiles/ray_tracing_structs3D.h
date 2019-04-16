/* 
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 *
 * 
 * Structues for use with the SHeM Ray Tracing Simulation.
 */

#ifndef _ray_tracing_structures_h
#define _ray_tracing_structures_h

#include "mex.h"
#include <gsl/gsl_rng.h>

/******************************************************************************/
/*                          Structure declarations                            */
/******************************************************************************/

/* 
 * A structure for holding information on a 3D sample surface constructed of
 * planar triangles.
 */
typedef struct _surface3d {
    int surf_index;      /* Index of the surface */
    int n_elements;      /* Number of elements in the surface */
    double *vertices;    /* Vertices of the surface */
    double *normals;     /* Normals to the elements of the surface */
    double *faces;       /* Faces of the surface. */
    double *composition; /* The type of scattering off the elements of this surface */
    double *scattering_parameters;
        /* An array of parameters for scattering off the surface  */
} Surface3D;

/* Information on the flat plate model of detection */
typedef struct _backWall {
    int surf_index;         /* Index of this surface */
    double *aperture_c;   /* The centre of the detector aperture */
    double *aperture_axes;/* The two axes of the elliptical aperture */
    double circle_plate_r;  /* The radius of the plate that the aperture sits on */
    double composition;
    double scattering_parameters;
    int plate_represent;    /* Should the plate be scattered off, 0 or 1 */
} BackWall;

/* Contains information on a whole series of back wall apertures */
typedef struct _nBackWall{
    int surf_index;
    int n_detect;
    double *aperture_c;
    double *aperture_axes;
    double circle_plate_r;
    double composition;
    double scattering_parameters;
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
    int surf_index;
    double scan_pos_x;             /* The x position of the scan */
    double scan_pos_z;             /* The z position of the scan */
    double dist_to_sphere;         /* The working distance */
    double sphere_r;               /* The radius of the sphere */
    double composition;     /* The type of scattering off the sphere */
    double scattering_parameters;  /* Any scattering parameters for the sphere */
    int make_sphere;               /* Does the sphere actually exsist? */
} AnalytSphere;


/* A structure for holding a single ray */
typedef struct _ray3d {
    double position[3];   /* First vertex of the surface element */
    double direction[3];  /* Second vertex of the surface element */
    int32_t nScatters;    /* The number of scattering events the ray has undergone */
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
Surface3D set_up_surface(double V[], double N[], double F[], double C[], 
        double P[], int ntraig, int surf_index);

/* Set up a Sphere struct */
AnalytSphere set_up_sphere(int make_sphere, double scan_pos_x, double scan_pos_z, 
        double dist_to_sample, double sphere_r, double sphere_scattering, 
        double sphere_parameters, int surf_index);

/* Sets up a struct containing information on all the rays to be simulated. */
Rays3D compose_rays3D(double ray_pos[], double ray_dir[], int nrays);

/* Updates the ray position */
void update_ray_position(Ray3D *the_ray, double new_pos[3]);

/* Gets an element of a triangulated surface */
void get_element3D(Surface3D *Sample, int n, double v1[3], double v2[3], 
        double v3[3], double nn[3]);

/* Cleans up the allocated memory inside the ray struct */
void clean_up_rays(Rays3D all_rays);

/* Get the final possitions */
void get_positions(Rays3D *all_rays, double *final_pos);

/* Get the final directions */
void get_directions(Rays3D *all_rays, double *final_dir);

/* Get the number of scattering events for the rays */
void get_scatters(Rays3D *all_rays, int32_t *nScatters);

/* Prints information on a ray to the terminal */
void print_ray(Ray3D *the_ray);

/* Prints information on a BackWall to the terminal */
void print_BackWall(BackWall *wall);

/* Prints the centres and axes of all the apertures in the NBackWall struct */
void print_nBackWall(NBackWall *all_apertures);

/* Gets information on one aperure out of a series of apertures */
void get_nth_aperture(int n, NBackWall *allApertures, BackWall *this_wall);

/* Creates a ray in the pinhole */
Ray3D create_ray_source(double pinhole_r, double *pinhole_c, double theta_max, 
        double init_angle, int source_model, gsl_rng *my_rng, double sigma);

#endif
