/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the  
 * GNU/GPL-3.0-or-later.
 * 
 * Contains structs used in the 2D ray tracing off surfaces.
 */
#include "mtwister.h"

#ifndef _ray_tracing_structs2D_h
#define _ray_tracing_structs2D_h

/* A structure for holding information on a 2D sample surface */
typedef struct _surface2d {
    int n_elements; /* Number of elements in the surface */
    double *V;      /* Vertices of the surface */
    double *N;      /* Normals to the elements of the surface */
    int scattering; /* The type of scattering off the surface */
    double *scattering_parameters;
        /* An array of parameters for scattering off the surface */
    MTRand* my_rng;/* A random number generator for scattering off the surface */
} Surface2D;

/* A structure for holding a single ray */
typedef struct _ray2d {
    double position[2];   /* Position of the ray */
    double direction[2];  /* Direction of the ray */
    int32_t nscatters;    /* The number of scattering events the ray has undergone */
    int on_element;       /* The index of the surface element that the ray is on */
} Ray2D;

/* A structure to hold an array of structs */
typedef struct _rays2d {
    Ray2D *rays;  /* Pointer to an array of Ray structs */
    int nrays;    /* The number of rays that are stored in the struct */
} Rays2D;

#endif
