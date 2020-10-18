/*
 * Copyright (c) 2020, Dan Seremet.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 * Functions to extract data from MATLAB mxArray structures into
 * the C structures for ray tracing.
 */

#ifndef _extract_matlab_inputs
#define _extract_matlab_inputs

#include <mex.h>

#include "ray_tracing_structs3D.h"

/*
 * Take the elements from a MATLAB cell array of strings
 * and put them in a C matrix of characters.
 * Returns how many strings were extracted.
 *
 * NB: this function will not allocate space for the character matrix.
 * You should allocate space in the same function that frees it, i.e. the caller.
 */
int get_string_cell_arr(const mxArray * cell_array, char ** strings);

/*
 * Extract the sphere C struct from the MATLAB struct array.
 * INPUTS: theSphere = mxArray containing ONE sphere struct
 *         index = the surface_index of the sphere surface
 */
AnalytSphere get_sphere(const mxArray * theSphere, int index);

/* Extract array of material structs. Return how many were extracted. */
int get_materials(const mxArray * mat, char ** names, Material * target);

/**
 * Package the names, function names and parameters in the input into
 * an array of Material pointers. Return the number of materials found.
 */
int get_materials_array(const mxArray * names, const mxArray * functions,
                        const mxArray * params, Material * materials);

/* Extract source properties from a MATLAB array and write them to the given pointers */
void get_source(const mxArray * source, int source_model, SourceParam * Source);

/*
 * Extract the plate properties from MATLAB cell array of length 5
 * and write to NBackWall struct.
 *
 * INPUTS:
 * - plate_opts = mxArray containing options
 * - plate_index = the surface index of the plate
 */
NBackWall get_plate(const mxArray * plate_opts, int plate_index);

#endif
