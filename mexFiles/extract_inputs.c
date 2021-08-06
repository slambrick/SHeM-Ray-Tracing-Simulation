/*
 * Copyright (c) 2020, Sam Lambrick, Dan Serment.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 */
#include "extract_inputs.h"

/*
 * Take the elements from a MATLAB cell array of strings
 * and put them in a C matrix of characters.
 * Returns how many strings were extracted.
 *
 * NB: this function will not allocate space for the character matrix.
 * You should allocate space in the same function that frees it, i.e. the caller.
 */
int get_string_cell_arr(const mxArray * cell_array, char ** strings) {
    int num = mxGetN(cell_array);
    mxArray * cell;

    for(int icell = 0; icell < num; icell++) {
        cell = mxGetCell(cell_array, icell);
        if (!mxIsChar(cell))
            mexErrMsgIdAndTxt("AtomRayTracing:get_string_cell_arr:strings",
                              "Each cell must be a char array. In get_string_cell_arr.");
        if (mxGetM(cell) != 1)
            mexErrMsgIdAndTxt("AtomRayTracing:get_string_cell_arr:strings",
                              "Each cell must be a row vector. In get_string_cell_arr.");

        strings[icell] = mxArrayToString(cell);
    }
    return num;
}


/*
 * Extract the sphere C struct from the MATLAB struct array.
 * INPUTS: theSphere = mxArray containing ONE sphere struct
 *         index = the surface_index of the sphere surface
 */
AnalytSphere get_sphere(const mxArray * theSphere, int index) {
    // check if sphere is one struct
    if(!mxIsStruct(theSphere))
        mexErrMsgIdAndTxt("AtomRayTracing:get_sphere:theSphere",
                          "Must be struct array. In get_sphere.");
    if(mxGetN(theSphere) != 1 || mxGetM(theSphere) != 1)
        mexErrMsgIdAndTxt("AtomRayTracing:get_sphere:theSphere",
                          "Must be one struct. In get_sphere.");

    // get the centre
    mxArray * field = mxGetField(theSphere, 0, "c");
    if (!mxIsDouble(field) || mxGetN(field) != 3)
        mexErrMsgIdAndTxt("AtomRayTracing:get_sphere:theSphere",
                          "Centre must be array of 3 doubles. In get_sphere.");
    double * sphere_c = mxGetDoubles(field);

    // get the radius - one double
    field = mxGetField(theSphere, 0, "r");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("AtomRayTracing:get_sphere:theSphere",
                          "Sphere radius must be scalar. In get_sphere.");
    double sphere_r = mxGetScalar(field);

    // get the 'make' - one integer
    field = mxGetField(theSphere, 0, "make");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("AtomRayTracing:get_sphere:theSphere",
                          "Sphere make must be scalar. In get_sphere.");
    int make_sphere = (int)mxGetScalar(field);

    // now extract the material -- a nested struct
    mxArray * sph_material = mxGetField(theSphere, 0, "material");
    Material mat;
    char *names[] = {"sphere"};
    get_materials(sph_material, names, &mat);

    AnalytSphere sph;
    set_up_sphere(make_sphere, sphere_c, sphere_r, mat, index, &sph);
    return sph;
}

Circle get_circle(const mxArray * theCircle, int index) {
    // check if sphere is one struct
    if(!mxIsStruct(theCircle))
        mexErrMsgIdAndTxt("AtomRayTracing:get_circle:theCircle",
                          "Must be struct array. In get_circle.");
    if(mxGetN(theCircle) != 1 || mxGetM(theCircle) != 1)
        mexErrMsgIdAndTxt("AtomRayTracing:get_circle:theCircle",
                          "Must be one struct. In get_circle.");
    // get the centre
    mxArray * field = mxGetField(theCircle, 0, "c");
    if (!mxIsDouble(field) || mxGetN(field) != 3)
        mexErrMsgIdAndTxt("AtomRayTracing:get_circle:theCircle",
                          "Centre must be array of 3 doubles. In get_circle.");
    double * circle_c = mxGetDoubles(field);
    
    // get the radius - one double
    field = mxGetField(theCircle, 0, "r");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("AtomRayTracing:get_circle:theCircle",
                          "Sphere radius must be scalar. In get_circle.");
    double circle_r = mxGetScalar(field);
    
    // get the normal
    field = mxGetField(theCircle, 0, "n");
    if (!mxIsDouble(field) || mxGetN(field) != 3)
        mexErrMsgIdAndTxt("AtomRayTracing:get_circle:theCircle",
                          "Normal must be array of 3 doubles. In get_circle.");
    double * circle_n = mxGetDoubles(field);

    // get the 'make' - one integer
    field = mxGetField(theCircle, 0, "make");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("AtomRayTracing:get_circle:theCircle",
                          "Sphere make must be scalar. In get_circle.");
    int make_circle = (int)mxGetScalar(field);

    // now extract the material -- a nested struct
    mxArray * circ_material = mxGetField(theCircle, 0, "material");
    Material mat;
    char *names[] = {"circle"};
    get_materials(circ_material, names, &mat);
    
    Circle circ;
    set_up_circle(make_circle, circle_c, circle_r, circle_n, mat, index, &circ);
    return circ;
}

/* Extract array of material structs. Return how many were extracted. */
int get_materials(const mxArray * mat, char ** names, Material * target) {
    mxArray * field;
    if(!mxIsStruct(mat))
        mexErrMsgIdAndTxt("AtomRayTracing:get_materials:material",
                          "Must be struct array. In get_materials.");
    int n_materials = mxGetN(mat);

    for(int imat = 0; imat < n_materials; imat++) {
        field = mxGetField(mat, imat, "params");
        if(!mxIsDouble(field))
            mexErrMsgIdAndTxt("AtomRayTracing:get_materials:material",
                            "params must be doubles array. In get_materials.");
        int num_params = mxGetN(field);
        double * params = mxGetDoubles(field);

        field = mxGetField(mat, imat, "function");
        if(!mxIsChar(field))
            mexErrMsgIdAndTxt("AtomRayTracing:get_materials:material",
                            "function must be char array. In get_materials.");
        char * function = mxArrayToString(field);

        set_up_material(names[imat], function, params, num_params, &target[imat]);
    }
    return n_materials;
}


/**
 * Package the names, function names and parameters in the input into
 * an array of Material pointers. Return the number of materials found.
 */
int get_materials_array(const mxArray * names, const mxArray * functions,
                        const mxArray * params, Material * materials) {

    if(!mxIsCell(names) || !mxIsCell(functions) || !mxIsCell(params))
        mexErrMsgIdAndTxt("AtomRayTracing:get_materials_array:materials",
                          "Materials arrays must be cell arrays. In get_materials_array.");

    unsigned int num = mxGetN(names);
    if(mxGetN(functions) != num || mxGetN(params) != num)
        mexErrMsgIdAndTxt("AtomRayTracing:get_materials_array:materials",
                          "Dimensions of materials arrays must be the same. In get_materials_array.");

    // use cell array extraction for names and functions
    char ** names_extr = calloc(num, sizeof(char*));
    char ** functions_extr = calloc(num, sizeof(char*));

    get_string_cell_arr(names, names_extr);
    get_string_cell_arr(functions, functions_extr);

    for(unsigned int idx = 0; idx < num; idx++) {
        mxArray * params_cell = mxGetCell(params, idx);
        int num_params = mxGetN(params_cell);

        if(!mxIsDouble(params_cell))
            mexErrMsgIdAndTxt("AtomRayTracing:get_materials_array:materials",
                              "Parameters of each material must be double array. In get_materials_array.");
        double * params_extr = mxGetDoubles(params_cell);
        set_up_material(names_extr[idx], functions_extr[idx], params_extr, num_params,
        		&materials[idx]);
    }

    free(names_extr);
    free(functions_extr);
    return num;
}

/* Extract source properties from a MATLAB array and write them to the given pointers */
void get_source(const mxArray * source, int source_model, SourceParam * Source) {

    const unsigned int n_params = 7;
    if(!mxIsDouble(source))
        mexErrMsgIdAndTxt("AtomRayTracing:get_source:source",
                            "Source parameters must be double array. In get_source.");
    if(mxGetN(source) != n_params)
        mexErrMsgIdAndTxt("AtomRayTracing:get_source:source",
                          "Source must have exactly %d parameters. In get_source.", n_params);

    double * source_parameters = mxGetDoubles(source);

    Source->pinhole_r = source_parameters[0];      // pinhole radius
    Source->pinhole_c[0] = source_parameters[1];
    Source->pinhole_c[1] = source_parameters[2];
    Source->pinhole_c[2] = source_parameters[3];    // pinhole centre
    Source->theta_max = source_parameters[4];      // source max angle
    Source->init_angle = source_parameters[5];     // source initial angle
    Source->sigma = source_parameters[6];          // source std dev
    Source->source_model = source_model;
}


/*
 * Extract the plate properties from MATLAB struct of length 5
 * and write to NBackWall struct.
 *
 * INPUTS:
 * - plate_opts = mxArray containing options
 * - plate_index = the surface index of the plate
 */
NBackWall get_plate(const mxArray * plate_opts, int plate_index) {
    NBackWall plate;
    mxArray * field;

    if(!mxIsStruct(plate_opts))
        mexErrMsgIdAndTxt("AtomRayTracing:get_plate:plate_opts",
                          "Must be struct array. In get_plate.");
    if(mxGetN(plate_opts) != 1 || mxGetM(plate_opts) != 1)
        mexErrMsgIdAndTxt("AtomRayTracing:tracingMex:plate_opts",
                          "Must be one struct. In get_plate.");
        
    // Get whether to represent the pinhole plate as scattering
    field = mxGetField(plate_opts, 0, "plate_represent");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("AtomRayTracing:get_plate:plate_opts",
                          "Representation of pinhole plate must be scalar. In get_plate.");
    int plate_represent = (int)mxGetScalar(field);
    
    // Get the number of detectors
    field = mxGetField(plate_opts, 0, "n_detectors");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("AtomRayTracing:get_plate:plate_opts",
                          "Number of detectors must be scalar. In get_plate.");
    int n_detect = (int)mxGetScalar(field);
    
    // Get the radius of the pinhole plate
    field = mxGetField(plate_opts, 0, "circle_plate_r");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("AtomRayTracing:get_plate:plate_opts",
                          "Radius of pinhole plate must be scalar. In get_plate.");
    double circle_plate_r = mxGetScalar(field);
    
    // Get the axes of the apertures
    field = mxGetField(plate_opts, 0, "aperture_axes");
    if (!mxIsDouble(field) || mxGetN(field) != (unsigned int)2*n_detect)
        mexErrMsgIdAndTxt("AtomRayTracing:get_plate:plate_opts",
                          "Aperture axes must be array of 2xnumber of detectors doubles. In get_plate.");
    double * aperture_axes = mxGetDoubles(field);
    
    // Get the centres of the apertures
    field = mxGetField(plate_opts, 0, "aperture_c");
    if (!mxIsDouble(field) || mxGetN(field) !=  (unsigned int)2*n_detect)
        mexErrMsgIdAndTxt("AtomRayTracing:get_plate:plate_opts",
                          "Aperture centres must be array of 2xnumber of detectors doubles. In get_plate.");
    double * aperture_c = mxGetDoubles(field);
        
    // now extract the material -- a nested struct
    mxArray * circ_material = mxGetField(plate_opts, 0, "material");
    Material mat;
    char *names[] = {"plate"};
    get_materials(circ_material, names, &mat);

    set_up_plate(plate_represent, n_detect, circle_plate_r, aperture_axes, aperture_c,
        mat, plate_index, &plate);
    
    return plate;
}

/*
 * Extract the plate properties from MATLAB struct of length 5
 * and write to NBackWall struct.
 *
 * INPUTS:
 * - plate_opts = mxArray containing options
 * - plate_index = the surface index of the plate
 */
AbstractHemi get_abstract(const mxArray * plate_opts, int plate_index) {
    AbstractHemi plate;
    mxArray * field;

    if(!mxIsStruct(plate_opts))
        mexErrMsgIdAndTxt("AtomRayTracing:get_abstract:plate_opts",
                          "Must be struct array. In get_abstract.");
    if(mxGetN(plate_opts) != 1 || mxGetM(plate_opts) != 1)
        mexErrMsgIdAndTxt("AtomRayTracing:tracingMex:plate_opts",
                          "Must be one struct. In get_abstract.");

    
    // Get the half cone angle
    field = mxGetField(plate_opts, 0, "half_cone_angle");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("AtomRayTracing:get_abstract:plate_opts",
                          "Half cone angle must be scalar. In get_abstract.");
    double half_cone_angle = mxGetScalar(field);
    
    // Get the axes of the apertures
    field = mxGetField(plate_opts, 0, "direction");
    if (!mxIsDouble(field) || mxGetN(field) != 3)
        mexErrMsgIdAndTxt("AtomRayTracing:get_abstract:plate_opts",
                          "Aperture axes must be array of 3. In get_abstract.");
    double * det_dir = mxGetDoubles(field);

    plate.half_cone_angle = half_cone_angle;
    for (int i = 0; i < 3; i++)
        plate.det_dir[i] = det_dir[i];
    plate.surf_index = plate_index;
    
    return plate;
}
