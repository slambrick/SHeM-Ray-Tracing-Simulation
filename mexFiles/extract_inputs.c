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
            mexErrMsgIdAndTxt("MyToolbox:tracingMex:strings",
                              "Each cell must be a char array.");
        if (mxGetM(cell) != 1)
            mexErrMsgIdAndTxt("MyToolbox:tracingMex:strings",
                              "Each cell must be a row vector.");

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
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:theSphere",
                          "Must be struct array");
    if(mxGetN(theSphere) != 1 || mxGetM(theSphere) != 1)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:theSphere",
                          "Must be one struct");

    // get the centre
    mxArray * field = mxGetField(theSphere, 0, "c");
    if (!mxIsDouble(field) || mxGetN(field) != 3)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:theSphere",
                          "Centre must be array of 3 doubles.");
    double * sphere_c = mxGetDoubles(field);

    // get the radius - one double
    field = mxGetField(theSphere, 0, "r");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("MuToolbox:tracingMex:theSphere",
                          "Sphere radius must be scalar");
    double sphere_r = mxGetScalar(field);

    // get the 'make' - one integer
    field = mxGetField(theSphere, 0, "make");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("MuToolbox:tracingMex:theSphere",
                          "Sphere make must be scalar");
    int make_sphere = (int)mxGetScalar(field);

    // now extract the material -- a nested struct
    mxArray * sph_material = mxGetField(theSphere, 0, "material");
    if(!mxIsStruct(sph_material))
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:sphereMaterial",
                          "Must be struct array");
    if(mxGetN(sph_material) != 1 || mxGetM(sph_material) != 1)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:sphereMaterial",
                          "Must be one struct");
    int num_params = mxGetN(mxGetField(sph_material, 0, "params"));
    double * params = mxGetDoubles(mxGetField(sph_material, 0, "params"));
    char * function = mxArrayToString(mxGetField(sph_material, 0, "function"));

    Material mat = set_up_material("sphere", function, params, num_params);

    AnalytSphere sph = set_up_sphere(make_sphere, sphere_c, sphere_r, mat, index);
    return sph;
}


/**
 * Package the names, function names and parameters in the input into
 * an array of Material pointers. Return the number of materials found.
 */
int get_materials(const mxArray * names, const mxArray * functions,
                  const mxArray * params, Material * materials) {

    if(!mxIsCell(names) || !mxIsCell(functions) || !mxIsCell(params))
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:materials",
                          "Materials arrays must be cell arrays.");

    int num = mxGetN(names);
    if(mxGetN(functions) != num || mxGetN(params) != num)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:materials",
                          "Dimensions of materials arrays must be the same.");

    // use cell array extraction for names and functions
    char ** names_extr = mxCalloc(num, sizeof(char*));
    char ** functions_extr = mxCalloc(num, sizeof(char*));

    get_string_cell_arr(names, names_extr);
    get_string_cell_arr(functions, functions_extr);

    for(int idx = 0; idx < num; idx++) {
        mxArray * params_cell = mxGetCell(params, idx);
        int num_params = mxGetN(params_cell);

        if(!mxIsDouble(params_cell))
            mexErrMsgIdAndTxt("MyToolbox:tracingMex:materials",
                              "Parameters of each material must be double array.");
        double * params_extr = mxGetDoubles(params_cell);
        materials[idx] = set_up_material(
            names_extr[idx], functions_extr[idx], params_extr, num_params);
    }

    mxFree(names_extr); mxFree(functions_extr);
    return num;
}

/* Extract source properties from a MATLAB array and write them to the given pointers */
void get_source(const mxArray * source, double * pinhole_r, double * pinhole_c,
                double * theta_max, double * init_angle, double * sigma) {
    double * source_parameters = mxGetDoubles(source);

    *pinhole_r = source_parameters[0];      // pinhole radius
    pinhole_c[0] = source_parameters[1];
    pinhole_c[1] = source_parameters[2];    // pinhole centre
    *theta_max = source_parameters[3];      // source max angle
    *init_angle = source_parameters[4];     // source initial angle
    *sigma = source_parameters[5];          // source std dev
}


/*
 * Extract the plate properties from MATLAB cell array of length 5
 * and write to NBackWall struct.
 *
 * INPUTS:
 * - plate_opts = mxArray containing options
 * - plate_index = the surface index of the plate
 */
NBackWall get_plate(const mxArray * plate_opts, int plate_index) {
    NBackWall plate;

    if(!mxIsCell(plate_opts))
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:plate",
                          "Must be cell array");
    if(mxGetN(plate_opts) != 5 && mxGetM(plate_opts) != 5)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:plate",
                          "Must be of length 5");

    plate.plate_represent = (int)mxGetScalar(mxGetCell(plate_opts, 0));
    plate.n_detect = (int)mxGetScalar(mxGetCell(plate_opts, 1));
    plate.circle_plate_r = mxGetScalar(mxGetCell(plate_opts, 2));
    plate.aperture_axes = mxGetDoubles(mxGetCell(plate_opts, 3));
    plate.aperture_c = mxGetDoubles(mxGetCell(plate_opts, 4));
    plate.surf_index = plate_index;

    return plate;
}
