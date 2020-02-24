function mexCompile()
    % must include the GSL libraries
    % The -ffast-math improves speed

    % For stl pinhole plate and sample
    % mex -lgsl -lm -lgslcblas CFLAGS="\$CFLAGS -Wall" ...
    %     mexFiles/tracingMex.c            mexFiles/trace_ray.c ...
    %     mexFiles/tracing_functions.c ...
    %     mexFiles/intersect_detection.c          mexFiles/distributions.c ...
    %     mexFiles/ray_tracing_structs3D.c mexFiles/small_functions3D.c ...
    %     mexFiles/common_helpers.c
    % mex -lgsl -lm -lgslcblas CFLAGS="\$CFLAGS -Wall" ...
    %     mexFiles/tracingGenMex.c         mexFiles/trace_ray.c ...
    %     mexFiles/tracing_functions.c ...
    %     mexFiles/intersect_detection.c          mexFiles/distributions.c ...
    %     mexFiles/ray_tracing_structs3D.c mexFiles/small_functions3D.c ...
    %     mexFiles/common_helpers.c

    % For a simple model of the pinhole plate
    % mex -lgsl -lm -lgslcblas CFLAGS="\$CFLAGS -Wall" ...
    %     mexFiles/tracingSimpleMex.c      mexFiles/trace_ray.c ...
    %     mexFiles/tracing_functions.c ...
    %     mexFiles/intersect_detection.c          mexFiles/distributions.c ...
    %     mexFiles/ray_tracing_structs3D.c mexFiles/small_functions3D.c ...
    %     mexFiles/common_helpers.c
    % mex -lgsl -lm -lgslcblas CFLAGS="\$CFLAGS -Wall" ...
    %     mexFiles/tracingSimpleGenMex.c   mexFiles/trace_ray.c ...
    %     mexFiles/tracing_functions.c ...
    %     mexFiles/intersect_detection.c          mexFiles/distributions.c ...
    %     mexFiles/ray_tracing_structs3D.c mexFiles/small_functions3D.c ...
    %     mexFiles/common_helpers.c
    mex -g -R2018a -lgsl -lm -lgslcblas CFLAGS="\$CFLAGS -Wall" ...
        mexFiles/tracingMultiGenMex.c ...
        mexFiles/extract_inputs.c ...
        mexFiles/ray_tracing_structs3D.c ...
        mexFiles/small_functions3D.c ...
        mexFiles/common_helpers.c ...
        mexFiles/intersect_detection.c ...
        mexFiles/trace_ray.c ...
        mexFiles/tracing_functions.c ...
        mexFiles/distributions.c

    % For distribution or trac scattering just off a sample
    % mex -lgsl -lm -lgslcblas CFLAGS="\$CFLAGS -Wall" ...
    %     mexFiles/distributionCalcMex.c   mexFiles/trace_ray.c ...
    %     mexFiles/tracing_functions.c ...
    %     mexFiles/intersect_detection.c          mexFiles/distributions.c ...
    %     mexFiles/ray_tracing_structs3D.c mexFiles/small_functions3D.c ...
    %     mexFiles/common_helpers.c

    %mex -lgsl -lgslcblas x mexFiles/noPlateMex.c ...
    %    mexFiles/tracing_functions.c mexFiles/small_functions.c
    %mex -lgsl -lgslcblas CFLAGS="\$CFLAGS -ffast-math" mexFiles/cosineMex.c ...
    %    mexFiles/small_functions.c
    %mex -lgsl -lgslcblas CFLAGS="\$CFLAGS -ffast-math" mexFiles/abstractPlateMex.c ...
    %    mexFiles/tracing_functions.c mexFiles/small_functions.c

    % A binning function I wrote.
    mex mexFiles/binMyWayMex.c
end

