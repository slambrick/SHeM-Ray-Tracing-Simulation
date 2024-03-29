function mexCompile(recompile, library)
    if nargin == 0
        recompile = true;
    elseif nargin == 1
        library = 'none';
    elseif nargin == 2
        recompile = false;
    end
    
    if ~ispc && ~isunix
        warning('You appear to be using a system that is neither Microsoft Windows or unix, this is untested.');
    end
    
    % Create a directory for the compiled mex binaries
    % TODO: create a seperate directory for object files, make this not
    % clash with ones used by the stand alone make files
    if ~exist('bin', 'dir')
        mkdir('bin')
    end

    %% Compile the main 3D ray tracing library
    if ispc
        atom_obj = 'atom_ray_tracing_library/atom_ray_tracing3D.obj';
    else
        atom_obj = 'atom_ray_tracing_library/atom_ray_tracing3D.o';
    end
    if ~exist(atom_obj, 'file') || recompile || strcmp(library, 'atom3D')
        mex -c -R2018a CFLAGS='$CFLAGS -std=c99 -Imtwister -Iatom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
            -outdir atom_ray_tracing_library ...
            atom_ray_tracing_library/atom_ray_tracing3D.c
    end
    
    %% Compile the 2D ray tracing library
    % THIS IS BROKEN!!!!
    if ispc
        atom_obj2 = 'atom_ray_tracing_library/atom_ray_tracing2D.obj';
    else
        atom_obj2 = 'atom_ray_tracing_library/atom_ray_tracing2D.o';
    end
    if (~exist(atom_obj2, 'file') || recompile || strcmp(library, 'atom2D')) && true
        % TODO: serious errors here
        %mex -c -R2018a CFLAGS='$CFLAGS -std=c99 -Imtwister -Iatom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
        %    -outdir atom_ray_tracing_library ...
        %    atom_ray_tracing_library/atom_ray_tracing2D.c
    end
    
    %% Compile the mtwister library
    if ispc
        mtwist_obj = 'mtwister/mtwister.obj';
    else
        mtwist_obj = 'mtwister/mtwister.o';
    end
    if ~exist(mtwist_obj, 'file') || recompile || strcmp(library, 'mtwister')
        mex -c -R2018a CFLAGS='$CFLAGS -std=c99 -Imtwister -Iatom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
            -outdir mtwister ...
            mtwister/mtwister.c
    end
    
    %% For stl pinhole plate and sample
    if ispc
        tracingMex = 'bin/tracingMex.mexw64';
    else
        tracingMex = 'bin/tracingMex.mexa64';
    end
    if ~exist(tracingMex, 'file') || recompile
        if ispc
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/tracingMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.obj ...
                mtwister/mtwister.obj
        else
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/tracingMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.o ...
                mtwister/mtwister.o
        end
    end
    if ispc
        tracingGenMex = 'bin/tracingGenMex.mexw64';
    else
        tracingGenMex = 'bin/tracingGenMex.mexa64';
    end
    if ~exist(tracingGenMex, 'file') || recompile
        if ispc
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/tracingGenMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.obj ...
                mtwister/mtwister.obj 
        else
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/tracingGenMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.o ...
                mtwister/mtwister.o
        end
    end
    
    %% For a simple model of the pinhole plate
    if ispc
        multiMex = 'bin/tracingMultiMex.mexw64';
    else
        multiMex = 'bin/tracingMultiMex.mexa64';
    end
    if ~exist(multiMex, 'file') || recompile
        if ispc
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/tracingMultiMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.obj ...
                mtwister/mtwister.obj
        else
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/tracingMultiMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.o ...
                mtwister/mtwister.o
        end
    end
    if ispc
        multiGenMex = 'bin/tracingMultiGenMex.mexw64';
    else
        multiGenMex = 'bin/tracingMultiGenMex.mexa64';
    end
    if ~exist(multiGenMex, 'file') || recompile
        if ispc
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/tracingMultiGenMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.obj ...
                mtwister/mtwister.obj
        else
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/tracingMultiGenMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.o ...
                mtwister/mtwister.o
        end
    end
    
    %% For an abstract model of the pinhole plate and a sample
    if ispc
        multiGenMex = 'bin/tracingAbstractGenMex.mexw64';
    else
        multiGenMex = 'bin/tracingAbstractGenMex.mexa64';
    end
    if ~exist(multiGenMex, 'file') || recompile
        if ispc
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/tracingAbstractGenMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.obj ...
                mtwister/mtwister.obj
        else
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/tracingAbstractGenMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.o ...
                mtwister/mtwister.o
        end
    end
    
    %% For distribution or trace scattering just off a sample
    if ispc
        distCalcMex = 'bin/distributionCalcMex.mexw64';
    else
        distCalcMex = 'bin/distributionCalcMex.mexa64';
    end
    if ~exist(distCalcMex, 'file') || recompile
        if ispc
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/distributionCalcMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.obj ...
                mtwister/mtwister.obj
        else
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/distributionCalcMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.o ...
                mtwister/mtwister.o
        end
    end
    
    %% For 2D scattering calculations
    % 2D scattering seriously broken!!!!
    if ispc
        scat2dMex = 'bin/scatterRaysMex2D.mexw64';
    else
        scat2dMex = 'bin/scatterRaysMex2D.mexa64';
    end
    if (~exist(scat2dMex, 'file') || recompile) && true
        if ispc
            %mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
            %    -outdir bin ...
            %    mexFiles/scatterRaysMex2D.c ...
            %    atom_ray_tracing_library/atom_ray_tracing2D.obj ...
            %    mtwister/mtwister.obj
        else
            %mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
            %    -outdir bin ...
            %    mexFiles/scatterRaysMex2D.c ...
            %    atom_ray_tracing_library/atom_ray_tracing2D.o ...
            %    mtwister/mtwister.o
        end
    end
    
    %% For testing scattering distributions
    if ispc
        distTestMex = 'bin/distributionTestMex.mexw64';
    else
        distTestMex = 'bin/distributionTestMex.mexa64';
    end
    if ~exist(distTestMex, 'file') || recompile
        if ispc
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/distributionTestMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.obj ...
                mtwister/mtwister.obj
        else
            mex -R2018a CFLAGS='$CFLAGS -std=c99 -I mtwister -I atom_ray_tracing_library -Wall -pedantic -Wextra -O3   ' ...
                -outdir bin ...
                mexFiles/distributionTestMex.c ...
                mexFiles/extract_inputs.c ...
                atom_ray_tracing_library/atom_ray_tracing3D.o ...
                mtwister/mtwister.o
        end
    end
    
    %% A binning function I wrote.
    if ispc
        binMex = 'bin/binMyWayMex.mexw64';
    else
        binMex = 'bin/binMyWayMex.mexa64';
    end
    if ~exist(binMex, 'file') || recompile
        mex -outdir bin mexFiles/binMyWayMex.c CFLAGS='$CFLAGS -std=c99 -Wall -pedantic -Wextra -O3  '
    end
end

