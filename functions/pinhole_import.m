function [pinhole_surface, thePlate, aperture_abstract] = pinhole_import(pinhole_plate_inputs, sample_surface)
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    switch pinhole_plate_inputs.pinhole_model
        case 'stl'
            pinhole_surface = import_plate(pinhole_plate_inputs.plate_accuracy);

            % Plot if using a graphical window
            if ~isOctave
                if feature('ShowFigureWindows')
                    sample_surface.patchPlot(true);
                    pinhole_surface.patchPlot(false);
                    view([-5 -5 5]);
                end
            else
                sample_surface.patchPlot(true);
                pinhole_surface.patchPlot(false);
                view([-5 -5 5]);
            end

            % To pass to the functions
            thePlate = PinholeModel();
            aperture_abstract = 0;
        case 'new'
            pinhole_surface = import_newPlate(pinhole_plate_inputs.plate_accuracy);

            % Plot if using a graphical window
            if ~isOctave
                if feature('ShowFigureWindows')
                    sample_surface.patchPlot(true);
                    pinhole_surface.patchPlot(false);
                    view([-5 -5 5])
                end
            else
                sample_surface.patchPlot(true);
                pinhole_surface.patchPlot(false);
                view([-5 -5 5])
            end

            % To pass to the functions
            thePlate = PinholeModel();
            aperture_abstract = 0;
        case 'new_micro'
            % TODO
            error('Not written this bit of code yet...');
        otherwise
            % Do not have a CAD repreentation of the pinhole plate

            % Create an empty TriagSurface as the pinhole plate
            pinhole_surface = TriagSurface();

            % Struct with the information about the plate in
            thePlate = PinholeModel(pinhole_plate_inputs.plate_represent, ...
                                    pinhole_plate_inputs.n_detectors, ...
                                    pinhole_plate_inputs.circle_plate_r, ...
                                    pinhole_plate_inputs.aperture_axes, ...
                                    pinhole_plate_inputs.aperture_c);
            aperture_abstract = {pinhole_plate_inputs.aperture_theta, ...
                pinhole_plate_inputs.aperture_phi, pinhole_plate_inputs.aperture_half_cone};
    end
end