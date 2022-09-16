function material = parse_scattering(model, reflectivity, sigma)
    switch model
        case 'cosine'
            material.function = 'cosine';
            material.params = 0;
            material.color = [0.8 0.8 1.0];
        case 'uniform'
            material.function = 'uniform';
            material.params = 0;
            material.color = [0.8 1.0 0.8];
        case 'specular'
            material.function = 'pure_specular';
            material.params = 0;
            material.color = [1.0 0.8 0.8];
        case 'broad_specular'
            material.function = 'broad_specular';
            material.params = [reflectivity, sigma];
            material.color = [1.0 0.4 1.0];
        case 'backscattering'
            material.function = 'backscattering';
            material.params = [reflectivity, sigma]
            material.color = [1.0 1.0 0.4];
        otherwise
            error('Input not reconised for the scattering');
    end
end

