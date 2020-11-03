function diffuse = parse_scattering(model, reflectivity, sigma)
    switch model
        case 'cosine'
            diffuse = [1, 0];
        case 'uniform'
            diffuse = [2, 0];
        case 'specular'
            diffuse = [1 - reflectivity, 0];
        case 'broad_specular'
            diffuse = [3, sigma*pi/180];
        otherwise
            error('Input not reconised for the scattering');
    end
end

