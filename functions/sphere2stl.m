function sphere2stl(sphere_r, sphere_c, dist_to_sample)
    sphere_disp = -sphere_c(2) - dist_to_sample;
    disp(sphere_disp)
    x = linspace(-sphere_r*2, sphere_r*2, 100);
    z = linspace(-sphere_r*2, sphere_r*2, 100);
    
    [xx, zz] = meshgrid(x, z);
    yy = sphere_eval(xx, zz);
    
    figure;
    surf(xx, zz, yy);
    axis('equal')
    
    surf2stl('test.stl', xx, zz, yy);
    
    
    
    function y = sphere_eval(x, z)  
        ind = x.^2 + z.^2 < sphere_r^2 - sphere_disp^2;
        y = -dist_to_sample*ones(size(x));
        y(ind) = sqrt(sphere_r^2 - x(ind).^2 - z(ind).^2) + sphere_c(2);
    end
end


