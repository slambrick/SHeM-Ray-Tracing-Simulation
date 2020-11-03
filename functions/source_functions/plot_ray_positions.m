function plot_ray_positions(ray_pos, thePath)
        figure
        plot(ray_pos(:,1), ray_pos(:,3), '.', 'color', [0.8 0.2 0])
        xlabel('x/mm')
        ylabel('z/mm')
        axis('equal')
        xlim([-1.5*pinhole_r 1.5*pinhole_r])
        ylim([-1.5*pinhole_r 1.5*pinhole_r])
        saveas(gcf, [thePath '/starting_positions1'], 'epsc')
        
        figure
        plot(ray_pos(:,1), ray_pos(:,2), '.', 'color', [0.8 0.2 0])
        hold on
        plot([-pinhole_r*sqrt(2), pinhole_r*sqrt(2)], [0, 0], 'color', ...
            [0.1 0.7 0.2], 'linewidth', 2)
        xlabel('x/mm')
        ylabel('y/mm')
        axis('equal')
        xlim([-1.5*pinhole_r 1.5*pinhole_r])
        ylim([-1.5*pinhole_r 1.5*pinhole_r])
        saveas(gcf, [thePath '/starting_positions2'], 'epsc')
        hold off
        
        figure
        plot(ray_pos(:,2), ray_pos(:,3), '.', 'color', [0.8 0.2 0])
        hold on
        plot([0,0], [-pinhole_r, pinhole_r], 'color', [0.1 0.7 0.2], 'linewidth', 3)
        xlabel('y/mm')
        ylabel('z/mm')
        axis('equal')
        xlim([-1.5*pinhole_r 1.5*pinhole_r])
        ylim([-1.5*pinhole_r 1.5*pinhole_r])
        saveas(gcf, [thePath '/starting_positions3'], 'epsc')
end

