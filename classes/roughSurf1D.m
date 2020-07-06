classdef roughSurf1D
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Grouping class
    end
    
    methods
        function [f,x,Z] = rsgeng(N,rL,h,cl)
            %
            % [f,x] = rsgeng1D(N,rL,h,cl) 
            %
            % generates a 1-dimensional random rough surface f(x) with N surface points. 
            % The surface has a Gaussian height distribution function and a Gaussian 
            % autocovariance function, where rL is the length of the surface, h is the 
            % RMS height and cl is the correlation length.
            %
            % Input:    N   - number of surface points
            %           rL  - length of surface
            %           h   - rms height
            %           cl  - correlation length
            %
            % Output:   f   - surface heights
            %           x   - surface points
            %
            % Last updated: 2010-07-26 (David Bergström).  
            %

            format long;

            x = linspace(-rL/2,rL/2,N);

            Z = h.*randn(1,N); % uncorrelated Gaussian random rough surface distribution
                                 % with mean 0 and standard deviation h

            % Gaussian filter
            F = exp(-x.^2/(cl^2/2));

            % correlation of surface using convolution (faltung), inverse
            % Fourier transform and normalizing prefactors
            f = sqrt(2/sqrt(pi))*sqrt(rL/N/cl)*ifft(fft(Z).*fft(F));
        end
        
        function [f,x] = rsgene(N,rL,h,cl)
            %
            % [f,x] = rsgene1D(N,rL,h,cl) 
            %
            % generates a 1-dimensional random rough surface f(x) with N surface points.
            % The surface has a Gaussian height distribution and an
            % exponential autocovariance function, where rL is the length of the surface, 
            % h is the RMS height and cl is the correlation length.
            %
            % Input:    N   - number of surface points
            %           rL  - length of surface
            %           h   - rms height
            %           cl  - correlation length
            %
            % Output:   f   - surface heights
            %           x   - surface points
            %
            % Last updated: 2010-07-26 (David Bergström).  
            %

            format long;

            x = linspace(-rL/2,rL/2,N);

            Z = h.*randn(1,N); % uncorrelated Gaussian random rough surface distribution
                               % with rms height h

            % Exponential filter
            % NOTE: I think this means that the correlation function remains unnormalized 
            F = exp(-abs(x)/cl);

            % correlated surface generation including convolution (faltning) and inverse
            % Fourier transform and normalizing prefactors
            f = ifft(fft(Z).*fft(F));
            f2 = sqrt(2)*sqrt(rL/N/cl)*ifft(fft(Z).*fft(F));
        end
        
        function [hdf,bc,h] = hdf(f,b,opt)
            %
            % [hdf,bc,h] = hdf1D(f,b,opt)
            %
            % calculates the height distribution function and rms height of a 
            % 1-d surface profile f(x) using b bins.
            %
            % Input:    f   - surface heights
            %           b   - number of bins
            %           opt - optional parameter (type 'hist' for drawing histogram,
            %                 type 'plot' for continuous plot) 
            %
            % Output:   hdf - height distribution function
            %           bc  - bin centers
            %           h   - rms height
            %
            % Last updated: 2009-03-11 (David Bergström)
            %

            format long

            h = std(f); % rms height

            [hdf,bc] = hist(f,b); 
            hdf = hdf/sum(hdf); % normalization to get height distribution function

            % optional plotting
            if nargin<3 || isempty(opt)
                return;
            end
            if nargin==3
                if ischar(opt)
                    if strcmp(opt,'hist')
                        bar(bc,hdf);
                        xlabel('surface height')
                        ylabel('probability')
                        title('Histogram of height distribution function for height of surface profile y=f(x)')
                    elseif strcmp(opt,'plot');
                        plot(bc,hdf);
                        xlabel('surface height')
                        ylabel('probability')
                        title('Plot of height distribution function for height of surface profile y=f(x)')
                    else
                        fprintf('%s is not a valid option. Type \''help pmf1D\'' for further details.\n',opt);
                    end
                else
                    fprintf('Option must be a string. Type \''help pmf1D\'' for further details.\n');
                end
            end
        end
        
        function [acf,cl,lags] = acf(f,x,opt)
            %
            % [acf,cl,lags] = acf1D(f,x)
            %
            % calculates the autocovariance function and correlation length of a 
            % 1-d surface profile f(x).
            %
            % Input:    f    - surface heights
            %           x    - surface points
            %           opt - optional parameter (type 'plot' for plotting the 
            %                 normalized autocovariance function)
            %
            % Output:   acf  - autocovariance function
            %           cl   - correlation length
            %           lags - lag length vector (useful for plotting the acf)
            %
            % Last updated: 2010-07-26 (David Bergstrom)
            %

            format long

            N = length(x); % number of sample points
            lags = linspace(0,x(N)-x(1),N); % lag lengths

            % autocovariance function calculation
            c=xcov(f,'coeff'); % the autocovariance function
            acf=c(N:2*N-1); % right-sided version

            % correlation length calculation
            k = 1;
            while (acf(k) > 1/exp(1))
                k = k + 1;
            end
            cl = 1/2*(x(k-1)+x(k)-2*x(1)); % the correlation length

            % optional plotting
            if nargin<3 || isempty(opt)
                return;
            end
            if nargin==3 
                if ischar(opt)
                    if strcmp(opt,'plot')
                        plot(lags,acf);
                        xlabel('lag length')
                        ylabel('acf (normalized)')
                        title('Plot of the normalized acf')
                    else
                        fprintf('%s is not a valid option. Type \''help acf1D\'' for further details.\n',opt);
                    end
                else
                    fprintf('Option must be a string. Type \''help acf1D\'' for further details.\n');
                end
            end
        end
    end
end

