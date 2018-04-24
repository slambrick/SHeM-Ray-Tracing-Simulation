function histRay = binMyWay(numScattersRay, maxScatter)
% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
% 
% Wrapper function for binMyWayMex.
% Produces a histogram of the provided array. Giving the number of occurances of
% each number. Goes upto the provided maximum value. Should be used only for the
% specific role it was designed for, for flexibility use 'hist' or 'histogram'
% (which one?).
%
% INPUTS:
%  numScattersRay - the list of values to bin
%  maxScatter     - value to bin upto
%
% OUTPUT:
%  histRay - histogram of the number of occurances of each number

    % Need to give binMyWayMex an array of more than 2, deal with the edge
    % cases
    if isempty(numScattersRay)
        histRay = zeros(1, maxScatter);
    elseif length(numScattersRay) == 1
        histRay = zeros(1, maxScatter);
        histRay(numScattersRay) = 1;
    elseif length(numScattersRay) == 2
        histRay = zeros(1, maxScatter);
        histRay(numScattersRay(1)) = 1;
        histRay(numScattersRay(2)) = histRay(numScattersRay(2)) + 1;
    else
        histRay = binMyWayMex(numScattersRay, maxScatter);
    end
end

