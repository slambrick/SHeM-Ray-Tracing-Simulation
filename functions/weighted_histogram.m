% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the Sub-beam Ray Tracing Simulation, subject to the  
% GNU/GPL-3.0-or-later.
%
% weughted_histogram.m
%
% Creates a weighted histogram of the data provided. Normalises the result such
% that the integral of the histogram is 1.
%
% INPUTS:
%  values  - The values to be binned
%  weights - The weights of all the values, a vector of the same length as
%            'values'.
%  nbins   - The number of bins to put the data into.
%  normalisation - optional, the normalisation of the histogram, 'count' or
%                  'pdf', defaults to 'pdf'
%  new_fig - create a new figure? defaults to true
%
% OUTPUTS:
%  binned   - The value of each bin 
%  binWidth - The width of bin chosen.
%  binEdges - An array of the edges of the bins chosen.
function [binned, binWidth, binEdges] = weighted_histogram(values, weights, ...
        nbins, normalisation, new_fig)
    % Default value for normalisation
    if ~exist('normalisation', 'var')
        normalisation = 'pdf';
    elseif ~strcmp(normalisation, 'pdf') && ~strcmp(normalisation, 'count')
        error('Incorrect normalisation for wieghted histogram.')
    end
    
    if ~exist('new_fig', 'var')
        new_fig = false;
    end
    
    % The edges and widths of the bins
    binEdges = linspace(min(values), max(values), nbins + 1);
    binWidth = binEdges(2) - binEdges(1);
    
    % Loop through each of the bins summing the contribution to that bin from
    % all the wieghted values
    binned = zeros(1, nbins);
    for i_=1:nbins
        binned(i_) = sum(weights.*((values >= binEdges(i_)) & ...
            (values < binEdges(i_+1))));
    end
    
    if strcmp(normalisation, 'pdf')
        % Normalise such that the integral of the histogram is 1, i.e. the results
        % is normalised as a PDF
        binned = binned/sum(binned*binWidth);
    end
    
    % Produce a figure and add a plot of the data to it
    if new_fig
        figure
    end
    bar(binEdges(1:end-1) + binWidth/2, binned, 1)
end

