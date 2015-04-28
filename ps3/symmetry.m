function [V] = symmetry(b, bids, bwidth)
    % SYMMETRY Estimation for 1(3)(b), assuming independent, symmetric values
    %   Given a matrix of bids with bidders in columns and auctions in rows,
    %   determine how much each bid is shaded under equilibrium play in a first
    %   price auction. Uses normal kernel.
    % Input arguments:
    %   b = scalar at which to evaluate cdf
    %   bids = N by I vector of bids for I players in N auctions
    %   bwidth = bandwidth parameter for kernel function
    % Outputs:
    %   V = scalar value for cdf evaluated at b

    [T, N] = size(bids);  % number of auctions and bidders
    V = zeros(T, 1);
    kern = @(x) (1/sqrt(2*pi)) .* exp(-0.5 .* x.^2);
    
    bidstack = bids(:);
    G = mean(bidstack <= b);
    g = mean(cell2mat(arrayfun(@(x) kern((b - x) ./ bwidth), bidstack, ...
                               'UniformOutput', false))) ./ bwidth;
    V = b + G / ((N-1)*g);
end
