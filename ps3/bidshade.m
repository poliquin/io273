function [V] = bidshade(bids, maxopp, bwidth)
    % BIDSHADE Calculate amount of shading in first price auction.
    %   Given a vector of bids and a vector of maximum opponent bids, figure
    %   out how much each bid is shaded under equilibrium play in a first
    %   price auction. Uses normal kernel.
    % Input arguments:
    %   bids = N by 1 vector of bids for player i, N is number of auctions
    %   maxopp = N by 1 vector of maximum bids among bidders other than i
    %   bwidth = bandwidth parameter for kernel function
    % Outputs:
    %   V = N by 1 vector of psuedo values for player i

    T = size(bids, 1);  % number of auctions
    V = zeros(T, 1);

    for t = 1:T  % each auction-bid combo is estimated separately
        kern = makedist('Normal', 'mu', bids(t) / bwidth, 'sigma', 1);
        % find auctions in which max opponent bid is less than player i bid
        mask = maxopp < bids(t);
        mbid = bids(mask);
        % evaluate kernel function at these bids
        G = arrayfun(@(x) pdf(kern, x / bwidth), mbid, 'UniformOutput', false);
        G = sum(cell2mat(G)) / (T * bwidth);
        g = arrayfun(@(x, y) prod(pdf(kern, [x, y] ./ bwidth)), bids, ...
                     maxopp, 'UniformOutput', false);
        g = sum(cell2mat(g)) / (T * bwidth);
        % calculate amount bid is shaded and player i valuation in auction t
        shade = G / g;
        V(t) = bids(t) + shade;
    end
end
