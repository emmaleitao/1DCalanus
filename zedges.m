function edges = zedges(wcdepth,mld)
%
% zedges
% F.Maps 2009
%
% This function computes the edges of the vertical layers
% of the model given the mixed layer depth & the water column depth.
% 
% Usage : [edges] = zedges(wcdepth,mld)
% 
% edges(timesteps,13)
% wcdepth
% mld(timesteps)

% Number of timesteps
tmax = length(mld);

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!                          !!!
%!!!   Must verify what the   !!!
%!!!   values of the first    !!!
%!!!   column should be...    !!!
%!!!                          !!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

edges = zeros(tmax,13);

% Mixed layer depth split in two
edges(:,2) = mld * 0.5;

% Lower layers equally spaced
for t = 1:tmax
    edges(t,3:13) = linspace( mld(t), wcdepth, 11 );
end

% Last layer = maximum depth
edges(:,13) = wcdepth;

end