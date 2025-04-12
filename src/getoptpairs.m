%%%%%%%%%%%%%%%%%%%%%%%%%% getoptpairs.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program realized the generation of 
% L-th order m sequence optimum pairs.
%
% date : 2025.4.2  GuRX
%
% [num, optpairs] = getoptpairs(L)
%
% ***********************************************
% L       : order of m sequence
% num     : number of optimum pairs
% optpair : all Lth-order optimum pairs
% ***********************************************

function [num, optpair] = getoptpairs(L)

num     = 0;                        % number of optimum pairs
ppoly   = gfprimfd(L, 'all');       % set of primitive polynomials
octPoly = str2double(compose("%o", bi2de(ppoly)));
C       = PEIJI(L);
idx     = nchoosek(1:length(octPoly), 2);
pairs   = [octPoly(idx(:, 1), :), octPoly(idx(:, 2), :)];

for i = 1 : size(pairs, 1)
    r = zeros(size(C, 1), 1);
    for ii = 1 : size(C, 1)
        mspair  = msgen(pairs(i, :), L);
        bmspair = 1 - 2*mspair;
        a       = bmspair(1, :);
        b       = circshift(bmspair(2, :), C{ii,1}(1));
        r(ii)   = a*b.';
    end
    if length(unique(r)) <= 3
        num              = num + 1;
        optpair(num, :) = pairs(i, :);
    end
end


function [set] = PEIJI(L)

% This function is used to construct PEIJI.
%
% date : 2025.4.2  GuRX
%
% [set] = PEIJI(L)
%
% ***********************************************
% L   : order of m sequence
% set : set of PEIJI corresponding order
% ***********************************************

P      = 2^L - 1;
C{1,1} = 0;
i      = 3;
a      = 1:P-1;
C{2,1} = 2.^(0:L-1);
a      = setdiff(a, C{2,1}, 'stable');

while ~isempty(a)
    C{i,1} = unique(mod(a(1).*C{2,1}, P));
    a      = setdiff(a, C{i,1}, 'stable');
    i      = i+1;
end

set = C;