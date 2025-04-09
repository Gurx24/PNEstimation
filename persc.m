%%%%%%%%%%%%%%%%%%%%%%%%%%%%% persc.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameters estimation and reconstruction of spreading code.
%
% This program implements the identification of 
% the spreading code type, generator polynomial, 
% and initial state,as well as 
% the reconstruction of the spreading code.
%
% date : 2025.3.26  GuRX
%
% [reseq, poly, state] = pessc(seq)
%
% *************************************************************
% seq     : estimated spreading code sequence
% reseq   : reconstructed sequence without error bits
% poly    : estimated generator polynomial of speading sequence
% state   : estimated initial state of speading sequence
% *************************************************************

function [reseq, poly, state] = persc(seq)

% load("parameters.mat", "pairs");

N     = length(seq);           % length of ss sequence
L     = log2(N+1);             % number of LFSR stages
a     = 2 * seq - 1;           % bipolar Encoding
Tm    = 0.4*N;                 % threshold of m sequence
Tg    = 0.4*N;                 % threshold of Gold sequence
ppoly = gfprimfd(L, 'all');    % set of Lth-order m sequence polynomials
rs    = [zeros(1, L-1) 1];     % random state
flag  = false;                 % execution flag of the algorithm

% Traverse local m sequence generator polynomial set
for i = 1:size(ppoly, 1)
    tempseq = generator(ppoly(i, :), rs);
    b       = 2*tempseq - 1;
    R       = zeros(1, N);
    for t = 1 : N
        R(t) = a.' * circshift(b, t-1);
    end
    [mv, idx] = max(R);
    if abs(mv) > Tm
        disp("m sequence");
        poly  = ppoly(i, :);
        state = [];
        reseq = circshift(tempseq, idx-1);
        flag  = true;
        break;
    end
end

% Traverse local Gold sequence generator polynomial set
if ~flag
    disp('non m sequence');
    [num, pairs] = getoptpairs(L); % obtain all optimal pairs of L-th order
    fprintf('obtained %d optimal pairs of %d order.\n', num, L);
    for j = 1:size(pairs, 1)
        mpair    = dec2bin(base2dec(num2str(pairs(j, :).'), 8))-'0';
        gpoly    = conv(mpair(1, :), mpair(2, :));
        eststate = GoldStateEst(seq, gpoly);
        goldgen  = generator(gpoly, eststate);
        goldgen  = goldgen(1:N);
        r = (2*seq-1).' * (2*goldgen-1);
        if abs(r) >= Tg
            disp("Gold sequence");
            poly  = mpair;
            state = eststate;
            reseq = goldgen;
            flag  = true;
            break;
        end
    end
end

if flag == false
    disp('recognition and reconstruction failed');
end