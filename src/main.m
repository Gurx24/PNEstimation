%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generating m and Gold sequence through simulation 
% tests the feasibility of the algorithm.
%
% date : 2025.3.26  GuRX
%

clc; clear; close all;

%************************ Parameters part *************************

load("..\Data\parameters.mat", "est_PN_code", "pairs");

L      = 10;                % order of PN sequence
N      = 2^L - 1;           % length of PN sequence
ber    = 0.005;             % bit error rate
mgen   = 2011;              % generator polynomial of m sequence
minist = 1001;              % initial state of m sequence
ggen   = [2011 2157];       % optimal pair of m sequence
ginist = [1001 1010];       % initial state of optimal pair

%********************** Generate PN sequence ***********************

ms     = msgen(mgen, L, minist).';          % generate m sequence
mspair = msgen(ggen, L, ginist);            % generate m sequence pair
gs     = xor(mspair(1, :), mspair(2, :)).'; % genetare Gold sequence

%********************** Obtain optimal pairs ***********************

% [num, pairs] = getoptpairs(L);
% fprintf('obtained %d optimal pairs of %d order.\n', num, L);

%********************** Generate error bits ************************

mask = rand(1, N) < ber;       % used for generating bit errors
rb   = xor(gs, mask.');        % received bits with errors

%****************** Estimation and reconstruction ******************

tic
[reseq, estpoly, eststate] = persc(rb, pairs);
toc

%************************** end of file ***************************