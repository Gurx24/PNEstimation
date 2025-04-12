clc; clear; close all;

load(".\Data\parameters.mat", "est_PN_code", "pairs");

% [num, pairs] = getoptpairs(L); % obtain all optimal pairs of L-th order
% fprintf('obtained %d optimal pairs of %d order.\n', num, L);

tic
[reseq, poly, state] = persc(est_PN_code, pairs);
toc