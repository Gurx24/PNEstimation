clc; clear; close all;

load("parameters.mat", "est_PN_code");

tic
[reseq, poly, state] = persc(est_PN_code);
toc