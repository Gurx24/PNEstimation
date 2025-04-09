clc; clear; close all;

load("parameters.mat", "est_PN_code");

[reseq, poly, state] = persc(est_PN_code);