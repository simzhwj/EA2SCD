The corresponding paper is:

Unified robust network embedding framework for community detection via extreme adversarial attacks, Information Sciences, 2023, 643: 119200.


# ea2scd
EA2SCD

run EA2SCD.m with matlab

function [C_ea2synmf, P_ea2synmf, cost_ea2synmf] = EA2SCD(A, ncls, lambda, maxin, maxout, tol, epsilon)

%-----

% inputs:

% A: adjacent matrix

% ncls: no. of cmmunities

% parameters: maxin, maxout, tol, epsilon

%-------

Set the following parameters for Dolphins dataset, and the parameters can be fine-tuned on other datasets.

lambda = 1.2;

maxin = 1000;

maxout = 100; 

tol = 1e-2;

epsilon = 1e-5;


