%% Example: Dense SDP relaxation for random binary quadratic programming (BQP)

clc; clear; close all; restoredefaultpath; % start clean

mosekpath = '../../../mosek';
addpath(genpath(pwd))
addpath(genpath(mosekpath))

%% Generate random binary quadratic program
d       = 20; % BQP with d variables
x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
Q       = randn(d,d); Q = Q + Q'; % a random symmetric matrix
c       = randn(d,1);
f       = x'*Q*x + c'*x; % objective function of the BQP
h       = x.^2 - 1; % equality constraints of the BQP (binary variables, x(i)^2-1 = 0)

%% Relax BQP into an SDP
problem.vars            = x;
problem.objective       = f;
problem.equality        = h; 
kappa                   = 2; % relaxation order
[SDP,info]              = dense_sdp_relax(problem,kappa);

%% Solve using MOSEK, should be slow
prob       = convert_sedumi2mosek(SDP.sedumi.At,...
                                  SDP.sedumi.b,...
                                  SDP.sedumi.c,...
                                  SDP.sedumi.K);
[~,res]    = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);

%% Compute certificate of global optimality



