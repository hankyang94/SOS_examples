%% Example: Dense second-order moment-SOS relaxation for random quartic optimization over the sphere
clc; clear; close all; restoredefaultpath; % start clean

mosekpath = '../../../mosek'; % replace this with path to MOSEK in your computer
addpath(genpath(pwd))
addpath(genpath(mosekpath))

%% Generate random quartic optimization over the unit sphere
d       = 10; % d variables
x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
x4      = monomials(x,4); % list of all monomials up to degree 4
c       = randn(length(x4),1);
f       = c'*x4; % objective function: a random degree-4 polynomial
h       = x'*x - 1; % equality constraints (unit sphere, x'*x - 1 = 0)

%% Relax BQP into an SDP
problem.vars            = x;
problem.objective       = f;
problem.equality        = h; 
kappa                   = 2; % relaxation order (Gomans-Williamson corresponds to kappa=1)
[SDP,info]              = dense_sdp_relax(problem,kappa); % generate the SDP data

%% Solve using MOSEK
prob       = convert_sedumi2mosek(SDP.sedumi.At,...
                                  SDP.sedumi.b,...
                                  SDP.sedumi.c,...
                                  SDP.sedumi.K);
[~,res]    = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
X = Xopt{1};

%% Compute certificate of global optimality
lower_bound = obj(1); % SDP relaxation provides a lower bound
feasible_sol = X(2:d+1,1); feasible_sol = feasible_sol / norm(feasible_sol); % Simply take the order-one monomials and round it to be feasible
upper_bound = double(subs(f,x,feasible_sol)); % evaluate the objective at the feasible sol to get an upper bound
gap = abs(lower_bound - upper_bound) / (1 + abs(lower_bound) + abs(upper_bound)); % relative suboptimality gap
fprintf("Relative suboptimality gap is %3.2e.\n",gap);
