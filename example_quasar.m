%% Example: Dense second-order moment-SOS relaxation for a computer vision problem

clc; clear; close all; restoredefaultpath; % start clean

mosekpath = '../../../mosek'; % replace this with path to MOSEK in your computer
addpath(genpath(pwd))
addpath(genpath(mosekpath))

%% Create nonconvex problem data
N       = 10; % number of measurements, if N > 10, gets very slow
outrate = 0.6; % outlier rate
problem.N               = N;
problem.Covariance      = 0*eye(3); % no noise on inliers
problem.v1_distribution = 'uniform';
problem.nrOutliers      = round(N*outrate);

[a, b, R_gt, problem]   = createWahbaProblem(problem);
betasq = 0.1;

%% Define polynomial optimization problem
n = 4 + N;
q = msspoly('q',4); % unit quaternion
Rq = Rofq(q); % convert quaternion to rotation matrix
th = msspoly('th',N); % binary variables

% equality constraints
h = [q'*q - 1; th.^2 - 1];

% objective function
f = 0;
for i = 1:N
    bi = b(:,i); ai = a(:,i);
    thi = th(i);
    pi = (1+thi)/2 * (bi'*bi + ai'*ai - 2*bi'*Rq*ai) + (1-thi)/2 * betasq;
    f = f + pi;
end

pop.vars            = [q;th];
pop.objective       = f;
pop.equality        = h; 
kappa                   = 2;
[SDP,info]              = dense_sdp_relax(pop,kappa); % generate the SDP data

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
Xq = X(2:5,2:5); [V,D] = sorteig(Xq); 
feasible_q = V(:,1); feasible_q = feasible_q / norm(feasible_q);
feasible_th = sign( X(6:5+N,1) );
upper_bound = double( subs(f,[q;th],[feasible_q;feasible_th]) );
gap = abs(lower_bound - upper_bound) / (1 + abs(lower_bound) + abs(upper_bound)); % relative suboptimality gap
fprintf("Relative suboptimality gap is %3.2e.\n",gap);
