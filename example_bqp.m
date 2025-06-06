%% Example: Dense second-order moment-SOS relaxation for random binary quadratic programming (BQP)

clc; clear; close all; restoredefaultpath; % start clean

mosekpath = '../../../mosek'; % replace this with path to MOSEK in your computer
addpath(genpath(pwd))
addpath(genpath(mosekpath))

%% Generate random binary quadratic program
% d       = 5; % BQP with d variables (when d=20 it is already pretty slow)
% x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
% 
% % Q       = randn(d,d); Q = Q + Q'; % a random symmetric matrix
% u = randn(d,1); v = randn(d, 1);
% Q = u * u'; % rank-one symmetric matrix
% Q = ones(d,d);
% 
% c       = randn(d,1);
% f       = x'*Q*x + c'*x; % objective function of the BQP
% h       = x.^2 - 1; % equality constraints of the BQP (binary variables, x(i)^2-1 = 0)

load("bqp_loose.mat")
f       = x'*Q*x + c'*x; % objective function of the BQP

%% Relax BQP into an SDP
% problem.vars            = x;
% problem.objective       = f;
% problem.equality        = h; 
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
feasible_sol = sign(X(2:d+1,1)); 
% feasible_sol = [1;1;-1;-1;1]; % for bqp loose
feasible_sol = [-1;1;-1;1;-1]; % for bqp loose 1
upper_bound = double(subs(f,x,feasible_sol)); % evaluate the objective at the feasible sol to get an upper bound
gap = abs(lower_bound - upper_bound) / (1 + abs(lower_bound) + abs(upper_bound)); % relative suboptimality gap
fprintf("Relative suboptimality gap is %3.2e.\n",gap);


%% find extreme rays through iterative SDP

load("MomCone.mat")

N = 21;
d = 5;
n_points = 8;
alpha = rand(n_points,1);
point = randn(d,n_points);
X = zeros(N,N);
for i = 1:n_points
    point_i = point(:,i);
    point_basis = double(subs(info.v,problem.vars,point_i));
    X = X + rand * (point_basis * point_basis');
end

[V,D] = eig(X);
r = n_points;
N = 21;
kernel = V(:,1:(N-r));
for i = 1:(N-r)
    vi = kernel(:,i);
    Ai = vi * vi';
    MomCone.At = [MomCone.At, sparse(Ai(:))];
    MomCone.b  = [MomCone.b; sparse(1,1)];
end
Cnew = randn(N,N);
MomCone.c = Cnew(:);

prob       = convert_sedumi2mosek(MomCone.At,...
                                  MomCone.b,...
                                  MomCone.c,...
                                  MomCone.K);
[~,res]    = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
Xnew = Xopt{1};
figure; bar(eig(Xnew));
norm(Xnew - X, 'fro')

[Vnew,~] = eig(Xnew);
vnew = Vnew(:,end);
V_nz = V(:,N-r+1:end);

%{
%% find extreme rays through nonlinear programming
test = [1;2;3;4;5];
v = double(subs(info.v,x,test));
r = 6;

alpha0 = randn(r,1);
X0 = -1 + 2 * rand(d, r);
z0 = [alpha0; X0(:)];
% z0 = zeros(size(z0));

% error = match_pseudomoment(z0,r,d,info.v,problem.vars,X)

myfun = @(z) match_pseudomoment(z,r,d,info.v,problem.vars,X);

options = optimoptions('fsolve',...
    'Display','iter',...
    'Algorithm','levenberg-marquardt',...
    'MaxFunctionEvaluations',1e6);

zopt = fsolve(myfun,z0,options);
alphaopt = zopt(1:r).^2;
supportopt = reshape(zopt(r+1:end),d,r);
rayopt = 0;
for i = 1:r
    ai = alphaopt(i);
    xi = supportopt(:,i);
    xi_basis = double(subs(info.v,problem.vars,xi));
    rayopt = rayopt + ai * (xi_basis*xi_basis');
end

%}



function error = match_pseudomoment(z,r,d,basis,var,target)

alpha = z(1:r);
X = reshape(z(r+1:end), d, r);

M = 0;
for i = 1:r
    ai = alpha(i)^2;
    xi = X(:,i);
    xi_basis = double(subs(basis,var,xi));
    M = M + ai * (xi_basis*xi_basis');
end

mask  = triu(true(size(M)));

diff = M - target;

error = [diff(mask); 1 - alpha'*alpha];

end




