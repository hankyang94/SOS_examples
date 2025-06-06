
clc; clear; close all; restoredefaultpath; % start clean

mosekpath = '../../../mosek'; % replace this with path to MOSEK in your computer
addpath(genpath(pwd))
addpath(genpath(mosekpath))

n = 5;
d = 2;
N = nchoosek(n+d,d);
MomCone = genMomCone(n,d);

x = msspoly('x',n);
v = [1;x;monomials(x,2:d)];

n_points = 7;
alpha = rand(n_points,1);
point = randn(n,n_points);
X = zeros(N,N);
for i = 1:n_points
    point_i = point(:,i);
    point_basis = double(subs(v,x,point_i));
    X = X + alpha(i) * (point_basis * point_basis');
end

out = MomConeExtremeRay(X,MomCone);

alpha

out.weights

