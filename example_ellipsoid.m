clc; clear; close all;

mosekpath = '../../../mosek'; % replace this with path to MOSEK in your computer
addpath(genpath(pwd))
addpath(genpath(mosekpath))

d = 500;
n = d^2 / 8;

x = randn(d,n);

At = zeros(d^2, n);
for i = 1:n
    xixi = x(:,i) * x(:,i)';
    At(:,i) = xixi(:);
end

b = ones(n,1);
c = zeros(d^2,1);
K.s = d;
fprintf('Convert sedumi data to MOSEK.\n')
prob       = convert_sedumi2mosek(At,b,c,K);
fprintf('Done.\n')


[~,res]    = mosekopt('minimize info',prob);
blk{1,1} = 's';
blk{1,2} = d;
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,blk);



