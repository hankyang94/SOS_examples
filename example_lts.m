clc; clear; close all; restoredefaultpath; % start clean

mosekpath = '../../../mosek'; % replace this with path to MOSEK in your computer
addpath(genpath(pwd))
addpath(genpath(mosekpath))

N = 12;
N_out = 5;
K = 6;
theta = [1,1];
x = randn(N,1);
y = theta(1) * x + theta(2) + 0.05*randn(N,1);
y(1:N_out) = randn(N_out,1);

figure;
scatter(x,y,50,'blue','filled');
hold on;

beta = msspoly('beta',2);
th = msspoly('th',N);

f = 0;
for i = 1:N
    fi = th(i) * (y(i) - beta(1)*x(i) - beta(2))^2;
    f = f + fi;
end

h = [sum(th) - K;
     th.^2 - th];
g = [5 - beta'*beta];

problem.vars            = [beta;th];
problem.objective       = f;
problem.equality        = h; 
problem.inequality      = g;
kappa                   = 2;
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
beta_hat = X(2:3,1);
res_sq = zeros(N,1);
for i = 1:N
    res = (y(i) - beta_hat(1)*x(i) - beta_hat(2))^2;
    res_sq(i) = res;
end
res_sorted = sort(res_sq,'ascend');
upper_bound = sum(res_sorted(1:K));

gap = abs(lower_bound - upper_bound) / (1 + abs(lower_bound) + abs(upper_bound)); % relative suboptimality gap
fprintf("Relative suboptimality gap is %3.2e.\n",gap);

xx = linspace(min(x)-0.5, max(x)+0.5);
yy = beta_hat(1) * xx + beta_hat(2);
plot(xx,yy,'LineWidth',3);
xlabel('x','FontSize',24)
ylabel('y','FontSize',24)
ax = gca;
ax.FontSize = 16;




