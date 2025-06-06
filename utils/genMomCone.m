function MomCone = genMomCone(n,d)

x       = msspoly('x',n); % symbolic decision variables using SPOTLESS
f = zeros(1,n) * x;
problem.vars            = x;
problem.objective       = f;
[SDP,~]                 = dense_sdp_relax(problem,d);

MomCone     = SDP.sedumi;
MomCone.blk = SDP.blk;

end
