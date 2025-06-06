function out = MomConeExtremeRay(M,MomCone)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decompose a pseudomoment matrix as a convex combination of extreme rays
% M: the given pseudomoment matrix
% MomCone: sedumi format describing the pseudo-moment cone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(M,1);

weights = [];
ranks = [];
rays = {};

while true
    M = (M + M')/2;
    [V,Lam] = sorteig(M);
    lam = diag(Lam);
    r = sum(lam > 1e-6); % rank of M
    lam_pd = lam(lam >= 1e-6);
    range = V(:,1:r);
    Lam(Lam < 1e-6) = 0;
    % M = V * Lam * V'; % clean the spectrum of M

    if r == 1
        fprintf("terminate: a rank-one extreme ray.\n");
        break
    end
    
    %% minimize a random linear function on the face
    ker = V(:,r+1:end);
    Face = MomCone;
    for i = 1:size(ker,2)
        vi = ker(:,i);
        Ai = vi * vi';
        Face.At = [Face.At, sparse(Ai(:))];
        Face.b  = [Face.b; sparse(1,1)];
    end
    Face.c = randn(N^2,1);

    prob   = convert_sedumi2mosek(Face.At,...
                                  Face.b,...
                                  Face.c,...
                                  Face.K);
    [~,res] = mosekopt('minimize info',prob);
    [Xopt,~,~,~] = recover_mosek_sol_blk(res,Face.blk);
    ray = Xopt{1};

    if norm(M - ray, 'fro') < 1e-3
        fprintf("terminate: cannot be further decomposed.\n")
        break
    end

    %% subtract the ray to decrease rank
    [Vr,Lamr] = sorteig(ray);
    lamr = diag(Lamr);
    rr = sum(lamr > 1e-3); % rank of ray
    Lamr(Lamr < 1e-3) = 0;
    % ray = Vr * Lamr * Vr'; % clean the spectrum of ray
    ray = (ray + ray')/2;


    Lamhalf = diag(lam_pd .^ (-0.5));
    ray_transform = Lamhalf * (range' * ray * range) * Lamhalf;
    [~,sigs] = sorteig(ray_transform);
    step = 1 / sigs(1,1);

    weights = [weights; step];
    ranks = [ranks; rr];
    rays = [rays; {ray}];

    %% update M
    M = M - step * ray;
end

out.weights = weights;
out.ranks = ranks;
out.rays = rays;
end





