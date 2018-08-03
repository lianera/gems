function Dout = refine(I, D, CamInfo, weights)    
    [h,w] = size(D);
    
    %% pre process
    D = smoothdepth(D, 7, 0.3);
    Ig = rgb2gray(I);

    %% points, normal, mask
    [DEPTH, mask] = D2N(D, CamInfo);
    nhood = ones(3,3);
    mask(1,:) = 0;
    mask(end,:) = 0;
    mask(:,1) = 0;
    mask(:,end) = 0;
    mask = imerode(mask, nhood);

    %% estimate hamonics
    VN = reshape(DEPTH.N, h*w, 3);
    A = harmonics(VN);
    MaskA = all(isfinite(A), 2);
    A(~MaskA, :) = [];
    VIg = reshape(Ig, h*w, 1);
    VIg(~MaskA, :) = [];
    l = A \ VIg;

    S = nan(h, w);
    S(MaskA) = A*l;
    R = Ig./S;

    Dout = optimize( Ig, D, R, l, CamInfo, mask, weights);
end