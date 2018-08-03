function [ DEPTH, mask] = D2N( D, CamInfo)
    [h,w,~] = size(D);
    [i,j] = meshgrid(1:w,1:h);
    
    fx = CamInfo.fx;
    fy = CamInfo.fy;
    cx = CamInfo.cx;
    cy = CamInfo.cy;
    
    px = D.*(i-cx) / fx;
    py = D.*(j-cy) / fy;
    pz = D;
    P = cat(3, px, py, pz);
    
    D_im1 = [nan(h,1), D(:,1:end-1)];
    D_jm1 = [nan(1,w); D(1:end-1,:)];
    
    vx = D_jm1.*(D-D_im1)/fy;
    vy = D_im1.*(D-D_jm1)/fx;
    vz = vx.*(cx-i)/fx + vy.*(cy-j)/fy - D_im1.*D_jm1/(fx*fy);
    V = cat(3, vx, vy, vz);
    
    LEN = sqrt(sum(V.^2, 3));
    N = V./cat(3, LEN, LEN, LEN);

    threshold_angle = 0.8 * pi/2; % between normal
    mask = isfinite(N(:,:,3)) & (abs(N(:,:,3)) > cos(threshold_angle));
    
    Nimg = -N;
    Nimg(:,:,1) = Nimg(:,:,1)*0.5 + 0.5;
    Nimg(:,:,2) = Nimg(:,:,2)*0.5 + 0.5;

    
    DEPTH.D = D;
    DEPTH.P = P;
    DEPTH.V = V;
    DEPTH.LEN = LEN;
    DEPTH.N = N;
    DEPTH.Nimg = Nimg;
end

