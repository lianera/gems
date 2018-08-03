function [ D ] = optimize( I, Dinit, k, l, CamInfo, mask, weights)
    [h, w] = size(Dinit);
    k = ones(h,w);
    D = Dinit;    
    
    for i = 1:10
        [DEPTH, ~] = D2N( D, CamInfo);
        J = jaco( DEPTH, k, l, CamInfo, mask, weights);
        F = Fterms( DEPTH, k, l, I, Dinit, mask, weights); 
        E = sum(F.*F);
        disp(E);
        A = J'*J;
        b = -J'*F;
        [mh, mw] = size(mask);  
        pixel_num = mh*mw;
        tolerance = sqrt(pixel_num*0.001^2); % norm(b-A*x0)
        delta = pcg(A, b, tolerance, 5);    
        delta = reshape(delta, h, w);
        delta(~mask) = 0;
        D = D+delta;
    end
end

function [ F ] = Fterms( DEPTH, k, l, I, Dinit, mask, weights)
    ws = 0.25;
    
    D = DEPTH.D;
    N = DEPTH.N;
    P = DEPTH.P;

    [h, w] = size(D);
    VN = reshape(N, h*w, 3);
    A = harmonics(VN);
    VK = reshape(k, h*w, 1);
    B = VK.*(A*l);
    B = reshape(B, h, w);
    
    B_ip1 = MatOffset(B, 1, 0);
    I_ip1 = MatOffset(I, 1, 0);
    
    r1n = B - B_ip1 - (I - I_ip1);
    r1n = reshape(r1n, h*w, 1)*weights.wg;
    
    B_jp1 = MatOffset(B, 0, 1);
    I_jp1 = MatOffset(I, 0, 1);
    r2n = B - B_jp1 - (I - I_jp1);
    r2n = reshape(r2n, h*w, 1)*weights.wg;
    
    P_ip1 = MatOffset(P, 1, 0);
    P_jp1 = MatOffset(P, 0, 1);
    P_im1 = MatOffset(P, -1, 0);
    P_jm1 = MatOffset(P, 0, -1);
    Es = P - ws*(P_ip1+P_jp1+P_im1+P_jm1);
    r3n = reshape(Es(:,:,1), h*w, 1)*weights.ws;
    r4n = reshape(Es(:,:,2), h*w, 1)*weights.ws;
    r5n = reshape(Es(:,:,3), h*w, 1)*weights.ws;
    
    r6n = D - Dinit;
    r6n = reshape(r6n, h*w, 1)*weights.wp;
    
    F = [r1n;r2n;r3n;r4n;r5n;r6n];
    
    vmask = reshape(mask, h*w, 1);
    rm = repmat(~vmask, [6, 1]);
    F(rm) = 0;
    F = sparse(F);
end

function [ J ] = jaco( DEPTH, k, l, CamInfo, mask, weights)
    D = DEPTH.D;
    [h,w] = size(D);
    
    [gr1, gr2] = gradient_r1r2(DEPTH, k, l, CamInfo);
    [gr3, gr4, gr5] = gradient_r3r4r5(D, CamInfo);
    gr6 = gradient_r6(D);
    
    grs = [gr1,gr2, gr3, gr4, gr5, gr6];
    params = [weights.wg, weights.wg, ...
        weights.ws, weights.ws, weights.ws, ...
        weights.wp];
    lines = 0;
    for gr = grs
        [~, c] = size(gr.offset);
        lines = lines + c;
    end
        
    vmask = mask(:);

    [~,rnum] = size(grs);
    n = h*w;            % col num  
    m = n*rnum;            % row num
    enum = sum(vmask);    % valid element num
    
    rows_idx = nan(enum, lines);
    cols_idx = nan(enum, lines);
    vals = nan(enum, lines);
    line = 0;
    for k = 1:rnum        
        gr = grs(k);
        [~, c] = size(gr.offset);
        for i = 1:c;
            line = line+1;
            offset = gr.offset(i);
            fullrows = (k-1)*n+1 : k*n;
            fullcols = (1:n) + offset;
            fullval = gr.vals{i};
            fullval = fullval(:)*params(k);
            
            rows_idx(:,line) = fullrows(vmask);
            cols_idx(:,line) = fullcols(vmask);
            vals(:, line) = fullval(vmask);
        end
    end
    
    J = sparse(rows_idx(:), cols_idx(:), vals(:), m, n);
end

function [gr1, gr2] = gradient_r1r2(DEPTH, k, l, CamInfo)
    D = DEPTH.D;
    [h,w] = size(D);
    
    [dv, dv_dD_im1, dv_dD_jm1, dv_ip1, dv_jp1] = dv_dD(D, CamInfo);

    dB = dB_dD(k, l, dv, DEPTH);
    
    k_ip1 = MatOffset(k, 1, 0);
    k_jp1 = MatOffset(k, 0, 1);
    DEPTH_ip1 = DepthOffset(DEPTH, 1, 0);
    DEPTH_jp1 = DepthOffset(DEPTH, 0, 1);

    dB_ip1 = dB_dD(k_ip1, l, dv_ip1, DEPTH_ip1);
    dB_jp1 = dB_dD(k_jp1, l, dv_jp1, DEPTH_jp1);
    
    dr1 = dB-dB_ip1;
    dr2 = dB-dB_jp1;
    
    dB_dD_im1 = dB_dD(k, l, dv_dD_im1, DEPTH);
    dB_dD_jm1 = dB_dD(k, l, dv_dD_jm1, DEPTH);

    dr1_dD_im1 = dB_dD_im1;
    dr1_dD_jm1 = dB_dD_jm1;
    dr2_dD_im1 = dB_dD_im1;
    dr2_dD_jm1 = dB_dD_jm1;    
    
    
    gr1.vals = {dr1, dr1_dD_im1, dr1_dD_jm1};
    gr2.vals = {dr2, dr2_dD_im1, dr2_dD_jm1};
    gr1.offset = [0, -h, -1];
    gr2.offset = [0, -h, -1];    
end

function [gr3, gr4, gr5] = gradient_r3r4r5(D, CamInfo)
    [h,w] = size(D);
    
    ws = 0.25;
    
    fx = CamInfo.fx;
    fy = CamInfo.fy;
    cx = CamInfo.cx;
    cy = CamInfo.cy;
    
    [i,j] = meshgrid(1:w,1:h);
    
    dr3 = (i-cx)/fx;
    dr3_dD_ip1 = -ws*(i+1-cx)/fx;
    dr3_dD_im1 = -ws*(i-1-cx)/fx;
    dr3_dD_jp1 = -ws*(i-cx)/fx;
    dr3_dD_jm1 = -ws*(i-cx)/fx;
    
    dr4 = (j-cy)/fy;
    dr4_dD_ip1 = -ws*(j-cy)/fy;
    dr4_dD_im1 = -ws*(j-cy)/fy;
    dr4_dD_jp1 = -ws*(j+1-cy)/fy;
    dr4_dD_jm1 = -ws*(j-1-cy)/fy;
    
    dr5 = ones(h,w);
    dr5_dD_ip1 = -ws*ones(h,w);
    dr5_dD_im1 = -ws*ones(h,w);
    dr5_dD_jp1 = -ws*ones(h,w);
    dr5_dD_jm1 = -ws*ones(h,w);

    gr3.vals = {dr3, dr3_dD_ip1, dr3_dD_im1, dr3_dD_jp1, dr3_dD_jm1};
    gr4.vals = {dr4, dr4_dD_ip1, dr4_dD_im1, dr4_dD_jp1, dr4_dD_jm1};
    gr5.vals = {dr5, dr5_dD_ip1, dr5_dD_im1, dr5_dD_jp1, dr5_dD_jm1};
    gr3.offset = [0, h, -h, 1, -1];
    gr4.offset = [0, h, -h, 1, -1];
    gr5.offset = [0, h, -h, 1, -1];

end

function gr6 = gradient_r6(D)
    [h, w] = size(D);
    dr6 = ones(h,w);
    gr6.vals = {dr6};
    gr6.offset = 0;
end

function [dv, dv_dD_im1, dv_dD_jm1, dv_ip1, dv_jp1] = dv_dD(D, CamInfo)
    [h,w] = size(D);
    [i,j] = meshgrid(1:w,1:h);
    
    fx = CamInfo.fx;
    fy = CamInfo.fy;
    cx = CamInfo.cx;
    cy = CamInfo.cy;    
    
    D_im1 = MatOffset(D, -1, 0);
    D_jm1 = MatOffset(D, 0, -1);
    D_ip1 = MatOffset(D, 1, 0);
    D_jp1 = MatOffset(D, 0, 1);
    D_ip1_jm1 = MatOffset(D, 1, -1);
    D_im1_jp1 = MatOffset(D, -1, 1);

    dvx = D_jm1/fy;
    dvy = D_im1/fx;
    dvz = (cx-i)/fx.*dvx + (cy-j)/fy.*dvy;  
    dv = cat(3, dvx, dvy, dvz);
    
    dvx_dD_im1 = -D_jm1/fy;    
    dvy_dD_im1 = (D - D_jm1)/fx;
    dvz_dD_im1 = (D_jm1.*(i+j-cx-cy-1)+D.*(cy-j))/(fx*fy);
    dv_dD_im1 = cat(3, dvx_dD_im1, dvy_dD_im1, dvz_dD_im1);

    dvx_dD_jm1 = (D - D_im1)/fy;
    dvy_dD_jm1 = -D_im1/fx;
    dvz_dD_jm1 = (D_im1.*(i+j-cx-cy-1)+D.*(cx-i))/(fx*fy);
    dv_dD_jm1 = cat(3, dvx_dD_jm1, dvy_dD_jm1, dvz_dD_jm1);

    dvx_ip1 = -D_ip1_jm1/fy;
    dvy_ip1 = (D_ip1-D_ip1_jm1)/fx;
    dvz_ip1 = (D_ip1_jm1.*(i+j-cx-cy)+D_ip1.*(cy-j))/(fx*fy);
    dv_ip1 = cat(3, dvx_ip1, dvy_ip1, dvz_ip1);

    dvx_jp1 = (D_jp1 - D_im1_jp1) / fy;
    dvy_jp1 = -D_im1_jp1 / fx;
    dvz_jp1 = (D_im1_jp1.*(i+j-cx-cy)+D_jp1.*(cx-i))/(fx*fy);
    dv_jp1 = cat(3, dvx_jp1, dvy_jp1, dvz_jp1);
end

function dB = dB_dD(k, l, dv, DEPTH)
    N = DEPTH.N;
    V = DEPTH.V;
    nx = N(:,:,1);
    ny = N(:,:,2);
    nz = N(:,:,3);
    vx = V(:,:,1);
    vy = V(:,:,2);
    vz = V(:,:,3);
    
    lencb = DEPTH.LEN.^3;
    
    dvx = dv(:,:,1);
    dvy = dv(:,:,2);
    dvz = dv(:,:,3);
    
    dnx = ((vy.*vy+vz.*vz).*dvx -vx.*vy.*dvy -vx.*vz.*dvz)./lencb;
    dny = (-vx.*vy.*dvx +(vx.*vx+vz.*vz).*dvy -vy.*vz.*dvz)./lencb;
    dnz = (-vx.*vz.*dvx -vy.*vz.*dvy +(vx.*vx+vy.*vy).*dvz)./lencb;

    % h1=1, h2=ny, h3=nz, h4=nx
    % h5=nx*ny, h6=ny*nz
    % h7=-nx*nx-ny*ny+2nz*nz
    % h8=nz*nx, h9=nx*nx-ny*ny
    dB_dnx = k.*(l(4)+l(5)*ny-l(7)*2*nx+l(8)*nz+l(9)*2*nx);
    dB_dny = k.*(l(2)+l(5)*nx+l(6)*nz-l(7)*2*ny-l(9)*2*ny);
    dB_dnz = k.*(l(3)+l(6)*ny+l(7)*4*nz+l(8)*nx);
    dB = dB_dnx.*dnx + dB_dny.*dny + dB_dnz.*dnz;
end


function [DEPTH_out] = DepthOffset(DEPTH, i, j)
    DEPTH_out.D = MatOffset(DEPTH.D, i, j);
    DEPTH_out.P = MatOffset(DEPTH.P, i, j);
    DEPTH_out.V = MatOffset(DEPTH.V, i, j);
    DEPTH_out.LEN = MatOffset(DEPTH.LEN, i, j);
    DEPTH_out.N = MatOffset(DEPTH.N, i, j);
end

function D = MatOffset(S, i, j)
    [h, w, c] = size(S);
    D = nan(h,w,c);
    dx = clamp(1-i, 1, w):clamp(w-i, 1, w);
    sx = clamp(1+i, 1, w):clamp(w+i, 1, w);
    dy = clamp(1-j, 1, h):clamp(h-j, 1, h);
    sy = clamp(1+j, 1, h):clamp(h+j, 1, h);
    D(dy, dx, :) = S(sy, sx, :); 
end

function y = clamp(x,minim, maxim)
    y = min(max(x, minim), maxim);
end