function [ RGBout, Dout ] = kinect_calib( RGB, D)
    % reference: http://www.mdpi.com/1424-8220/15/11/27569/pdf
        
    % unit millimeter
    
    rgb_flen = 3.291;
    rgb_sensor_w = 6.00;
    rgb_sensor_h = 3.38;
    rgb_w = 1920;
    rgb_h = 1080;
    
    d_flen = 3.657;
    d_sensor_w = 5.12;
    d_sensor_h = 4.24;
    d_w = 512;
    d_h = 424;
    
    
    baseline = 52.5;

    % cast parameters
    rgb2d = [1,0,0,-baseline;...
             0,1,0,0;...
             0,0,1,0;...
             0,0,0,1];
    d2rgb = inv(rgb2d);
        
    rgbfx = rgb_flen*rgb_w/rgb_sensor_w;
    rgbfy = rgb_flen*rgb_h/rgb_sensor_h;
    rgbcx = rgb_w/2;
    rgbcy = rgb_h/2;
    
    dfx = d_flen*d_w/d_sensor_w;
    dfy = d_flen*d_h/d_sensor_h;
    dcx = d_w/2;
    dcy = d_h/2;
    
    % project depth to world
    [i,j] = meshgrid(1:d_w,1:d_h);
    dz = double(D);
    dx = dz.*(i-dcx)/dfx;
    dy = dz.*(j-dcy)/dfy;
    dp = reshape(cat(3, dx, dy, dz, ones(d_h,d_w)), d_h*d_w, 4)';
    
    % camera transform
    rgbp = (d2rgb * dp)';
        
    % world to rgb camera sensor
    rgbi = rgbp(:,1)*rgbfx./rgbp(:,3)+rgbcx;
    rgbj = rgbp(:,2)*rgbfy./rgbp(:,3)+rgbcy;
    rgbz = rgbp(:,3);
        
    % rasterlize
    rgbi = round(rgbi);
    rgbj = round(rgbj);
    mask = isfinite(rgbi) & rgbi>=1 & rgbi<=rgb_w ...
         & isfinite(rgbj) & rgbj>=1 & rgbj<=rgb_h;
     
    rgbi = rgbi(mask);
    rgbj = rgbj(mask);
    rgbz = rgbz(mask);

    % z-buffer sort
    [~,order] = sort(rgbz, 'descend');
    rgbi = rgbi(order);
    rgbj = rgbj(order);

    % get color for depth
    
    idx = (rgbi-1)*rgb_h+rgbj;
    R(order) = RGB(idx);
    G(order) = RGB(idx+rgb_w*rgb_h);
    B(order) = RGB(idx+rgb_w*rgb_h*2);
    
    Rfull = nan(d_h, d_w);
    Rfull(mask) = R;
    Gfull = nan(d_h, d_w);
    Gfull(mask) = G;
    Bfull = nan(d_h, d_w);
    Bfull(mask) = B;
    
    RGBout = cat(3, Rfull, Gfull, Bfull);
    
    % unit meter
    Dout = D / 1000;
end

