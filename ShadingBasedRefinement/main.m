CamInfo.focal_len = 0.035;
CamInfo.sensor_w = 0.036;
CamInfo.sensor_h = 0.024;

I = im2double(imread('data/rgb.png'));
D = double(hdrread('data/depth.hdr'));
D = D(:,:,1);
D(D>100) = nan;

[h,w,~] = size(I);
CamInfo.w = w;
CamInfo.h = h;
CamInfo.cx = w/2;
CamInfo.cy = h/2;

CamInfo.fx = w*CamInfo.focal_len / CamInfo.sensor_w;
CamInfo.fy = h*CamInfo.focal_len / CamInfo.sensor_h;

% node: w_here = sqrt(w_paper)
weights.wg = 1;
weights.ws = sqrt(400);
weights.wp = sqrt(10);

[DEPTH, mask] = D2N(D, CamInfo);
imshow(DEPTH.Nimg);
PC = pointCloud(reshape(DEPTH.P, [h*w, 3]));
pcwrite(PC,'origin','PLYFormat','binary');

%--- refine ---
D_refined = refine(I, D, CamInfo, weights);

% unit: meter
D_refined(~isfinite(D_refined)) = 0;
hdrwrite(cat(3,D_refined,D_refined,D_refined), 'refined_depth.hdr');
DEPTH_refined = D2N(D_refined, CamInfo);
figure;
imshow(DEPTH_refined.Nimg);
PC_refined = pointCloud(reshape(DEPTH_refined.P, [h*w, 3]));
pcwrite(PC_refined,'refined','PLYFormat','binary');