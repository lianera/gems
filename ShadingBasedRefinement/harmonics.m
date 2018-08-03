function [ A ] = harmonics( VN )
%HARMONICS 此处显示有关此函数的摘要
%   此处显示详细说明
    [rows, ~] = size(VN);
    
    nx = VN(:, 1);
    ny = VN(:, 2);
    nz = VN(:, 3);
    h0 = ones(rows, 1);
    h1 = ny;
    h2 = nz;
    h3 = nx;
    h4 = nx.*ny;
    h5 = ny.*nz;
    h6 = -nx.*nx - ny.*ny + 2*nz.*nz;
    h7 = nz.*nx;
    h8 = nx.*nx - ny.*ny;
    
    A = [h0, h1, h2, h3, h4, h5, h6, h7, h8];
end

