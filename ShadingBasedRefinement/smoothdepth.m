function [ D ] = smoothdepth( S, ksize, depth_threshold)
    T = sm1d(S', ksize, depth_threshold);
    D = sm1d(T', ksize, depth_threshold);
end

function D = sm1d(S, ksize, depth_threshold)
    [rows, cols] = size(S);
    D = nan(rows, cols);
    sigma = ksize / 4;
    kernel = repmat(fspecial('gaussian', [ksize,1], sigma)', rows, 1);
    halfsize = (ksize-1)/2;
    for i = 1:cols
        cov = zeros(rows, 1);
        cur = S(:,i);
        for k = -halfsize:halfsize        
            j = i+k;            
            if j < 1 || j > cols
                c = cur;
            else
                c = S(:,j);
                nfidx = find(~isfinite(c) | abs(c-cur)>depth_threshold);
                c(nfidx) = cur(nfidx);
            end
            cov = cov + c.*kernel(:,k+halfsize+1);
        end
        D(:, i) = cov;
    end
    
end