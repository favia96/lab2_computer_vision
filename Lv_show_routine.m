function [pixels,Lx,Ly] = Lv_show_routine(img,shape,dxmask,dymask,threshes,smooth_flag,vars)

    if (nargin < 2)
        shape = 'same';
    end
    
    smoothimg = zeros(size(img,1),size(img,2),6);
    Lx = zeros(size(img,1),size(img,2),6);
    Ly = zeros(size(img,1),size(img,2),6);
    pixels = zeros(size(img,1),size(img,2),6);

    %figure()
    %subplot(2,3,1); showgrey(img); title('original');
    
    figure()
    if(smooth_flag == true)
        for i = 1 : length(vars)
            if i == 1 
                smoothimg(:,:,i) = img;
            else
                smoothimg(:,:,i) = gaussfft(img,vars(i));
            end
            
            Lx(:,:,i) = conv2(smoothimg(:,:,i), dxmask, shape);
            Ly(:,:,i) = conv2(smoothimg(:,:,i), dymask, shape);
            pixels(:,:,i) = sqrt(Lx(:,:,i).^2 + Ly(:,:,i).^2);
            
            if i == 1
                subplot(2,3,i); showgrey(pixels(:,:,i)); title('grad magnitude');
            else
                subplot(2,3,i); showgrey((pixels(:,:,i) - threshes(5)) > 0); title(sprintf('sigma^2 = %.4f', vars(i)));
            end
        end
    end
    sgtitle(sprintf('Different Gaussian sigma^2, thresh = %.1f',threshes(5)));
    
%     Lx = conv2(img_1, dxmask, shape);
%     subplot(2,3,2); showgrey(Lx); title('1st deriv x');
    
%     Ly = conv2(img_1, dymask, shape);
%     subplot(2,3,3); showgrey(Ly); title('1st deriv y');
    
%     pixels_1 = sqrt(Lx.^2 + Ly.^2);
%     figure()
%     subplot(2,3,1); showgrey(pixels_1); title('grad magnitude');
    
%     subplot(2,3,5); g = histogram(pixels); title('hist grad magnitude');
    
%     subplot(2,3,2); showgrey((pixels_1 - thresh_1) > 0); title('thresh = 10')
%     subplot(2,3,3); showgrey((pixels_1 - thresh_2) > 0); title('thresh = 25')
%     subplot(2,3,4); showgrey((pixels_1 - thresh_3) > 0); title('thresh = 50')
%     subplot(2,3,5); showgrey((pixels_1 - thresh_4) > 0); title('thresh = 75')
%     subplot(2,3,6); showgrey((pixels_1 - thresh_5) > 0); title('thresh = 100')

end
