%% LAB 2 - COMPUTER VISION, November 2019
%% by Martin De Pellegrini, Federico Favia

%% Initialization
clear ; close all; clc
scale = [0.0001, 1.0, 4.0, 16.0, 64.0]; % variances of gaussian smoothing

%% Part 1: Difference operators
% simple difference operator
deltax = [-1, 0, 1]; deltay = deltax';

% central difference operator
deltax_c = [-0.5, 0, 0.5]; deltay_c = deltax_c';

% diagonal Roberts operator
rob_pos_dig = [-1 0; 
                0 1];
rob_neg_dig = [0 -1; 
               1 0];

%% load image and compute 1st ord deriv
tools = few256;
figure('name','SDO, CDO, and grad. magn for tools')
subplot(2,3,1); showgrey(tools); title('image');

dxtools = conv2(tools, deltax, 'same'); % try with valid to see differences of dimensions
dytools = conv2(tools, deltay, 'same');

subplot(2,3,2); showgrey(dxtools); title('simple, x direction');
subplot(2,3,3); showgrey(dytools); title('simple, y direction');
dxtools = conv2(tools, deltax_c, 'same'); 
dytools = conv2(tools, deltay_c, 'same');
subplot(2,3,4); showgrey(dxtools); title('central, x direction');
subplot(2,3,5); showgrey(dytools); title('central, y direction');

% Part 2: Pointwise thresholding of gradient magnitudes
gradmagntools = sqrt(dxtools .^2 + dytools .^2);
subplot(2,3,6); showgrey(gradmagntools); title('grad magnitude');
sgtitle('SDO, CDO, and grad. magn for tools');

figure('name','Hist of grad. magn and thresholded for tools')
subplot(1,2,1); g = histogram(gradmagntools); title('hist grad magnitude'); 
thresh = 20; % thresh to obtain thin edges
subplot(1,2,2); showgrey((gradmagntools - thresh) > 0); title('thresh grad magnitude');
sgtitle('Hist of grad. magn and thresholded for tools');

%% procedure to simplify process
vars = [0,0.45,0.55,0.6,0.8,0.95]; % trial-and-error
threshes_q2 = [10, 25, 50, 75, 100];
threshes_q3 = [10,30,50,75,90];
house = godthem256;
[pixels,Lx,Ly] = Lv_show_routine(house,'same',deltax_c,deltay_c,threshes_q3,true,vars); 
% 5 threshes, flag true if gaussian smoothing yes, 5 variances
% if We don't pass 'same' problem of matrix dimension because of conv2

%% Part 4: Differential geometry descriptors
% checking functions
dx = [0 0 0 0 0; 0 0 0 0 0; 0 -0.5 0 0.5 0; 0 0 0 0 0; 0 0 0 0 0];
dy = dx';
dxx = [0 0 0 0 0; 0 0 0 0 0; 0 1 -2 1 0; 0 0 0 0 0; 0 0 0 0 0];
dyy = dxx';
% dxy = filter2(dx, dy, 'same');
dxy = conv2(dy, dx, 'same');
dxxx = conv2(dx, dxx, 'same');
dxxy = conv2(dxx, dy, 'same');
dxyy = conv2(dx, dyy, 'same');
dyyy = conv2(dy, dyy, 'same');

[x y] = meshgrid(-5 : 5);

figure('name', 'Test dxxx, dxx, dxxy')
subplot(1,3,1); showgrey(filter2(dxxx,x.^3,'valid'));
subplot(1,3,2); showgrey(filter2(dxx,x.^3,'valid'));
subplot(1,3,3); showgrey(filter2(dxxy,x.^2 .* y,'valid'));
sgtitle('Test dxxx, dxx, dxxy');

%% plotting
house = godthem256;

figure('name', 'Lvv for house')
for i = 1 : length(scale)
    subplot(2,3,i);
    contour(Lvvtilde(discgaussfft(house,scale(i)),'same'),[0 0]); % check when lvv == 0
    if i == 1
            title(sprintf('Lvv house, sigma^2 = %.4f', scale(1))); % 4 decimal digits needed just once
    else
            title(sprintf('Lvv house, sigma^2 = %.1f', scale(i))); % enough 1 decimal digit
    end
    axis('image')
    axis('ij')
end
sgtitle('2nd derivative = 0');

tools = few256;
figure('name', 'Lvvv for tools')
for i = 1 : length(scale)
        subplot(2,3,i);
        showgrey(Lvvvtilde(discgaussfft(tools,scale(i)),'same') < 0);
        if i == 1
            title(sprintf('Lvvv tools, sigma^2 = %.4f', scale(1))); % 4 decimal digits needed just once
        else
            title(sprintf('Lvvv tools, sigma^2 = %.1f', scale(i))); % enough 1 decimal digit
        end
end
sgtitle('3rd derivative < 0');

%% quest. 6: try to combine Lvv to detect edges + improve results using Lvvv
figure('name', 'lvv + lvvv for house')
for i = 1 : length(scale)
    subplot(2,3,i)
    lvv = Lvvtilde(discgaussfft(house, scale(i)), 'same');
    lvvv = Lvvvtilde(discgaussfft(house, scale(i)), 'same');
%      m = lvv;
    m = lvv.*real(log(1+lvv));
    m(lvvv > 0) = NaN; % check only when 3rd is less than zero
    contour(m, [0 0]);
    axis('image') % without these the image is rotated, strange
    axis('ij')
    
    if i == 1
        title(sprintf('Lvv + Lvvv house, sigma^2 = %.4f', scale(1)));
    else
        title(sprintf('Lvv + Lvvv for house, sigma^2 = %.1f', scale(i)));
    end
end
sgtitle('Lvv + Lvvv for house');

%% Part 5: Extraction of edge segments
thresh_house = 3.5;
figure('name', 'Extracted edges for house')
for i = 1 : length(scale)
        subplot(2,3,i);
        extractedge(house, scale(i), thresh_house, 'same');
        if i == 1
            title(sprintf('Edges on house, sigma^2 = %.4f, thresh = %.1f', scale(1),thresh_house));
        else
            title(sprintf('Edges on house, sigma^2 = %.1f, thresh = %.1f', scale(i),thresh_house));
        end
end
sgtitle(sprintf('Extract. edges on house, thresh = %.1f',thresh_house));

thresh_tools = 8;
figure('name', 'Extracted edges for tools')
for i = 1 : length(scale)
        subplot(2,3,i);
        extractedge(tools, scale(i), thresh_tools, 'same');
        if i == 1
            title(sprintf('Edges on tools, sigma^2 = %.4f, thresh = %.1f', scale(1),thresh_tools)); 
        else
            title(sprintf('Edges on tools, sigma^2 = %.1f, thresh = %.1f', scale(i),thresh_tools));
        end
end
sgtitle(sprintf('Extract. edges on tools, thresh = %.1f',thresh_tools));

%% Part 6: Hough transform
triangle = triangle128;
mag = Lv(triangle, true, scale(3),'same');
[linepar, acc] = houghline(zerocrosscurves(triangle-128), mag, size(triangle,1), size(triangle,2),10,3); % thresh=10, 5 lines

%figure(); showgrey(acc); % show accumulator space
visual_output(triangle,linepar); % show edge lines

%houghedgeline(tr,scale(4),8,64,64,5,1);

%%  Question 8: hough edges
% testimage1 = triangle128;
% smalltest1 = binsubsample(testimage1);
% 
% [lptest1, acctest1] = houghline(zerocrosscurves(testimage1), Lv(testimage1, true, scale(3),'same'), size(testimage1,1), size(testimage1,2),10,5);
% visual_output(testimage1,lptest1)
% pause;
% showgrey(acctest1)
% 
% [lptestsmall1, acctestsmall1] = houghline(zerocrosscurves(smalltest1), Lv(smalltest1, true, scale(3),'same'), size(smalltest1,1), size(smalltest1,2),10,5);
% visual_output(smalltest1,lptestsmall1)
% pause;
% showgrey(acctestsmall1)
% 
% testimage2 = houghtest256;
% % smalltest2 = binsubsample(binsubsample(testimage2));
% 
% [lptest2, acctest2] = houghline(zerocrosscurves(testimage2), Lv(testimage2, true, scale(3),'same'), size(testimage2,1), size(testimage2,2),10,10);
% visual_output(testimage2,lptest2)
% figure;
% showgrey(acctest2)
% 
% [pos, value] = locmax8(acctest2);
% [dummy, indexvector] = sort(value);
% nmaxima = size(value, 1);

% [lptestsmall2, acctestsmall2] = houghline(zerocrosscurves(smalltest2), Lv(smalltest2, true, scale(3),'same'), size(smalltest2,1), size(smalltest2,2),10,20);
% visual_output(smalltest2,lptestsmall2)
% pause; showgrey(acctestsmall2)

%% Question 8:
testimage1 = triangle128;
testimage2 = houghtest256;

[lptest1, acctest1] = houghedgeline(testimage1,scale(3),5,size(testimage1,1)*2, size(testimage1,2)*2,3,2); %3 lines, verbose 2 2

[pos, value] = locmax8(acctest1);
[dummy, indexvector] = sort(value);
nmaxima = size(value, 1);

houghedgeline(testimage2,scale(3),5,size(testimage2,1)*4, size(testimage2,2)*4,10,2);

%% Question 9: number of cells in accumulator
phonecalc = phonecalc256;

for i = 1 : 3 % 3 changes in number of cells
    houghedgeline(phonecalc, scale(3), 10, size(phonecalc,1)*(i^2),size(phonecalc,2)*(i^2),10,2); % increment of number cells, ten lines, verbose 2
end

%% question 10: choice accumulator incrementation function (change inside houghline)
houghedgeline(few256, scale(3), 8, size(few256,1)*2,size(few256,2)*2,10,2); %ten lines, verbos 2


