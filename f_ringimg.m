function [ringimg,ringimg_filled]=f_ringimg()
addpath(genpath('W:\projects\alireza\core_toolbox'))

%% params
ringwidth=0.35;  % if possible run 0.3 too
ringsize=40; 
rings=[10 20 30];


[x,y]=meshgrid(-ringsize:ringsize,-ringsize:ringsize);
[~,r]=cart2pol(x,y);
ringimg=x*0;
[xc,yc]=find(r==0); %center
ring_im_size=length(r);

for ring=rings
    ringimg(ring-ringwidth < r & ring+ringwidth > r)=1;
end
ringimg_filled=x*0;
ringimg_filled(r<=rings(3))=1;
% imshow(ringimg_filled)
%% tile rings
ringimgs=[];
ringimgs_filled=[];
tilenum_ver=4;
tilenum_hor=8;
pts = hexagonalGrid([0 0 (tilenum_hor)*ring_im_size (tilenum_ver)*ring_im_size], [xc yc], ring_im_size);
pts(:,2)=ceil(pts(:,2));
pts(:,1)=ceil(pts(:,1));

for i=1:length(pts)
    ringimgs(pts(i,2)-(xc-1):pts(i,2)+(xc-1),pts(i,1)+26:pts(i,1)+(ring_im_size+25))=ringimg*0;
    ringimgs_filled(pts(i,2)-(xc-1):pts(i,2)+(xc-1),pts(i,1)+26:pts(i,1)+(ring_im_size+25))=ringimg_filled*0;
end

for i=1:length(pts)
    ringimgs(pts(i,2)-(xc-1):pts(i,2)+(xc-1),pts(i,1):pts(i,1)+(ring_im_size-1))=ringimgs(pts(i,2)-(xc-1):pts(i,2)+(xc-1),pts(i,1):pts(i,1)+(ring_im_size-1))+ringimg;
    ringimgs_filled(pts(i,2)-(xc-1):pts(i,2)+(xc-1),pts(i,1):pts(i,1)+(ring_im_size-1))=ringimgs_filled(pts(i,2)-(xc-1):pts(i,2)+(xc-1),pts(i,1):pts(i,1)+(ring_im_size-1))+ringimg_filled;
end

% imagesc(ringimgs_filled);
im_size   = size(ringimgs);
stimrect  = [1200 2400];
Scrn_size = [1200 1920];
portion   = floor(im_size.*(Scrn_size./stimrect));

ringimg = ringimgs(1:portion(1),1:portion(2));
ringimg_filled = ringimgs_filled(1:portion(1),1:portion(2));

% figure;imagesc(ringimg);
% y = linspace(1,length(ringimg(:,1)),length(ringimgs(:,1)));
% imagesc([1 26],[1 16],ringimg)
% x = linspace(1,length(ringimg(1,:)),length(ringimgs(1,:)));