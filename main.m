%% load image 
clc
clear all
close all
[ambimage,imarray,lightdirs] = LoadFaceImages();

% preprocessing
shape = size(imarray);
for i = 1:(shape(3))
    imarray(:,:,i) = (imarray(:,:,i)- ambimage) / 255; 
end
imarray(imarray < 0) = 0;

imarray_p = reshape(imarray,shape(1)*shape(2),shape(3));
t_lightdirs = transpose(lightdirs);

%% cf) answer
% [albedo_image, surface_normals] = photometric_stereo(imarray, t_lightdirs);
% imshow(albedo_image,[]);
% h = show_surfNorm(surface_normals,10);
% height_map = get_surface(surface_normals,[shape(1),shape(2)],'average');
% % height_map(height_map>1000)=0
% display_output(albedo_image, height_map);

%% compute albedo, normals
normaldirs = t_lightdirs\transpose(imarray_p);
% normaldirs = transpose(normaldirs);
albedo = (vecnorm(normaldirs));
albedo(albedo==0) =0.1;
normal_Vec = transpose((albedo.^-1) .* normaldirs);

albedo = albedo';
albedo_2d = reshape(albedo,shape(1),shape(2),1);
normal_2d = reshape(normal_Vec,shape(1),shape(2),3);
imshow(albedo_2d);


[X, Y] = meshgrid(1:10:shape(2), shape(1):-10:1);
U = normal_2d(1:10:shape(1), shape(2):-10:1, 1);
V = normal_2d(1:10:shape(1), shape(2):-10:1, 2);
W = normal_2d(1:10:shape(1), shape(2):-10:1, 3);
figure()
quiver3(X, Y, zeros(size(X)), U, V, W);
%%
normal_2d(normal_2d(:,:,3)==0)=1;
d_fx = normal_2d(:,:,1)./normal_2d(:,:,3);
d_fy = normal_2d(:,:,2)./normal_2d(:,:,3);
d_fx(abs(d_fx)>10) = 0.0001;
d_fy(abs(d_fy)>10) = 0.0001;
depth = zeros(shape(1),shape(2));
depth_1 = zeros(shape(1),shape(2));
depth_2 = zeros(shape(1),shape(2));
depth_3 = zeros(shape(1),shape(2));
depth_4 = zeros(shape(1),shape(2));
for i = 1:shape(1)
    for j = 1:shape(2)
        r_sum = cumsum(d_fy(1:i,1),1);
        c_sum = cumsum(d_fx(i,1:j),2); 
        depth_1(i,j) = (hypot(c_sum(end),r_sum(end)));
        
        r_sum = cumsum(d_fy(end:-1:i,1),1);
        c_sum = cumsum(d_fx(i,1:j),2); 
        depth_2(i,j) = (hypot(c_sum(end),r_sum(end)));
        
        r_sum = cumsum(d_fy(end:-1:i,1),1);
        c_sum = cumsum(d_fx(i,end:-1:j),2); 
        depth_3(i,j) = (hypot(c_sum(end),r_sum(end)));
        
        r_sum = cumsum(d_fy(1:i,1),1);
        c_sum = cumsum(d_fx(i,end:-1:j),2); 
        depth_4(i,j) = (hypot(c_sum(end),r_sum(end)));         
        depth(i,j) = 1/4* ( depth_1(i,j) + depth_2(i,j) +depth_3(i,j)+ depth_4(i,j));
%         depth(i,j) = hypot(depth_1(i,j),depth_2(i,j));
%         depth(i,j) = hypot(depth(i,j),depth_3(i,j));
%         depth(i,j) = hypot(depth(i,j),depth_4(i,j));
    end
end
% 6
% depth(depth>300)=0.01;
display_output(albedo_2d, depth);
