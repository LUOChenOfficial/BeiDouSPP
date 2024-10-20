clc;clear;close all;

%% read solution file
%fid=fopen("jfng0320.22o.pos");
fid=fopen("myPos.txt");
for i=1:8
header=fgetl(fid);
end
%sol_data=textscan(fid,"%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f");
sol_data=textscan(fid,"%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f");
fclose(fid);
sol_xyz=[sol_data{7},sol_data{8},sol_data{9}];
sol_real_xyz=[-2.27982915063712e+06 5.00470643812923e+06 3.21977736371519e+06];
error=zeros(size(sol_xyz,1),1);
det_enu=zeros(size(sol_xyz,1),3);
det_xyz=sol_xyz-sol_real_xyz;

%% calculate error for each epoch
for epochInd=1:size(sol_xyz,1)
    sol0_xyz=sol_xyz(epochInd,:);
    det0_xyz=sol0_xyz-sol_real_xyz;
    blh=xyz2blh(sol0_xyz);
    b=blh(1);
    l=blh(2);
    rotMat=[-sin(l),cos(l),0;-sin(b)*cos(l),-sin(b)*sin(l),cos(b);cos(b)*cos(l),cos(b)*sin(l),sin(b)];
    det_enu(epochInd,:)=(rotMat*det0_xyz.').';
    %error(epochInd,1)=sqrt(det_enu(epochInd,1)^2+det_enu(epochInd,2)^2+det_enu(epochInd,3)^2);
end

%% calculate rms
rms=sqrt(sum(det_enu.^2,1)/(size(sol_xyz,1)-1));

%% visualize error for whole epoches
epoch=(1:1:size(det_enu,1));
figure;
subplot(3,1,1);
plot(det_enu(:,1),"r");
ylabel("error(m)");
xlabel("epoch index");
title("east error for whole epoches");

subplot(3,1,2);
plot(det_enu(:,2),"g");
ylabel("error(m)");
xlabel("epoch index");
title("north error for whole epoches");

subplot(3,1,3);
plot(det_enu(:,3),"b");
ylabel("error(m)");
xlabel("epoch index");
title("up error for whole epoches");
