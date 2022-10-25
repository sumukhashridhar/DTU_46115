clear all; clc;

%% column descriptions

%      1       : y/delta Grid point in wall-normal direction, Normalized by channel half width 
%      2       : y^+ Grid point in wall-normal direction 
%      3       : U Mean profile of streamwise velocity 
%      4       : dU/dy Derivative of U in wall-normal direction 
%      5       : W Mean profile of spanwise velocity 
%      6       : P Mean profile of pressure 

%% load data

file_550 = load('data_Re550.csv');
file_180 = load('data_Re180.csv');

%% velocity

u_550 = file_550(:,3);
u_180 = file_180(:,3);

%% y+

yplus_180 = file_180(:,2);

yplus_550 = file_550(:,2);

figure
plot(yplus_180, u_180);
title('Re 180')
xlabel('y plus') 
ylabel('U') 

figure
plot(yplus_550, u_550);
title('Re 550')
xlabel('y plus') 
ylabel('U') 
