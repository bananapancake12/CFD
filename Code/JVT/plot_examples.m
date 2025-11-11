%   plot_examples                   
%                               
%   A script to create figures for 4A2 lecture course
%
%   Written by James Taylor                
%   September 2023

% Clear the workspace
clear; close all; clc;

% Figure directory
directory = 'C:\Users\jvt24\Documents\Teaching\4A2\Figures\';

% Load compressible flow tables
Comp_Tables


%% Draw different examples of mesh

% Create elliptical shape
th = linspace(0,90,101)'; e = 2; xy_surf = [-e * cosd(th) sind(th)]; n_surf = - m_norm(xy_surf);

% Open figure window
h.mesh = figure(); cols = lines(3); a = zeros(3,1);

% Plot the body multiple times
for n = 1:3
    a(n) = subplot(1,3,n,'visible','off'); hold on; axis equal; 
    axis([-2.5 -0.2 0.15 1.5]);
    plot(xy_surf(:,1),xy_surf(:,2),'k-')
end

% Plot the body conformal mesh
ni = 17; nj = 13; d = 1.5; si = linspace(0,1,ni)';
si_surf = m_dist(xy_surf,1); xy_1 = interp1(si_surf,xy_surf,si); 
n_1 = interp1(si_surf,n_surf,si); xy_2 = xy_1 + d * n_1;
sj = linspace(0,d,nj); xy = (1-sj) .* reshape(xy_1,[ni 1 2]) + sj .* reshape(xy_2,[ni 1 2]);
plot(a(1),xy(:,:,1),xy(:,:,2),'-','color',cols(1,:))
plot(a(1),xy(:,:,1)',xy(:,:,2)','-','color',cols(1,:))

% Plot an immersed boundary mesh
[x,y] = meshgrid(linspace(-2.6,0,17),linspace(0,1.6,13));
plot(a(3),x,y,'-','color',cols(3,:))
plot(a(3),x',y','-','color',cols(3,:))

% Plot an unstructured mesh
xy_list = [reshape(xy(:,1:3,1),[],1) reshape(xy(:,1:3,2),[],1)];
xy_rand = [xy_list ; -3 * rand(40,1) 2 * rand(40,1)];
q = inpolygon(xy_rand(:,1),xy_rand(:,2),[0 ; xy(:,3,1) ; 0],[0 ; xy(:,3,2) ; 0]);
xy_list = [xy_list ; xy_rand(q == 0,:)];
dt = delaunayTriangulation(xy_list);
axes(a(2)); triplot(dt,'color',cols(2,:))
%     -2.5 * ones(5,1) linspace(0,1.5,5)' ; linspace(-2.5,0,5)' -2.5 * ones(5,1)];

% Save the figure
% export_fig(h.mesh,[directory 'mesh_types.pdf'])


%% Fit surface to example data in IA question

% Input the data
[x,y] = ndgrid([1.9 2 2.1],[4.8 5 5.2]);
f = [112.37 116.09 119.97 ; 126.08 130 134.08 ; 140.99 145.11 149.39];


%% Plot 2nd and 4th order smoothing

% Plot part of a cosine wave and smoothing examples 
x = linspace(-0.2,0.8,5); y = sin(4*x); 
h.smooth = figure(); hold on; grid on; box on;
plot(x,y,'.','color',cols(1,:))

% Average values from second order smoothing
y2 = 0.5 * (y(1:end-2) + y(3:end));
plot(x(2:end-1),y2,'.','color',cols(2,:))

% Average values from fourth order smoothing
p = polyfit(x([1 2 4 5]),y([1 2 4 5]),3);
plot(x(3),polyval(p,x(3)),'.','color',cols(3,:))

% High resolution curves
x_interp = linspace(min(x),max(x),101);
y_act = sin(4*x_interp); plot(x_interp,y_act,'-','color',cols(1,:));
y_fit = polyval(p,x_interp); plot(x_interp,y_fit,'-','color',cols(3,:));

% Save the figure
% export_fig(h.smooth,[directory 'smooth_types.pdf'])

