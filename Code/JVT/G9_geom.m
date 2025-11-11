%   G9_geom
%
%   Assemble, plot and save the G9 turbine cascade used in the 4A3 lab

% Clear the workspace
close all; clear; clc;

% Read the raw geometry from file
load('cascade.mat')

% Plot the raw geometry
figure(); hold on; grid on; box on; axis equal;
plot(b.cfd.x,b.cfd.y,'.-')

% Cut out the trailing edge
i_te = [8 length(b.cfd.x)-6]; 
xy = [b.cfd.x(i_te(1):i_te(2)) b.cfd.y(i_te(1):i_te(2))];
plot(b.cfd.x(i_te),b.cfd.y(i_te),'o')

% Calculate cumilative distance along the curve
s = [0 ; cumsum(sum(diff(xy,1,1).^2,2).^0.5)]; s = s / s(end);

% Modify y-coordinates to create sharp trailing edge
d = cos(s * pi);
dxy_te = 0.5 * (xy(1,:) + xy(end,:)) - xy(1,:);
xy = xy + dxy_te .* d;

% Split into two sides at leading edge
[~,i_le] = min(xy(:,1));
xy_a = xy(i_le:end,:) / 10; 
xy_b = xy(i_le:-1:1,:) / 10;
plot(xy_a(:,1),xy_a(:,2),'.-')
plot(xy_b(:,1),xy_b(:,2),'.-')

% Save the geometry and pressures to file
filename = 'g9.raw';
fid = fopen(filename,'w'); fprintf(fid,'x_a y_a x_b y_b'); fclose(fid);
dlmwrite(filename,[xy_a xy_b],'-append','delimiter',' ','roffset',1,'precision','%10.8f')

