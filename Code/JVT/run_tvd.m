%   run_tvd
%                               
%   A script to investigate TVD methods for 4A2 course
%
%   Written by James Taylor                
%   September 2023

% Clear the workspace
clear; close all; clc;

% Setup the case geometry
ni = 200; nt = 401;
L = 1; cfl = 0.5; x_node = linspace(0,L,ni+1)'; x_cell = 0.5 * (x_node(1:end-1) + x_node(2:end));
a = 1; dx = L / ni; dt = dx * cfl / a;

% Indices to reference fluxes to cell values
ip1 = [2:ni 1:2]; i = [1:ni 1]; im1 = [ni 1:ni]; im2 = [ni-1:ni 1:ni-1];

% Create figure windows
cols = [lines(7) ; 0.5 0.5 0.5];
figure(); ax.phi = axes; hold on; grid on; box on; xlabel('x'); ylabel('\phi'); axis([0 L -0.1 1.1]);
figure(); ax.lmter = axes; hold on; grid on; box on; xlabel('x'); ylabel('\sigma'); axis([0 L -0.1 2.1]);

% Initial condition with step function
phi = zeros(ni,1); phi(70:130) = 1; phi_exact = phi;

% Initial condition with sine pulse
% phi = zeros(ni,1); phi(70:130) = 0.5 * (1-cos(linspace(0,2*pi,61))); phi_exact = phi;

% Plot the exact solution one period later
plot(ax.phi,x_cell,phi_exact,'k.-')

% Loop over all methods 
varnames = {'Upwind' 'Central' 'Lax-Wendroff' 'Beam-Warming' 'Fromm' 'Minmod' 'Superbee' 'van-Leer'};
for v = 1:length(varnames)
    
    % Time march
    phi = phi_exact;
    for n = 1:nt
        
        % Ratio of successive gradients
        r = (phi(im1) - phi(im2)) ./ (phi(i) - phi(im1)); 
        r(isnan(r)) = 0;
        
        % Define fluxes for each method
        if strcmp(varnames{v},'Upwind') == 1
            f = phi(im1);
        elseif strcmp(varnames{v},'Central') == 1
            f = 0.5 * (phi(im1) + phi(i));
        elseif strcmp(varnames{v},'Lax-Wendroff') == 1
            f = phi(im1) + 0.5 * (1 - a * dt / dx) * (phi(i) - phi(im1));
        elseif strcmp(varnames{v},'Beam-Warming') == 1
            f = phi(im1) + 0.5 * (1 - a * dt / dx) * (phi(im1) - phi(im2));
        elseif strcmp(varnames{v},'Fromm') == 1
            f = phi(im1) + 0.25 * (1 - a * dt / dx) * (phi(i) - phi(im2));
        elseif strcmp(varnames{v},'Minmod') == 1
            sig = max(0,min(1,r));
            f = phi(im1) + 0.5 * (1 - a * dt / dx) * sig .* (phi(i) - phi(im1));
        elseif strcmp(varnames{v},'Superbee') == 1
            sig = max(0,max(min(2*r,1),min(r,2)));
            f = phi(im1) + 0.5 * (1 - a * dt / dx) * sig .* (phi(i) - phi(im1));
        elseif strcmp(varnames{v},'van-Leer') == 1
            r(isnan(r)) = 0; r(isinf(r)) = 0;
            sig = (r + abs(r)) ./ (1 + abs(r));
            f = phi(im1) + 0.5 * (1 - a * dt / dx) * sig .* (phi(i) - phi(im1));
        end

        % Sum the fluxes to calculate cell changes
        phi = phi + a * dt * (f(1:end-1) - f(2:end)) / dx;
        
    end
    
    % Plot the result after one period
    plot(ax.phi,x_cell,phi,'.-','color',cols(v,:)) 
    if v > 5; plot(ax.lmter,x_node,sig,'.-','color',cols(v,:)); end;
    
end

% Add the legend
legend(ax.phi,['Exact' varnames])




