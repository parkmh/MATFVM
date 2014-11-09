function [cx cy hx hy] = generate_regular2D(xstart,xend,ystart,yend,Nx,Ny)
% GENERATE_REGULAR2D generates x and y values of center of cells
%
%   [CX CY HX HY] = GENERATE_REGULAR2D(XSTART,XEND,YSTART,YEND,NX,NY)
%
%   See also 

% Written by Minho Park
% Date : 9-5-2011

% fprintf('Generating the center point of cells on [%d %d]x[%d %d]\n',...
%     xstart,xend,ystart,yend);
Nxy = Nx*Ny;

x = linspace(xstart,xend,Nx+1);
y = linspace(ystart,yend,Ny+1);

hx = x(2)-x(1);
hy = y(2)-y(1);

cx = zeros(Nxy,1);
cy = zeros(Nxy,1);


for i = 1:Nxy
    ix = mod(i-1,Nx)+1;
    iy = floor((i-1)/(Nx))+1;
    cx(i)   = (x(ix)+x(ix+1))/2;
    cy(i)   = (y(iy)+ y(iy+1))/2;
end