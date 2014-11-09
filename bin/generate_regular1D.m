function [cx hx] = generate_regular1D(xstart,xend,Nx)
% GENERATE_REGULAR1D generates x values of center of cells
%
%   [CX CY HX HY] = GENERATE_REGULAR2D(XSTART,XEND,YSTART,YEND,NX,NY)
%
%   See also 

% Written by Minho Park
% Date : 5-9-2011

% fprintf('Generating the center point of cells on [%d %d]x[%d %d]\n',...
%     xstart,xend,ystart,yend);

x = linspace(xstart,xend,Nx+1);

hx = x(2)-x(1);

cx = zeros(Nx,1);


for i = 1:Nx
    cx(i)   = (x(i)+x(i+1))/2;
end
