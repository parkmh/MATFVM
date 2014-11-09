function [cx cy cz hx hy hz] = ...
    generate_regular3D(xstart,xend,ystart,yend,zstart,zend,Nx,Ny,Nz)

% fprintf('Generating the center point of cells on [%d %d]x[%d %d]x[%d %d]\n',...
%     xstart,xend,ystart,yend,zstart,zend);

Nxy = Nx*Nz;
Nxyz = Nxy*Nz;

x = linspace(xstart,xend,Nx+1);
y = linspace(ystart,yend,Ny+1);
z = linspace(zstart,zend,Nz+1);

hx = x(2)-x(1);
hy = y(2)-y(1);
hz = z(2)-z(1);

cx = zeros(Nxyz,1);
cy = zeros(Nxyz,1);
cz = zeros(Nxyz,1);

for i = 1 : Nxyz
   ix = mod(mod(i-1,Nxy),Nx)+1;
   iy = floor(mod(i-1,Nxy)/Nx)+1;
   iz = floor((i-1)/Nxy)+1;
%    fprintf('(%d %d %d)\n',ix,iy,iz);
   cx(i) = (x(ix) + x(ix+1))/2;
   cy(i) = (y(iy) + y(iy+1))/2;
   cz(i) = (z(iz) + z(iz+1))/2;
end