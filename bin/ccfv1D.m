function [A b] = ccfv1D(cx,hx,Nx,k,f,gD,gN,bnd)
% fprintf('\n[ 2D Cell Centered Finite Volume Method ]\n')
b = zeros(Nx,1);

IA = zeros(3*Nx,1);
JA = zeros(3*Nx,1);
AA = zeros(3*Nx,1);
l = 1;

% Left Boundary
center = 0;
if bnd(1) == 1
    center = center + k(1)*2/hx^2;
    b(1) = b(1) + feval(gD{1},cx(1) - hx/2)*k(1)*2/hx^2;
elseif bnd(1) == 2
    b(1) = b(1) - feval(gN{1},cx(1) - hx/2)/hx;
else
    errid = 'MATLAB:funfun:ccfv1D:notABnd';
            errmsg = sprintf('Invalid value for boundary : value must be ''1'' or ''2''.');
            error(errid,errmsg);
end

% Right
value = 2 / ((1/k(1)+1/k(1+1)) * hx^2);
center = center + value;
IA(l) = 1;
JA(l) = 2;
AA(l) = -value;
l = l + 1;
        
% Center
IA(l) = 1;
JA(l) = 1;
AA(l) = center;
l = l + 1;

for i = 2 : Nx - 1
    b(i) = feval(f,cx(i));
    
    center = 0;
    
    % Left
    value = 2 / ((1/k(i)+1/k(i-1)) * hx^2);
    center = center + value;
    IA(l) = i;
    JA(l) = i-1;
    AA(l) = -value;
    l = l + 1;
    
    % Right
    value = 2 / ((1/k(i)+1/k(i+1)) * hx^2);
    center = center + value;
    IA(l) = i;
    JA(l) = i+1;
    AA(l) = -value;
    l = l + 1;
    
    %Center
    IA(l) = i;
    JA(l) = i;
    AA(l) = center;
    l = l + 1;
end

center = 0;
if bnd(2) == 1
    center = center + k(Nx)*2/hx^2;
    b(Nx) = b(Nx) + feval(gD{2},cx(Nx) + hx/2)*k(Nx)*2/hx^2;
elseif bnd(2) == 2
    b(Nx) = b(Nx) - feval(gN{2},cx(Nx) + hx/2)/hx;
else
    errid = 'MATLAB:funfun:ccfv1D:notABnd';
            errmsg = sprintf('Invalid value for boundary : value must be ''1'' or ''2''.');
            error(errid,errmsg);
end

value = 2 / ((1/k(Nx)+1/k(Nx-1)) * hx^2);
center = center + value;
IA(l) = Nx;
JA(l) = Nx-1;
AA(l) = -value;
l = l + 1;
        
% Center
IA(l) = Nx;
JA(l) = Nx;
AA(l) = center;
l = l + 1;
A = sparse(IA(1:l-1),JA(1:l-1),AA(1:l-1), Nx,Nx);