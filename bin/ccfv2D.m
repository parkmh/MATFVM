function [A b] = ccfv2D(cx,cy,hx,hy,Nx,Ny,k,f,gD,gN,bnd)
% fprintf('\n[ 2D Cell Centered Finite Volume Method ]\n')
Nxy = Nx*Ny;
b = zeros(Nxy,1);

IA = zeros(5*Nxy,1);
JA = zeros(5*Nxy,1);
AA = zeros(5*Nxy,1);
l = 1;

for i = 1 : Nxy
    ix = mod(i-1,Nx)+1;
    iy = floor((i-1)/Nx)+1;

    b(i) = feval(f,cx(i),cy(i));
    
    center = 0;
    if ix == 1        
        if bnd(3) == 1  % Dirichlet
            center = center + k(ix,iy)*2/hx^2;
            b(i) = b(i) + feval(gD{3},cx(i) - hx/2,cy(i))*k(ix,iy)*2/hx^2;
        elseif bnd(3) == 2                
            b(i) = b(i) - feval(gN{3},cx(i) - hx/2,cy(i))/hx;
        else
            errid = 'MATLAB:funfun:ccfv2D:notABnd';
            errmsg = sprintf('Invalid value for boundary : value must be ''1'' or ''2''.');
            error(errid,errmsg);
        end
    else        
        value = 2 / ((1/k(ix,iy)+1/k(ix-1,iy)) * hx^2);
        center = center + value;
        IA(l) = i;
        JA(l) = i-1;
        AA(l) = -value;
        l = l + 1;
    end
    
    if ix == Nx
        if bnd(2) == 1
            center = center + k(ix,iy)*2/hx^2;
            b(i) = b(i) + feval(gD{2},cx(i) + hx/2,cy(i))*k(ix,iy)*2/hx^2;
        elseif bnd(2) == 2
            b(i) = b(i) - feval(gN{2},cx(i) + hx/2,cy(i))/hx;
        else
            errid = 'MATLAB:funfun:ccfv2D:notABnd';
            errmsg = sprintf('Invalid value for boundary : value must be ''1'' or ''2''.');
            error(errid,errmsg);
        end
    else
        value = 2 / ((1/k(ix,iy)+1/k(ix+1,iy)) * hx^2);
        center = center + value;
        IA(l) = i;
        JA(l) = i+1;
        AA(l) = -value;
        l = l + 1;
    end
    
    % South
    if iy == 1
        if bnd(4) == 1           
            center = center + k(ix,iy)*2/hy^2;
            b(i) = b(i) + feval(gD{4},cx(i),cy(i)-hy/2)*k(ix,iy)*2/hy^2;
        elseif bnd(4) == 2
            b(i) = b(i) - feval(gN{4},cx(i),cy(i)-hy/2)/hy;
        else
            errid = 'MATLAB:funfun:ccfv2D:notABnd';
            errmsg = sprintf('Invalid value for boundary : value must be ''1'' or ''2''.');
            error(errid,errmsg);
        end
    else
        value = 2 / ((1/k(ix,iy)+1/k(ix,iy-1)) * hy^2);
        
        center = center + value;
        IA(l) = i;
        JA(l) = i-Nx;
        AA(l) = -value;
        l = l + 1;        
    end
    
    % North
    if iy == Nx
        if bnd(1) == 1
            center = center + k(ix,iy)*2/hy^2;
            b(i) = b(i) + feval(gD{1},cx(i),cy(i)+hy/2)*k(ix,iy)*2/hy^2;
        elseif bnd(1) == 2
            b(i) = b(i) - feval(gN{1},cx(i),cy(i)+hy/2)/hy;
        else
            errid = 'MATLAB:funfun:ccfv2D:notABnd';
            errmsg = sprintf('Invalid value for boundary : value must be ''1'' or ''2''.');
            error(errid,errmsg);
        end
    else
        value = 2 / ((1/k(ix,iy)+1/k(ix,iy+1)) * hy^2);
        center = center + value;
        IA(l) = i;
        JA(l) = i+Nx;
        AA(l) = -value;
        l = l + 1;
    end
    
    % Center
    IA(l) = i;
    JA(l) = i;
    AA(l) = center;
    l = l + 1;
end

A = sparse(IA(1:l-1),JA(1:l-1),AA(1:l-1), Nxy,Nxy);