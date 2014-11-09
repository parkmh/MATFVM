function [A b] = ccfv3D(cx,cy,cz,hx,hy,hz,Nx,Ny,Nz,k,f,gD,gN,bnd)
% bnd = [Left(x<0) Right(x>0) Front(y<0) Back(y>0) Down(z<0) Up(z>0)]
% fprintf('\n[ 3D Cell Centered Finite Volume Method ]\n')

Nxy = Nx*Ny;
Nxyz = Nxy*Nz;
b = zeros(Nxyz,1);

IA = zeros(7*Nxyz,1);
JA = zeros(7*Nxyz,1);
AA = zeros(7*Nxyz,1);

l = 1;

for i = 1 : Nxyz
    ix = mod(mod(i-1,Nxy),Nx)+1;
    iy = floor(mod(i-1,Nxy)/Nx)+1;
    iz = floor((i-1)/Nxy)+1;
%     fprintf('\n(%d,%d,%d)\n')
    b(i) = feval(f,cx(i),cy(i),cz(i));
    center = 0;
    
%     fprintf('[ Left ]\n')
    if ix == 1
        if bnd(1) == 1 
%             fprintf('Dirichlet Boundary : ')
            center = center + k(ix,iy,iz)*2/hx^2;
            b(i) = b(i) + feval(gD{1},cx(i) - hx/2,cy(i),cz(i))*k(ix,iy,iz)*2/hx^2;
%             fprintf('feval at (%d,%d,%d)\n',cx(i)-hx/2,cy(i),cz(i));
        elseif bnd(1) == 2
%             fprintf('Neumann Bou dary : ')
            b(i) = b(i) - feval(gN{1},cx(i) - hx/2,cy(i),cz(i))/hx;
%             fprintf('feval at (%d,%d,%d)\n',cx(i)-hx/2,cy(i),cz(i));
        else
            errid = 'MATLAB:funfun:ccfv3D:notABnd';
            errmsg = sprintf('Invalid value for boundary : value must be ''1'' or ''2''.');
            error(errid,errmsg);
        end
    else
        value = 2 / ((1/k(ix,iy,iz)+1/k(ix-1,iy,iz)) * hx^2);
        center = center + value;
        IA(l) = i;
        JA(l) = i-1;
        AA(l) = -value;
        l = l + 1;
%         fprintf('%d : ',JA(l-1));
%         fprintf('harmonic average of %d and %d\n',k(ix,iy,iz),k(ix-1,iy,iz));
    end
    
%     fprintf('[ Right ]\n')
    if ix == Nx
        if bnd(2) == 1 
%             fprintf('Dirichlet Boundary : ')
            center = center + k(ix,iy,iz)*2/hx^2;
            b(i) = b(i) + feval(gD{2},cx(i) + hx/2,cy(i),cz(i))*k(ix,iy,iz)*2/hx^2;
%             fprintf('feval at (%d,%d,%d)\n',cx(i)+hx/2,cy(i),cz(i));
        elseif bnd(2) == 2
%             fprintf('Neumann Bou dary : ')
            b(i) = b(i) - feval(gN{2},cx(i) + hx/2,cy(i),cz(i))/hx;
%             fprintf('feval at (%d,%d,%d)\n',cx(i)+hx/2,cy(i),cz(i));
        else
            errid = 'MATLAB:funfun:ccfv3D:notABnd';
            errmsg = sprintf('Invalid value for boundary : value must be ''1'' or ''2''.');
            error(errid,errmsg);
        end
    else
        value = 2 / ((1/k(ix,iy,iz)+1/k(ix+1,iy,iz)) * hx^2);
        center = center + value;
        IA(l) = i;
        JA(l) = i+1;
        AA(l) = -value;
        l = l + 1;
%         fprintf('%d : ',JA(l-1));
%         fprintf('harmonic average of %d and %d\n',k(ix,iy,iz),k(ix+1,iy,iz));
    end
    
%     fprintf('[ Front ]\n')
    if iy == 1
        if bnd(3) == 1 
%             fprintf('Dirichlet Boundary : ')
            center = center + k(ix,iy,iz)*2/hy^2;
            b(i) = b(i) + feval(gD{3},cx(i),cy(i) - hy/2,cz(i))*k(ix,iy,iz)*2/hy^2;
%             fprintf('feval at (%d,%d,%d)\n',cx(i),cy(i) - hy/2,cz(i));
        elseif bnd(3) == 2
%             fprintf('Neumann Bou dary : ')
            b(i) = b(i) - feval(gN{3},cx(i),cy(i) - hy/2,cz(i))/hy;
%             fprintf('feval at (%d,%d,%d)\n',cx(i),cy(i) - hy/2,cz(i));
        else
            errid = 'MATLAB:funfun:ccfv3D:notABnd';
            errmsg = sprintf('Invalid value for boundary : value must be ''1'' or ''2''.');
            error(errid,errmsg);
        end
    else
        value = 2 / ((1/k(ix,iy,iz)+1/k(ix,iy-1,iz)) * hy^2);
        center = center + value;
        IA(l) = i;
        JA(l) = i-Nx;
        AA(l) = -value;
        l = l + 1;
%         fprintf('%d : ',JA(l-1));
%         fprintf('harmonic average of %d and %d\n',k(ix,iy,iz),k(ix,iy-1,iz));
    end
    
%     fprintf('[ Back ]\n')
    if iy == Ny
        if bnd(4) == 1 
%             fprintf('Dirichlet Boundary : ')
            center = center + k(ix,iy,iz)*2/hy^2;
            b(i) = b(i) + feval(gD{4},cx(i),cy(i) + hy/2,cz(i))*k(ix,iy,iz)*2/hy^2;
%             fprintf('feval at (%d,%d,%d)\n',cx(i),cy(i) + hy/2,cz(i));
        elseif bnd(4) == 2
%             fprintf('Neumann Bou dary : ')
            b(i) = b(i) - feval(gN{4},cx(i),cy(i) + hy/2,cz(i))/hy;
%             fprintf('feval at (%d,%d,%d)\n',cx(i),cy(i) + hy/2,cz(i));
        else
            errid = 'MATLAB:funfun:ccfv3D:notABnd';
            errmsg = sprintf('Invalid value for boundary : value must be ''1'' or ''2''.');
            error(errid,errmsg);
        end
    else
        value = 2 / ((1/k(ix,iy,iz)+1/k(ix,iy+1,iz)) * hy^2);
        center = center + value;
        IA(l) = i;
        JA(l) = i+Nx;
        AA(l) = -value;
        l = l + 1;
%         fprintf('%d : ',JA(l-1));
%         fprintf('harmonic average of %d and %d\n',k(ix,iy,iz),k(ix,iy+1,iz));
    end
    
%     fprintf('[ Down ]\n')
    if iz == 1
        if bnd(5) == 1 
%             fprintf('Dirichlet Boundary : ')
            center = center + k(ix,iy,iz)*2/hz^2;
            b(i) = b(i) + feval(gD{5},cx(i),cy(i),cz(i)- hz/2)*k(ix,iy,iz)*2/hz^2;
%             fprintf('feval at (%d,%d,%d)\n',cx(i),cy(i),cz(i)- hz/2);
        elseif bnd(5) == 2
%             fprintf('Neumann Bou dary : ')
            b(i) = b(i) - feval(gN{5},cx(i),cy(i) ,cz(i)- hy/2)/hz;
%             fprintf('feval at (%d,%d,%d)\n',cx(i),cy(i),cz(i)- hz/2);
        else
            errid = 'MATLAB:funfun:ccfv3D:notABnd';
            errmsg = sprintf('Invalid value for boundary : value must be ''1'' or ''2''.');
            error(errid,errmsg);
        end
    else
        value = 2 / ((1/k(ix,iy,iz)+1/k(ix,iy,iz-1)) * hz^2);
        center = center + value;
        IA(l) = i;
        JA(l) = i-Nxy;
        AA(l) = -value;
        l = l + 1;
%         fprintf('%d : ',JA(l-1));
%         fprintf('harmonic average of %d and %d\n',k(ix,iy,iz),k(ix,iy,iz-1));
    end
    
%     fprintf('[ Up ]\n')
    if iz == Nz
        if bnd(6) == 1 
%             fprintf('Dirichlet Boundary : ')
            center = center + k(ix,iy,iz)*2/hz^2;
            b(i) = b(i) + feval(gD{6},cx(i),cy(i),cz(i)+ hz/2)*k(ix,iy,iz)*2/hz^2;
%             fprintf('feval at (%d,%d,%d)\n',cx(i),cy(i),cz(i)+ hz/2);
        elseif bnd(6) == 2
%             fprintf('Neumann Bou dary : ')
            b(i) = b(i) - feval(gN{6},cx(i),cy(i) ,cz(i)+ hy/2)/hz;
%             fprintf('feval at (%d,%d,%d)\n',cx(i),cy(i),cz(i)+ hz/2);
        else
            errid = 'MATLAB:funfun:ccfv3D:notABnd';
            errmsg = sprintf('Invalid value for boundary : value must be ''1'' or ''2''.');
            error(errid,errmsg);
        end
    else
        value = 2 / ((1/k(ix,iy,iz)+1/k(ix,iy,iz+1)) * hz^2);
        center = center + value;
        IA(l) = i;
        JA(l) = i+Nxy;
        AA(l) = -value;
        l = l + 1;
%         fprintf('%d : ',JA(l-1));
%         fprintf('harmonic average of %d and %d\n',k(ix,iy,iz),k(ix,iy,iz+1));
    end
    
    % Center
    IA(l) = i;
    JA(l) = i;
    AA(l) = center;
    l = l + 1;
end

A = sparse(IA(1:l-1),JA(1:l-1),AA(1:l-1),Nxyz,Nxyz);