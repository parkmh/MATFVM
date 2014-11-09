% For more information, see <a href="matlab: 
% web('http://www.grandmaster.colorado.edu/~parkmh')">Minho Park's Web site</a>.
cd bin
Ver = fvmversion;
user = 'Minho Park';
email = 'min.park@nottingham.ac.uk';
cd ..

clc
fprintf(' ********************************************\n')
fprintf('\n')
fprintf(' Matlab Finite Volume Toolbox  %s \n',Ver) 
fprintf('\n')
fprintf('%20s %s\n','Written by', user);
fprintf('%13s %s\n','email :',email);
fprintf(' ********************************************\n')


cwd = pwd;
matfvmroot = pwd;

% Add path
% Generate amgpath.m file
fid = fopen([pwd filesep 'bin' filesep 'fvmpath.m'],'w');
fprintf(fid,'function matfvmpath = fvmpath\n');

fprintf('\n1. Add path\n%s\n',matfvmroot)

addpath(fullfile(matfvmroot,'bin'));
fprintf(fid,'matfvmpath = ''%s'';\n',matfvmroot);

fclose(fid);
savepath

clear all

