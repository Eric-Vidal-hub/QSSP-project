%% START DEFAULT COMMANDS OF *.m SCRIPTS
close all;       % Close all figures if any
clear all;       % Clear all variables/functions in memory
clc;             % Clear screen in the command window
source('mystartdefaults.m'); % Contains SI physical sonstants

tic

%% UNITS FOR FREE ELECTRONS IN VACUUM
recipunit = 1.0E+10;
ekinscale = ((hbar * recipunit)^2 / (2.0 * elm)) / qel

%% COMPUTING THE ELECTRON BAND STRUCTURE OF FCC SEMICONDUCTORS
%% USING EMPIRICAL PSEUDOPOTENTIALS
dispersion_relation=true;
compute_dos=false;
CohenBergstresser1966

%% ORTHONORMAL BASIS VECTORS IN 3D RECIPROCAL SPACE
g=zeros(4, 3);
% Defininition according (eq 1.17)
g(1:3, 1)=cross(a(:,1), a(:,2)) / cell_volume;
g(1:3, 2)=cross(a(:,2), a(:,3)) / cell_volume;
g(1:3, 3)=cross(a(:,3), a(:,1)) / cell_volume;
% and square for each one
for i=1:3
    g(4, i)=g(1:3, i)' * g(1:3, i);
end

%% RECIPROCAL LATTICE VECTORS
min_norm =sqrt(min(g(4,:)));
nstep= ceil(sqrt(sqrt(cutoff))/min_norm);
nodes=(2*nstep+1)^3;
G=zeros(5,nodes);
n=0;
for ii=-nstep:nstep
  for jj=-nstep:nstep
    for kk=-nstep:nstep
      n=n+1;
      G(1:3,n)=ii*g(1:3,1)+jj*g(1:3,2)+kk*g(1:3,3);
      G(5,n)=G(1:3,n)'*G(1:3,n);
      G(4,n)=sqrt(G(5,n));
    end
  end
end

GT=sortrows(G',4);
G=GT';


toc