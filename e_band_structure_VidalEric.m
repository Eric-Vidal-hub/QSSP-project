%% START DEFAULT COMMANDS OF *.m SCRIPTS
close all;       % Close all figures if any
clear all;       % Clear all variables/functions in memory
clc;             % Clear screen in the command window
source('mystartdefaults.m'); % Contains SI physical sonstants

tic

%% UNITS FOR FREE ELECTRONS IN VACUUM
recipunit = 1.0E+10;
ekinscale = ((hbar * recipunit)^2 / (2.0 * elm)) / qel;

%% COMPUTING THE ELECTRON BAND STRUCTURE OF FCC SEMICONDUCTORS
%% USING EMPIRICAL PSEUDOPOTENTIALS
dispersion_relation=true;
compute_dos=false;
% Materials: 1'Si';2'Ge';3'Sn';4'GaP';5'GaAs';6'AlSb';7'InP';8'GaSb';
% 9'InAs';10'InSb';11'ZnS';12'ZnSe';13'ZnTe';14'CdTe';15'Empty lattice'
semiconductor='CdTe'    % Choose the semiconductor (SC)
CohenBergstresser1966   % Load the parameters of the pseudopotential for SC
spacing=ls(m);          % Spacing in the reciprocal space for SC
ekinunit=ekinscale * (2*pi / spacing)^2; % Energy unit for SC

%% ORTHONORMAL BASIS VECTORS IN 3D RECIPROCAL SPACE
g=zeros(4, 3);
% Defininition according (eq 1.17), except 2pi factor:
g(1:3, 1)=cross(a(:,1), a(:,2)) / cell_volume;
g(1:3, 2)=cross(a(:,2), a(:,3)) / cell_volume;
g(1:3, 3)=cross(a(:,3), a(:,1)) / cell_volume;
for i=1:3
  g(4, :)=vecnorm(g(1:3, :));   % Squared norm of each vector
end

%% RECIPROCAL LATTICE VECTORS
min_norm = sqrt(min(g(4,:)));
nstep = ceil(cutoff^(1/4.) / min_norm); % With cutoff of exercise 3.8.7
nodes = (2 * nstep + 1)^3;  % Number of nodes in the grid
G=zeros(4, nodes);          % All possible G vectors in the grid
n=0;
for i=-nstep:nstep
  for j=-nstep:nstep
    for k=-nstep:nstep
      n=n+1;
      G(1:3,n)=i*g(1:3,1) + j*g(1:3,2) + k*g(1:3,3);
      G(4,n)=sum(G(1:3,n).^2);  % Squared norm of each vector
    end
  end
end
G=sortrows(G',4)';  % Sort G vectors by norm
G=G(:, G(4,:) <= cutoff);  % Remove G vectors with norm > cutoff
[~,Gcut]=size(G);  % Number of G vectors fulfilling the cutoff

%% COMPUTE THE PSEUDOPOTENTIALS
V_G=zeros(Gcut,1);
cases = [0, 3, 4, 8, 11];  % The cases you want to handle
for Gnum = find(ismember(G(4,:), cases))
  switch G(4,Gnum)
    case 0
      V_S=ff(m,1);
      V_A=0;
    case 3
      V_S=ff(m,2);
      V_A=ff(m,5);
    case 4
      V_S=0;
      V_A=ff(m,6);
    case 8
      V_S=ff(m,3);
      V_A=0;
    case 11
      V_S=ff(m,4);
      V_A=ff(m,7);
  end
  G_s=2*pi * dot(G(1:3,Gnum), tau(1:3,1));
  V_G(Gnum)=(cos(G_s)*V_S - 1i*sin(G_s)*V_A) * Rydberg;   % (eq 3.107)
end

%% COMPUTE THE HAMILTONIAN FULFILLING GS_MAX
% CONDITION: |G|^2 of highest non zero Fourier coefficients
% in expanding potential [2*pi/spacing]^2
Ham=zeros(Gcut,Gcut); % Defined by Hamiltonian matrix (eq 3.20)
G_aux=zeros(4,1);     % Auxiliar vector
for jj=1:Gcut
  for ii=1:Gcut
    G_aux(1:3)=G(1:3,ii)-G(1:3,jj);
    G_aux(4)=sum(G_aux(1:3).^2);
    if(G_aux(4)<=Gs_max)
      diff_norm = vecnorm(G_aux(1:3)-G(1:3,:));
      kk = find(abs(diff_norm) < tol);
      if ~isempty(kk)
        Ham(ii,jj)=V_G(kk(1));  % (eq 3.20)
      end
    end
  end
end

%% COMPUTE THE BZ PATH
[qpath,tix,til]=BZpath(BZstep,qs,qe,qs_str,qe_str);
[~,nqpath]=size(qpath);

%% COMPUTE THE BAND STRUCTURE
for qq=1:nqpath
  for ii=1:Gcut
    for jj=1:3  % qpath: Cartesian coordinates of wavevector i
      p(jj)=qpath(jj,qq)-G(jj,ii);
    end
    Ham(ii,ii)=ekinunit*dot(p,p)+V_G(1);  % Diagonal
  end
  [~,D]=eig(Ham);     % Diagonalization
  eigenvals=diag(D);  % Extract eigenvalues
  Eband(:,qq)=eigenvals(1:nband);
end
% Print the max value for the VB (sometimes band 4/5)
fprintf('\nVB maximum %.9f\n\n',Eband(4,27))

%% WRITE THE BAND STRUCTURE IN A FILE
filename=strcat(int2str(m), 'bandstructure.dat');
if(dispersion_relation)
  fid=fopen(filename,'w');
  fprintf(fid,'%d %d\n',nqpath,nband); % Write nqpath and nband at the beginning of the file
  for ii=1:nqpath
    for jj=1:nband
      fprintf(fid,'%f %f\n',qpath(5,ii),Eband(jj,ii));
    end
    fprintf(fid,'\n');
  end
  fclose(fid);
end

toc