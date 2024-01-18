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
CohenBergstresser1966
spacing=ls(m);  % Spacing between k-points in the reciprocal space for SC m
ekinunit=ekinscale * (2*pi / spacing)^2; % Energy unit for SC m

%% ORTHONORMAL BASIS VECTORS IN 3D RECIPROCAL SPACE
g=zeros(4, 3);
% Defininition according (eq 1.17), except 2pi factor:
g(1:3, 1)=cross(a(:,1), a(:,2)) / cell_volume;
g(1:3, 2)=cross(a(:,2), a(:,3)) / cell_volume;
g(1:3, 3)=cross(a(:,3), a(:,1)) / cell_volume;
for i=1:3
    g(4, i)=g(1:3, i)' * g(1:3, i);   % Squared norm of each vector
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
      G(4,n)=G(1:3,n)' * G(1:3,n);  % Squared norm of each vector
    end
  end
end
G=sortrows(G',4)';  % Sort G vectors by norm
G=G(:, G(4,:) <= cutoff);  % Remove G vectors with norm > cutoff
[~,Gcut]=size(G);  % Number of G vectors fulfilling the cutoff

%% COMPUTE THE PSEUDOPOTENTIALS
cvg=zeros(Gcut,1);
cases = [0, 3, 4, 8, 11];  % The cases you want to handle
for Gnum = find(ismember(G(4,:), cases))
  switch G(4,Gnum)
    case 0
      sym=ff(m,1);
      asym=0;
    case 3
      sym=ff(m,2);
      asym=ff(m,5);
    case 4
      sym=0;
      asym=ff(m,6);
    case 8
      sym=ff(m,3);
      asym=0;
    case 11
      sym=ff(m,4);
      asym=ff(m,7);
  end
  argu=2*pi * (G(1:3,Gnum)' * tau(1:3,1));
  cvg(Gnum)=(cos(argu)*sym - 1i*sin(argu)*asym) * Rydberg;
end

toc

%% COMPUTE THE BZ PATH AND HAMILTONIAN
[q,tix,til]=BZpath(BZstep,qs,qe,qs_str,qe_str);
[~,nq]=size(q);

H=zeros(Gcut,Gcut);
G_diff=zeros(5,1);
for jj=1:Gcut
  for ii=1:Gcut
    G_diff(1:3)=G(1:3,ii)-G(1:3,jj);
    G_diff(5)=G_diff(1:3)'*G_diff(1:3);
    if(G_diff(5)<=Gs_max)
      diff_norm = vecnorm(G_diff(1:3)-G(1:3,:), 2, 1);
      kk = find(abs(diff_norm) < tol);
      if ~isempty(kk)
        H(ii,jj)=cvg(kk(1));
      end
    end
  end
end

for iq=1:nq
  testo(iq)=iq;
  %%Kinetic energy
  for ii=1:Gcut
    for jj=1:3
      p(jj)=q(jj,iq)-G(jj,ii);
    end
    H(ii,ii)=ekinunit*(p*p')+cvg(1);  % Hamiltonian diagonal
  end
  [v,ev]=eig(H);
  E=real(diag(ev));
  [E,perm]=sort(E);
  v=v(:,perm);
  for ii=1:nband
        firstEs(ii,iq)=E(ii);
  end

end


for ii=1:nq
  separ(ii)=20*((-1)^q(6,ii)); %Function used to later put a line to separate paths
end
% Print the last value until 7th decimal
fprintf('%.9f\n',firstEs(4,27))
plot(q(5,:),firstEs,'-o','color', 'black','MarkerSize', 1, 'linewidth', 0.5,q(5,:),separ,'color','red');
ylabel('E(eV)','FontSize',18);
xlim([0,q(5,nq)]);
if(m==1)
  ylim([-5.5,6]);
  title ("Band structure for Si", "fontsize", 20);
  set(gca, 'ytick', -5:1:6);
elseif(m==2)
  ylim([-5,7]);
  title ("Band structure for Ge", "fontsize", 20);
  set(gca, 'ytick', -5:1:7);
elseif(m==3)
  ylim([-4,6]);
  title ("Band structure for Sn", "fontsize", 20);
  set(gca, 'ytick', -4:1:6);
elseif(m==4)
  ylim([-4,7]);
  title ("Band structure for GaP", "fontsize", 20);
  set(gca, 'ytick', -4:1:7);
elseif(m==5)
  ylim([-4,7]);
  title ("Band structure for GaAs", "fontsize", 20);
  set(gca, 'ytick', -4:1:7);
elseif(m==6)
  ylim([-4,7]);
  title ("Band structure for AlSb", "fontsize", 20);
  set(gca, 'ytick', -4:1:7);
elseif(m==7)
  ylim([-4,7]);
  title ("Band structure for InP", "fontsize", 20);
  set(gca, 'ytick', -4:1:7);
elseif(m==8)
  ylim([-3,6.4]);
  title ("Band structure for GaSb", "fontsize", 20);
  set(gca, 'ytick', -3:1:6);
elseif(m==9)
  ylim([-4,7]);
  title ("Band structure for InAs", "fontsize", 20);
  set(gca, 'ytick', -4:1:7);
elseif(m==10)
  ylim([-3,6]);
  title ("Band structure for InSb", "fontsize", 20);
  set(gca, 'ytick', -3:1:6);
elseif(m==11)
  ylim([-3,10]);
  title ("Band structure for ZnS", "fontsize", 20);
  set(gca, 'ytick', -3:1:10);
elseif(m==12)
  ylim([-3,9]);
  title ("Band structure for ZnSe", "fontsize", 20);
  set(gca, 'ytick', -3:1:9);
elseif(m==13)
  ylim([-3,8.5]);
  title ("Band structure for ZnTe", "fontsize", 20);
  set(gca, 'ytick', -3:1:8);
elseif(m==14)
  ylim([-3.5,8]);
  title ("Band structure for CdTe", "fontsize", 20);
  set(gca, 'ytick', -3:1:8);
end

set (gca,'xtick',tix);
set (gca,'xticklabel',til,'FontSize',18,'FontWeight','bold');

