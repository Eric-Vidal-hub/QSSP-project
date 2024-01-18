clear all
close all
clc
mystartdefaults
%%%%%%%%%% Pseudopotentials%%%%%%%%%%%%%%
dispersion_relation=true; compute_dos=false;
CohenBergstresser1966
recipunit=1e+10;
ekinscale=(hbar*recipunit)^2/(2*elm)/qel;
%%%%%%%%%%%%% Unit vectors in reciprocal space %%%%%%%%%%%%%%%

% 2pi missing in definition, and, <'* == .*>?

g=zeros(4,3);
g(1:3,1)=cross(a(:,2),a(:,3))/cell_volume;
g(1:3,2)=cross(a(:,3),a(:,1))/cell_volume;
g(1:3,3)=cross(a(:,1),a(:,2))/cell_volume;
for ii=1:3
  g(4,ii)=g(1:3,ii)' * g(1:3,ii);
end
%%%%%%%%%Build reciprocal lattice%%%%%%%%%%
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

%%%%%%%% only vectors inside cutoff %%%%%%%%%
kept=1;
for nn=2:nodes
  if(G(5,nn)<=cutoff)
    kept=kept+1;
  end
end



%%%%%%%%%%%%%%%%%%% pseudopotential %%%%%%%%%%%%%%%%%%%%%%%%%%

spacing=ls(m);
ekinunit=ekinscale*(2*pi/spacing)^2;

for nn=1:kept
  sym=0;
  asym=0;
  cvg(nn)=0; %%%%%% its zero unless its norm^2 is smaller than 11, which is checked right after
  if(G(5,nn)==0)
    sym=ff(m,1)*Rydberg;
    asym=0;
  end
  if(G(5,nn)==3)
    sym=ff(m,2)*Rydberg;
    asym=ff(m,5)*Rydberg;
  end
  if(G(5,nn)==4)
    sym=0;
    asym=ff(m,6)*Rydberg;
  end
  if(G(5,nn)==8)
    sym=ff(m,3)*Rydberg;
    asym=0;
  end
  if(G(5,nn)==11)
    sym=ff(m,4)*Rydberg;
    asym=ff(m,7)*Rydberg;
  end
  argu=2*pi*(G(1:3,nn)'*tau(1:3,1));
  cvg(nn)=cos(argu)*sym-1i*sin(argu)*asym;
end
%%%%%%%%%%Loop over wavevectors

[q,tix,til]=BZpath(BZstep,qs,qe,qs_str,qe_str);
[~,nq]=size(q);
for ii=1:nq
  separ(ii)=20*((-1)^q(6,ii)); %Function used to later put a line to separate paths
end

H=zeros(kept,kept);
G_diff=zeros(5,1);
for jj=1:kept
  for ii=1:kept
    G_diff(1:3)=G(1:3,ii)-G(1:3,jj);
    G_diff(5)=G_diff(1:3)'*G_diff(1:3);
    if(G_diff(5)<=Gs_max)
      for kk=1:kept
        if(norm(G_diff(1:3)-G(1:3,kk))<tol)
          H(ii,jj)=cvg(kk);
        end
      end
    end
  end
end


%shifts=[0,-9.48,-6.815,-9.238,-8.877,-6.978,-7.676,-6.964,-7.153,-6.073,-6.372,-6.128,-5.347,-4.272]; %shifts applied to the band structures so that the valence band is at 0


for iq=1:nq
  testo(iq)=iq;
  %%Kinetic energy
  for ii=1:kept
    for jj=1:3
      p(jj)=q(jj,iq)-G(jj,ii);
    end
    H(ii,ii)=ekinunit*(p*p')+cvg(1);  %+shifts(m); %diagonal + shift
  end
  [v,ev]=eig(H);
  E=real(diag(ev));
  [E,perm]=sort(E);
  v=v(:,perm);
  for ii=1:nband
        firstEs(ii,iq)=E(ii);
  end

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

