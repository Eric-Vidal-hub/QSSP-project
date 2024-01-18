
%% READ THE BAND STRUCTURE FROM A FILE
m=14; % Choose the material
filename = strcat(int2str(m),'bandstructure.dat');
fid=fopen(filename,'r');
for ii=1:nqpath
  for jj=1:nband
    Eband(jj,ii)=fscanf(fid,'%f',1);
    Eband(jj,ii)=fscanf(fid,'%f',1);
  end
  fscanf(fid,'%c',1);
end
fclose(fid);


%% PLOT THE BAND STRUCTURE
for ii=1:nqpath
  separ(ii)=20*((-1)^qpath(6,ii)); %Function used to later put a line to separate paths
end

plot(qpath(5,:),Eband,'-o','color', 'black','MarkerSize', 1, 'linewidth', 0.5,qpath(5,:),separ,'color','red');
ylabel('E(eV)','FontSize',18);
xlim([0,qpath(5,nqpath)]);
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