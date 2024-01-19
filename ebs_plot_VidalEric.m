%% START DEFAULT COMMANDS OF *.m SCRIPTS
close all;       % Close all figures if any
clear all;       % Clear all variables/functions in memory
clc;             % Clear screen in the command window

% READ THE BAND STRUCTURE FROM A FILE
m=14; % Choose the material
bs_step=4; % Choose the step for the Band Structure plot
filename = strcat(int2str(m),'bandstructure.dat');
fid=fopen(filename,'r');
nqpath=fscanf(fid,'%d',1);  % Read nqpath from the file
nband=fscanf(fid,'%d',1);   % Read nband from the file
Eband=zeros(nband,nqpath);  % Initialize Eband with the correct size
for ii=1:nqpath
  for jj=1:nband
    qpath(5,ii)=fscanf(fid,'%f',1);
    Eband(jj,ii)=fscanf(fid,'%f',1);
  end
  fscanf(fid,'%c',1);
end
fclose(fid);

%% PLOT THE BAND STRUCTURE
% Define the titles and y-limits for each material
titles = {"Band structure for Si", "Band structure for Ge", "Band structure for Sn", "Band structure for GaP", "Band structure for GaAs", "Band structure for AlSb", "Band structure for InP", "Band structure for GaSb", "Band structure for InAs", "Band structure for InSb", "Band structure for ZnS", "Band structure for ZnSe", "Band structure for ZnTe", "Band structure for CdTe"};
ylimits = [[-5.5,6]; [-5,7]; [-4,6]; [-4,7]; [-4,7]; [-4,7]; [-4,7]; [-3,6.4]; [-4,7]; [-3,6]; [-3,10]; [-3,9]; [-3,8.5]; [-3.5,8]];

% Plot the band structure
plot(qpath(5,1:bs_step:end),Eband(:,1:bs_step:end),'-o','color', 'black','MarkerSize', 2, 'linewidth', 0.5, 'MarkerFaceColor', 'white');
ylabel('E(eV)','FontSize',18);
xlim([0,qpath(5,nqpath)]);
ylim(ylimits(m,:));
title (titles{m}, "fontsize", 20);
set(gca, 'ytick', ylimits(m,1):1:ylimits(m,2));

set (gca,'xtick',tix);
set (gca,'xticklabel',til,'FontSize',18,'FontWeight','bold');