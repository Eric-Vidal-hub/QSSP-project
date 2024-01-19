%% START DEFAULT COMMANDS OF *.m SCRIPTS
close all;       % Close all figures if any
clear all;       % Clear all variables/functions in memory
clc;             % Clear screen in the command window

%% DEFINE THE BAND STRUCTURE PLOT PARAMETERS
semiconductor='Ge';   % Choose the semiconductor
bs_step=4;              % Choose the step for the Band Structure plot

%% READ THE BAND STRUCTURE FROM A FILE
materials = {'Si','Ge','Sn','GaP','GaAs','AlSb','InP','GaSb', ...
'InAs','InSb','ZnS','ZnSe','ZnTe','CdTe','Empty lattice'};
m = find(strcmp(materials,semiconductor))
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
tix = fscanf(fid, '%f', [1, inf]); % Read tix
til = {}; % Initialize til as an empty cell array
line = fgetl(fid);  % Read the first line
index = 1;          % Initialize index
while ischar(line)  % Continue until fgetl returns -1
  if ~isempty(line) % If the line is not empty
    if mod(index, 2) == 0 % If index is even
      line = sprintf('\n\n%s', line); % Add two new lines before the line
    end
    til{end+1} = line;  % Add the line to til
    index = index + 1;  % Increment index
  end
  line = fgetl(fid);    % Read the next line
end
fclose(fid);

%% PLOT THE BAND STRUCTURE
ylimits = [[-5.5,6]; [-5,7]; [-4,6]; [-4,7]; [-4,7]; [-4,7]; [-4,7];...
[-3,6.4]; [-4,7]; [-3,6]; [-3,10]; [-3,9]; [-3,8.5]; [-3.5,8]];

plot(qpath(5,1:bs_step:end),Eband(:,1:bs_step:end),'-o','color',...
'black','MarkerSize', 2, 'linewidth', 0.5, 'MarkerFaceColor', 'white');
ylabel('E   (eV)','FontSize',18);
xlim([0,qpath(5,nqpath)]);
ylim(ylimits(m,:));
title (materials{m}, "fontsize", 20);
set(gca, 'ytick', ylimits(m,1):1:ylimits(m,2));

set (gca,'xtick',tix);
set (gca,'xticklabel',til,'FontSize',18,'FontWeight','bold');