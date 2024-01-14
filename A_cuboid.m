function cuboid(a,varargin)

%%  3-DIMENSIONAL PLOT OF THE CUBOID HEXAHEDRON
%%  DEFINED BY BASIS VECTORS

%%  Input [unit of length]:

%  a = 3X3 array of basis vectors in column order
%      i.e: a(1:3,i) = basis vector i 
%      in any unit of length.
    
%% Output
    
%  3-dimensional plot of cuboid hexahedron in current figure
    
%% Recognized options in varargin 
% (uppercases for readability are optional): 

% if varargin{k} = 'Origin', 
%                   then varargin{k+1} = [x,y,z]' 
%                   = coordinates of origin of the basis vectors
% if varargin{k} = 'DisplayName', 
%                   then varargin{k+1} = Character string 
%                   that will appear in legend box
% if varargin{k} = 'ColorMap', 
%                   then varargin{k+1} = Character string 
%                   defining one of the pre-defined colormaps
% if varargin{k} = 'ShowLegend', 
%                   then varargin{k+1} = Boolean to show legend or not
% List of Octave pre-defined colormaps:
% https://octave.sourceforge.io/octave/function/colormap.html
    
%% DEFAUT VALUES OF OPTIONAL ARGUMENTS

default_color_map=true; show_legend=true;
origin=[0 0 0]';
displayname='Cuboid spanned by basis vectors';
    
%% PARSE OPTIONAL ARGUMENT LIST

name_value_pair=false;
for k = 1:length(varargin);
    if (name_value_pair)
        name_value_pair=false;
    else
        switch lower(varargin{k}) % varargin is a "cell array" 
          case {'origin'}
            origin=varargin{k+1}; name_value_pair=true;
          case {'displayname'}
            displayname=varargin{k+1}; name_value_pair=true;
          case {'colormap'}
            color_map=varargin{k+1}; default_color_map=false;
            name_value_pair=true;
          otherwise
            fprintf(['error in function cuboid: ',...
                     'option %s not recognized.\n'],...
                    varargin{k}); return;
        end
    end
end

%% CORE JOB

h = findall(0,'type','figure');
if isempty(h)
    fig=figure('NumberTitle', 'off','name','Cuboid hexahedron');
    k=0; % Legend item counter
else
    hleg=legend(); % Caution: requires at least one graphical object
                   % inside the current figure window created by the
                   % calling program (at least using initialize_legend3D.m)
    if isempty(hleg)
        k=0;
    else
        Legend = get(legend(),'string'); % Retrieve legend of current figure
        [k,~]=size(Legend); % Legend item counter
    end
end
Legend{k+1}=displayname;

if (default_color_map) 
    cmap=colormap('default');
else
    cmap=colormap(color_map);
end

vertex(1:3,1)  = origin(1:3);
vertex(1:3,2)  = origin(1:3) + a(1:3,1);
vertex(1:3,3)  = origin(1:3) + a(1:3,1) + a(1:3,2);
vertex(1:3,4)  = origin(1:3) + a(1:3,2);
vertex(1:3,5)  = origin(1:3) + a(1:3,3);
vertex(1:3,6)  = origin(1:3) + a(1:3,1) + a(1:3,3);
vertex(1:3,7)  = origin(1:3) + a(1:3,1) + a(1:3,2) + a(1:3,3);
vertex(1:3,8)  = origin(1:3) + a(1:3,2) + a(1:3,3);

face = [ 1, 2, 3, 4;
         5, 6, 7, 8;
         1, 2, 6, 5;
         4, 3, 7, 8;
         2, 3, 7, 6;
         4, 1, 5, 8 ];

%nw=10; cmap(end-nw+1:end,:)=[]; cmap(1:nw-1,:)=[]; % Avoid extreme colors
%step=rows(cmap)/(rows(face)-1); 
%color_scheme(1,1:3)=cmap(1,1:3);
%for i=2:rows(face)
%    color_scheme(i,1:3)=cmap(floor((i-1)*step),1:3);
%end

    
offset=rows(cmap)-10; % Avoid extreme colors when changing
                      % face. Nice when using the
                      % following corlormaps: viridis,
                      % gray,bone,pink,autumn,ocean,winter
for i=1:rows(face)
    color_scheme(i,1:3)=cmap(offset-i,1:3);
end

patch ('Vertices', vertex', 'Faces', face, ...
       'LineStyle','-',...
       'LineWidth', 1,...
       'EdgeColor',[0 0 0],...% Options: 'none', 'flat', 'interp', [R G B]
       'EdgeAlpha',1,...      % Edge transparency in range 0 (transparent) to 1 (opaque)
       'FaceColor', 'flat',...% Options: 'none', 'flat', 'interp', [R G B]
       'FaceAlpha',0.2,...   % Face transparency in range 0 (transparent) to 1 (opaque)
       'FaceVertexCData', color_scheme,...
       'DisplayName',Legend{k+1}
); hold on;

% Axes, Labels & Legend
    
axis equal; grid on; legend(Legend,'Location','EastOutside');

if(~show_legend)
    legend(gca,'off');
end

end % End of function cuboid
    