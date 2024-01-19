function [u,kpt,weight]=MonkhorstPack(g,n,varargin)

%%  MONKHORST-PACK's SAMPLING OF THE BRILLOUIN ZONE
%  
%  Reference:
%  [1] H.J. Monkhorst & J.D. Pack, Phys. Rev. B13, 5188-5192(1976).
          
%% Input
%  g = 3X3 array of reciprocal lattice basis vectors in column
%  order, i.e: g(1:3,i) = reciprocal basis vector i
%  in any unit.
%
%  n = integer q in equation (3) of reference [1]
    
%% Output
%  
%  u = n X 1 array such that u(i) = (2*i-n-1)/(2*n) [dimensionless]
%    = equation (3) of reference [1]
%
%  kpt = 3 X (n^3) array of k-points sampling the Brillouin Zone
%        in same unit as g    
%      = equation (4) of reference [1]
%  weight = n X 1 array, w(i)= weight of kpt(1:3,i)
    
%% Recognized options in varargin 
% (uppercases for readability are optional): 
  
% 'PlotRecipBasisVectors': plot basis vectors of reciprocal lattice
% 'PlotSampling': plot sampling points 
% 'PlotPrimitiveCell'  to plot the Brillouin Zone cuboid  
%                      supported by the reciprocal basis vectors  
%                      and centered at [0 0 0]
% if varargin{k} = 'Symmetry',
%                  then varargin{k+1} = Guideline for exploiting
%                  symmetry. Allowed values: 
%                  1 default, no symmetry used, weight of each
%                    point = 1/n^3
%                  8 to exploit 8-fold symmetry so that the grid is
%                    reduced to one octant of the BZ, each point
%                    being given a weight=8/n^3.
    
% if varargin{k} = 'File', 
%                  then varargin{k+1} = Character string 
%                                     = name of output file
% if varargin{k} = 'Color', then varargin{k+1 } = [R G B] color of
%                   the sampling points
% if varargin{k} = 'ColorMap', 
%                   then varargin{k+1} = Character string 
%                   defining one of the pre-defined colormaps
%                   used to plot the primitive cell's faces
% if varargin{k} = 'ShowLegend', 
%                   then varargin{k+1} = Boolean to show legend or not
    
%% DEFAUT VALUES OF OPTIONAL ARGUMENTS

color_kpt=[0.6350,0.0780,0.1840];
color_map='default'; show_legend=true;
plot_sampling = false; plot_primitive_cell = false;
plot_recip_basis_vectors=false;
fileID=[];
Title='Monkhorst-Pack sampling of the Brillouin zone';
    
%% PARSE OPTIONAL ARGUMENT LIST

name_value_pair=false;
for k=1:length(varargin);
    if (name_value_pair)
        name_value_pair=false;
    else
        switch lower(varargin{k}) % varargin is a "cell array"
          case {'plotrecipbasisvectors'}
            plot_recip_basis_vectors  = true;
          case {'plotsampling'}
            plot_sampling = true;
          case {'plotprimitivecell'}
            plot_primitive_cell = true; 
          case {'symmetry'}
              name_value_pair=true; foldsym=varargin{k+1};
          case {'file'}
            name_value_pair=true;
            fileID=varargin{k+1};
          case {'color'}
            color_kpt = varargin{k+1}; name_value_pair=true;
          case {'colormap'}
            color_map = varargin{k+1}; name_value_pair=true;
          case {'showlegend'}
            show_legend=varargin{k+1}; name_value_pair=true;
          otherwise
            error(['MonkhorstPack: ',...
                   'option %s not recognized.\n'],...
                  varargin{k});
        end
    end
end

%% PRECISION
tol = 1e-12; % Two floating-point numbers will be considered equal 
             % if the absolute value of their difference is < tol

%% CORE JOB

if (n<=0)
    error('Input error in function MonkhorstPack: n <=0')
end


if (foldsym==8 && rem(n,2) ~= 0)
    n=n-1;
    fprintf(['function MonkhorstPack: number of points in each ' ...
             'reciprocal basis vector direction changed to %d\n'],n)
end

for i=1:n
    u(i) = (2*i-n-1)/(2*n);
end


for i=1:3 % Adjusting upper limits of loops allows 1D and 2D 
          % sampling without duplicating points.
    if (norm(g(:,i)) < tol) % True for some direction(s) 
                            % in 1D and 2D lattices
        nn(i)=1;            % Avoid sampling along g(:,i)
    else
        nn(i)=n;            % Sampling along g(:,i)
    end
end

switch foldsym
  case 1 
    fprintf('Sampling of BZ not exploiting any symmetry\n')
    m=0;
    for i=1:nn(1)
        for j=1:nn(2)
            for k=1:nn(3)
                m=m+1;
                kpt(:,m)=u(i)*g(:,1)+u(j)*g(:,2)+u(k)*g(:,3);
            end
        end
    end
    weight(1:m)=1/m;
  case 8
    fprintf('Sampling of BZ exploiting 8-fold symmetry\n')
    m=0;
    for i=1:nn(1)
        for j=1:nn(2)
            for k=1:nn(3)
                if(u(i)>0 && u(j)>0 && u(k)>0)
                    m=m+1;
                    kpt(:,m)=u(i)*g(:,1)+u(j)*g(:,2)+u(k)*g(:,3);
                end
            end
        end
    end
    weight(1:m)=8/m;
  otherwise
    error('function MonkhorstPack: symmetry =%d not implemented',foldsym)
end

%% Formatted output

if(~isempty(fileID))
    printtable(kpt','LineNumber',true,'ArrayName','q',...
               'file',fileID,'Title',Title)
end

%% OPTIONAL PLOT

if (plot_sampling || plot_primitive_cell)
    
    h = findall(0,'type','figure');
    if isempty(h)
       fig=figure('NumberTitle', 'off','name',Title);
    end
    
    if (plot_recip_basis_vectors) 
        origin = [ 0 0 0 ]';
        color_scheme = get(gca,'colororder');
        gvec{1}='g_1';gvec{2}='g_2';gvec{3}='g_3';        
        for i=1:3
            k=k+1; Legend{k}=gvec{i};
            quiver3(origin(1),origin(2),origin(3),g(1,i),g(2,i),g(3,i),...
                    'autoscale','off','color',color_scheme(i+3,1:3),...
                    'maxheadsize',0.1,...
                    'DisplayName',Legend{k}); hold on;
        end
    end
    
    if (plot_sampling)
        plot3(kpt(1,:),kpt(2,:),kpt(3,:),...
              'o','color',color_kpt,...
              'DisplayName','Monkorst-Pack grid');
        hold on;        
    end
    
    if (plot_primitive_cell) 
        origin = -(g(1:3,1) + g(1:3,2) +  g(1:3,3))/2;
        cuboid(g,'Origin',origin,'DisplayName','Monkhorst-Pack cell',...
               'ColorMap',color_map)
    end
    
 
    
    axis equal; grid on; legend('Location','EastOutside');
    xlabel('q_1'); ylabel('q_2'); zlabel('q_3');
    
    if(~show_legend)
        legend(gca,'off');
    end
end
    
end % End of function MonkhorstPack
