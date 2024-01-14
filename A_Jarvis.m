function [vertex]=Jarvis(R,varargin)

%% JARVIS' MARCH ALGORITHM
%% FIND VERTICES OF THE CONVEX HULL OF A SET OF POINTS IN A PLANE
% 
% Ref: R.A. Jarvis, "On the identification of the convex hull 
%      of a finite set of points in the plane". 
%      Information Processing Letters vol 2, 18–21 (1973). 
%      doi:10.1016/0020-0190(73)90020-3
    
%% DEPENDENCE ON EXTERNAL FUNCTIONS
%  JarvisTurn.m
    
%% INPUT
%  
%  R = coordinates of n points in plane x-y  (2Xn array)
%
%  R(1,i) = x-coordinate of point i (i=1,...,n)
%  R(2,i) = y-coordinate of point i (i=1,...,n)
    
%% OUTPUT
%
%  vertex = coordinates of vertices of convex hull containing all
%           input points with the properties enabling plotting a 
%           closed hull as a polygon, i.e:
%             - the vertices are ordered anticlockwise;
%             - the last vertex = first vertex.
%  THEREFORE, THE CORRECT NUMBER OF VERTICES = columns(vertex)-1
%
%  vertex(1,i) = x-coordinate of vertex i  
%  vertex(2,i) = y-coordinate of vertex i  
    
%% Recognized options in varargin 
% (uppercases for readability are optional): 
    
% 'PlotPoints': plot input points
% 'PlotHull'  : plot polygon hull containing input points
% 'Verbose'   : display some computation details

%% DEFAUT VALUES OF OPTIONAL ARGUMENTS

plot_points=false; plot_hull=false; verbose=false;

%% PARSE OPTIONAL ARGUMENT LIST

name_value_pair=false;
for k = 1:length(varargin);
    if (name_value_pair)
        name_value_pair=false;
    else
        switch lower(varargin{k}) % varargin is a "cell array" 
                                  % => array index within braces
          case {'plotpoints'}
            plot_points=true; S=R; % Need a backup of R for plotting
          case {'plothull'}
            plot_hull=true;   S=R; % Need a backup of R for plotting
          case {'verbose'}
            verbose=true;
          otherwise
            error(['function Jarvis: ',...
                   'option %s not recognized.\n'],...
                  varargin{k});
        end
    end
end

%% PRECISION

tol = 1e-12; % Two floating-point numbers will be considered equal 
             % if the absolute value of their difference is < tol
    
%% CORE JOB
 
[~,n]=size(R);

if (n<3)
    error(['function Jarvis: there must be a least 3 ' ...
           'points.\n']);
end
x_min=floor(min(R(1,:)));x_max=ceil(max(R(1,:))); % Stored for plotting purpose
y_min=floor(min(R(2,:)));y_max=ceil(max(R(2,:))); % Stored for plotting purpose

[ymin,lowest] = min(R(2,:)); % Ordinate and index of lowest point

vertex=R(:,lowest); % Lowest point = first vertex

R = [R(:,1:lowest-1),R(:,lowest+1:n)]; % Remove first vertex from
                                       % list where to search the
                                       % second vertex

step=1;  a=vertex;  target=vertex; % Initialization
unclosed_loop=true;

while (unclosed_loop); % Looping until back to starting point
    
    b=R(:,1); k=1; % Start with first point of the updated search list
    
    [~,n]=size(R);
    for i=2:n % Search the point that is most on the right
        c=R(:,i);
        turn=JarvisTurn(a,b,c);
        switch turn
        % case -2 % Points b and c colinear in opposite
                  % directions can not happen in this algorithm 
                  % => No action.
          case -1 % If point c is on the right of point b,
            buf=b; b=c; c=buf; % permute b and c
            k=i; 
          case 0  % Identical points are duplicated in list of vertices
            buf=b; b=c; c=buf; % permute b and c
            k=i; 
        % case 1  % Keep b as possible vertex until end of for <=> No action.
          case 2  % Points b and c colinear in same direction
            if (norm(b-a) <= norm(c-a)) % If c most distant from a,
                buf=b; b=c; c=buf;      % permute b and c.
                k=i;
            end  
        end     
    end
    vertex=[vertex b]; % Next vertex found = last value of b
                       % = point most on the right
   
    a=b; % Last found vertex becomes the next pivot 
    
    if(norm(a-target) < tol)
        unclosed_loop=false; % Closed loop => end of while
    end
    
    R(:,k)=[]; % Remove the last found vertex from the search list
    if(step == 1)            % After finding the second vertex,
        R=[ vertex(:,1) R ]; % reintroduce first vertex as first
                             % element of the search
                             % list to enable to find it again 
                             % as target terminating the while loop.
    end
    step=step+1; % Increment the number of hull vertices
    
end % End of while (unclosed_loop)

if (verbose)
    % The number of vertices = step-1 since the last vertex = first vertex
    fprintf('Function Jarvis: convex hull = %d vertex polygon.\n',step-1)
end

%% OPTIONAL PLOTS

if (plot_points || plot_hull)

    h = findall(0,'type','figure');
    if isempty(h)    
        figure('NumberTitle', 'off',...
               'name','2D convex hull by the method of Jarvis');
    end
    
    if(plot_points)
        for j=1:n
            plot(S(1,j),S(2, j),'or'); hold on; % Use backup S for
                                                % plotting because
                                                % R was mofified by algorithm
        end
    end
    
    if (plot_hull)
        plot(vertex(1,:),vertex(2,:),'-r'); hold on;
    end
    
    xlim([x_min,x_max]); ylim([y_min,y_max]);
    axis equal; xlabel('x'); ylabel('y');
end

end % End of function Jarvis



