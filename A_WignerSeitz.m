function [vertex,edge,face,vn]=WignerSeitz(a,varargin)
    
%% COMPUTE AND PLOT WIGNER-SEITZ PRIMITIVE CELL
%  Wigner-Seitz Cell = primitive cell containing a single lattice node.
%  In reciprocal space, Wigner-Seitz cell = Brillouin Zone.
    
%% Dependence on external functions
%  Jarvis.m
%  JarvisTurn.m
%  lattice.m
%  CAUTION: For a correct management of the legends, if the plot
%           options should be added to a figure already opened by
%           the calling program, at least one graphical object
%           inside the current figure window must have been created by the
%           calling program, for example by calling program 
%           the function initialize_legend3D.m.
    
%% Inputs:
%  a = 3X3 array of lattice basis vectors in column order
%      and in arbitrary unit of length
%
%      a(:,1) = lattice basis vector 1,
%      a(:,2) = lattice basis vector 2,
%      a(:,3) = lattice basis vector 3.
%
    
%% Outputs
%  vertex = 3XM array = M vertices of Wigner-Seitz cell
%  vertex(1:3,i) = Cartesian coordinates of vertex i
%
%  edge   = 2XN array = N edges of Wigner-Seitz cell
%  edge(1,n) = i <=> n-th edge starts at vertex(:,i)
%  edge(2,n) = j <=> n-th edge   ends at vertex(:,j)
%
%  face = NF*MV array where NF = number of faces 
%         of Wigner-Seitz cell and MV = maximum of the numbers of
%         vertices of each face and such that 
%
%         vertex(1:3,face(i,j)) = vertex i of face j
%
%  NB: if face(i,j) = NaN, this means that the number of
%      vertices of face j is smaller than i
%      (feature needed for plotting: see Matlab/Octave manual)
%  
%  vn = 5XNF array = vectors normal to faces that are closest to
%       the origin
%  vn(1:3,i) Cartesian coordinates of vector normal to face i
%  vn(4,i) = norm(v(1:3,i))
%  vn(5,i) = v(4,i)^2
%
 
%% Recognized options in varargin 
% (uppercases for readability are optional): 
    
% 'Direct': work in direct space (default)
% 'Recip' : work in reciprocal space (
% 'PlotEdges': 3-D plot of edges of Wigner-Seitz cell
% 'PlotFaces': 3-D plot of faces of Wigner-Seitz cell
% 'PlotBasisVectors': plot basis vectors  
% 'PlotLatticeNodes': plot lattice nodes  
% 'PlotCellVertices': plot cell vertices 
% 'pdf': pdf output file
% 'latex': combined pdf & LaTeX output files
% 'Verbose': display some computation details
% if varargin{k} = 'LatticeType', then varargin{k+1 } = string
%                =  Description of lattice type used for:
%                   [-] title of plot window;
%                   [-] switching fast procedure (3D only).
% if varargin{k} = 'Colors', then varargin{k+1 } = 7x3 matrix
%                =  Color order of successive graphical objects
% if varargin{k} = 'ColorMap', 
%                   then varargin{k+1} = Character string 
%                   defining one of the pre-defined colormaps
% if varargin{k} = 'ShowLegend', 
%                   then varargin{k+1} = Boolean to show legend or not
% List of Octave pre-defined colormaps:
% https://octave.sourceforge.io/octave/function/colormap.html

    
%% DEFAUT VALUES OF OPTIONAL ARGUMENTS

default_color=true; default_color_map=true; show_legend=true;
lattice_type= ' ';  
plot_edges = false; plot_faces = false;
plot_basis_vectors = false;
plot_lattice_nodes = false;
plot_cell_vertices = false;
pdf_output=false; latex_output=false;
verbose = false;
space='direct'; Title='Wigner-Seitz cell in direct space';
bvec{1}='a_1';bvec{2}='a_2';bvec{3}='a_3';
coord{1}='x_1';coord{2}='x_2';coord{3}='x_3';

%% PARSE OPTIONAL ARGUMENT LIST

name_value_pair=false;
for k = 1:length(varargin);
   
    if (name_value_pair)
        name_value_pair=false;
    else
        switch lower(varargin{k}) 
          case {'direct'}
            space='direct'; Title='Wigner-Seitz cell in direct space';
            bvec{1} ='a_1';bvec{2} ='a_2'; bvec{3}='a_3';
            coord{1}='x_1';coord{2}='x_2';coord{3}='x_3';
          case {'Reciprocal','reciprocal','Recip','recip'}
            space='recip'; 
            Title='Wigner-Seitz cell in recip space: Brillouin zone';
            bvec{1} ='g_1';bvec{2} ='g_2'; bvec{3}='g_3';
            coord{1}='q_1';coord{2}='q_2';coord{3}='q_3';
          case {'plotedges'}
            plot_edges = true; 
          case {'plotfaces'}
            plot_faces = true;   
          case {'plotbasisvectors'}
            plot_basis_vectors = true; 
          case {'plotcellvertices'}
            plot_cell_vertices = true; 
          case {'plotlatticenodes'}
            plot_lattice_nodes = true; 
          case {'latticetype'}
            lattice_type = varargin{k+1};
            name_value_pair=true;
            Title=[Title ' of ' lattice_type ' lattice'];
          case {'colors'}
            color_scheme = varargin{k+1}; default_color=false; 
            name_value_pair=true;
          case {'colormap'}
            color_map=varargin{k+1}; default_color_map=false;
            name_value_pair=true;
          case {'showlegend'}
            show_legend=varargin{k+1}; name_value_pair=true;
          case {'pdf'}
            pdf_output=true;
          case {'latex'}
            latex_output=true;
          case {'verbose'}
            verbose=true;
          otherwise
            error(['function WignerSeitz: ',...
                   'option %s not recognized.\n'],...
                  varargin{k});
        end
    end
end

%% PRECISION

tol = 1e-12; % Two floating-point numbers will be considered equal 
             % if the absolute value of their difference is < tol

%% PARSE INPUT ARRAY TO DETERMINE DIMENSION

for i=1:3
    norma(i)=norm(a(:,i));
end
dimension=sum(norma>tol);

if (verbose)
    fprintf('\n*************************************************\n')
    fprintf(' WIGNER-SEITZ CELL OF ORDER %d IN %dD %s SPACE\n',1,dimension,upper(space))
    fprintf('*************************************************\n')
end

%% 1-DIMENSION

if (dimension==1) 
    
    [R,neighbour]=lattice(a,2);
    
    k=0; 
    for j=1:3
        if (norma(j) > tol)     % 1-DIMENSION ALONG a(:,j)
            i=j; 
        else
            k=k+1; vn(1:3,k)=a(1:3,j);   % Vectors normal to edge
            vn(5,k)=v(1,k)^2+v(2,k)^2+v(3,k)^2; v(4,1)=sqrt(v(5,k));
        end
    end 
    vertex(:,1)= -a(:,i)/2;
    vertex(:,2)=  a(:,i)/2;
    edge(1,1)=1; edge(2,1)=2; 
    fprintf('Function WignerSeitz in %s space: 2 vertices, 1 edge.\n',...
            space);
    face=[1 2 1];  % Face reduced to an edge    
end

%% 2-DIMENSION

if (dimension==2) 
    if (norm(a(:,3)) < tol )     % 2-DIMENSION IN PLANE x_1-x_2
        ii=1; jj=2; iz=3;
    end
    if (norm(a(:,1)) < tol )     % 2-DIMENSION IN PLANE x_2-x_3
        ii=2; jj=3; iz=1;
    end
    if (norm(a(:,2)) < tol )     % 2-DIMENSION IN PLANE x_1-x_3
        ii=1; jj=3; iz=2;
    end
    
    %% VECTORS NORMAL TO FACE
    
    TU=eye(3,3);
    vn(1:3,1)=TU(1:3,iz); vn(1:3,2)=-TU(1:3,iz); 
    for k=1:2
        vn(5,k)=vn(1,k)^2+vn(2,k)^2+vn(3,k)^2; vnn(4,1)=sqrt(vn(5,k));
    end
    
    %% WORK SET OF 2D LATTICE VECTORS
    
    [T,R,~]=lattice(a,1,1); 
    
    %% Projection to plane x_i-x_j
    
    for n=1:columns(T)
        TT(1,n)=T(ii,n); 
        TT(2,n)=T(jj,n);
    end
    
    for n=1:columns(R)
        RR(1,n)=R(ii,n); 
        RR(2,n)=R(jj,n);
    end
    
    %% Optimization of the work set RR
    %  
    %  If a vector TT(:,m) is not colinear 
    %  with any of the vectors of the current work set RR, 
    %  this means that the direction TT(:,m) is unexplored by the
    %  work set => update work set RR by including TT(:,m).
    %
    %  If a vector TT(:,m) is colinear with some vector RR(:,n)
    %  of the current work set RR and pointing in the same
    %  direction as RR(:,n), there is no need :
    %  [-] to update the work set RR;
    %  [-] to consider the TT(:,j) that are beyond (as viewed from
    %      the origin) the line normal to RR(:,n) and containing RR(:,n).
    
    nodes=columns(RR); newmax=columns(RR); exclus=[];
    for m=nodes+1:columns(TT)    % Search beyond first neighbours
        if (~ismember(m,exclus)) % Check if excluded by a previous step
            colinear=false;
            for n=1:newmax
                sp=TT(1,m)*RR(1,n)+TT(2,m)*RR(2,n);
                cosinus=sp/( sqrt(TT(1,m)^2+TT(2,m)^2) * sqrt(RR(1,n)^2+RR(2,n)^2) );
                if (abs(cosinus-1)<tol)
                    colinear=true; % TT(:,m) and RR(:n) are colinear and
                    break;         % pointing in the same direction
                end
            end
 
            if (colinear)
                for j=m+1:columns(TT)
                    % prp = position of TT(:,j) relatively to plane normal
                    % to RR(:,n) and including RR(:,n)
                    prp=R(1,n)*T(1,j)+R(2,n)*T(2,j)-RR(1,n)^2-RR(2,n)^2;
                    if (prp>0) % TT(:,j) beyond plane (as viewed from the origin)
                        exclus = [exclus j]; % Excluding TT(;,j)
                    end
                end
            else
                newmax=columns(RR)+1; % If TT(:,m) not in the same direction
                RR(:,newmax)=TT(:,m); % as any vector of the current work
                                      % set, include TT(:,m) in work set
            end
        end
    end
   
    fprintf(['Function WignerSeitz using %d nodes after searching ' ...
             'colinear vectors from the origin.\n'],columns(RR))
        
    %% VERTICES OF THE 2D WIGNER-SEITZ CELL
    
    c=RR/2; % The lines of the Wigner-Seitz unit cell pass through
                 % a subset of the points c=zone*RR/2
    
    for i=1:columns(RR) 
        d(i) = RR(1,i)*c(1,i) + RR(2,i)*c(2,i);
    end
    
    % The equation of a line orthogonal to RR(1:2,i) and containing
    % c(1:2,i) = RR(1:2,i)/2 reads:
    %
    % RR(1,i)*x(1) + RR(2,i)*x(2) = d(i) 
    % where
    % d(i) = RR(1,i)*c(1,i) + RR(2,i)*c(2,i)
    
    % The intersection between 2 lines i,j 
    % = a unit cell vertex = solution of the 2X2 linear system of equations:
    %
    % RR(1,i)*x(1) + RR(2,i)*x(2) = d(i) 
    % RR(1,j)*x(1) + RR(2,j)*x(2) = d(j) 
     
    
    %% Search possible vertices
    
    n=0;
    for i=2:columns(RR) % Consider all possible intersections of two lines
        for j=i+1:columns(RR)  
            b = [d(i), d(j)]';
            S(1,1:2) = RR(1:2,i);
            S(2,1:2) = RR(1:2,j);
            if( abs(det(S)) > tol ) 
                x = S\b; % Solve 2X2 linear system of equations
                n=n+1; u(1:2,n)=x(1:2);
            end
        end
    end
    
    if (verbose) 
        fprintf(['Function WignerSeitz: %d intersections as possible ' ...
                 'vertices.\n'],columns(u));
    end
    
    %% Search intersections closest to the origin 
    %% than to any other lattice vector

    k=0; 
    for i=1:n
        dorigin = u(1,i)^2 + u(2,i)^2 ; flag=true;
        for j=2:columns(RR)
            dr = (u(1,i)-RR(1,j))^2 + (u(2,i)-RR(2,j))^2;
            if( dr < dorigin-tol)
                flag=false; % Discard intersection closer to any lattice
                            % vector than to the origin
            end
        end
        if ( flag )
            k=k+1; v(1:2,k)=u(1:2,i);  
        end
    end
    
    if (verbose)  
        fprintf(['Function WignerSeitz: %d intersections closest ' ...
                 'to origin than to any other lattice vector.\n'],columns(v));
    end
    
    %% Discard eventual vertex duplicates
    
    k=0;
    for i=1:columns(v)
        check=true;
        for j=i+1:columns(v)
            if ( abs(v(1,i)-v(1,j)) < tol && abs(v(2,i)-v(2,j)) < tol)
                check=false; % Discard duplicate
            end
        end
        if (check) 
            k=k+1; vertex(1:2,k)=v(1:2,i); vertex(3,k)=0;
        end
    end
    
    if (verbose) 
        fprintf(['Function WignerSeitz: %d vertices after discarding ' ...
                 'duplicates.\n'], columns(vertex));
    end
    
    %% BACK TO 3-DIMENSION VECTORS
      
    buf=zeros(3,columns(vertex));
    for n=1:columns(vertex)
        buf(ii,n)=vertex(1,n); buf(jj,n)=vertex(2,n);
    end
    vertex = buf;
    
    for n=1:columns(RR)
        R(ii,n)=RR(1,n); R(jj,n)=RR(2,n); R(iz,n)=0;
    end
   
    %% EDGES OF 2D WIGNER-SEITZ CELL
     
    ptxy=vertex; ptxy(iz,:)=[];
    ppt=Jarvis(ptxy); 
    % Function Jarvis returns the convex hull of
    % a set of points in a plane, thereby providing the correct 
    % rotation order of the vertices defining a face 
    % that will be suitable for plotting by the function
    % "patch" (see Matlab/Octave documentation).
    
    for n=1:columns(vertex)
        vertex(ii,n)=ppt(1,n); vertex(jj,n)=ppt(2,n);
    end
    
    for n=1:columns(vertex) % Number of edges of a polygon
                            % = number of its vertices
        edge(1,n)=n; edge(2,n)=n+1;
        if (n==columns(vertex))
            edge(2,n)=1;
        end
        
    end
    
    %% Vectors normal to edges
    
    for n=1:columns(vertex)-1
        vne(1:3,n)=vertex(1:3,n)+vertex(1:3,n+1);
    end
    vne(1:3,columns(vertex))=vertex(1:3,columns(vertex))+vertex(1:3,1);
    
    %% SINGLE FACE OF 2D WIGNER-SEITZ CELL
    
    face=[ 1:1:columns(vertex) ];
    
    fprintf('Function WignerSeitz in 2-D %s space: %d vertices, %d edges.\n',...
            space,columns(vertex), columns(vertex));

end % End of 2-dimension
    
%% 3-DIMENSION

if (dimension==3)
    
    %% WORK SET OF 3D LATTICE VECTORS
    
    switch lattice_type
      case {'cP','cI','cF'}
        [~,R,~]=lattice(a,1,2); % Fast lane for high symmetry lattices
        T=R;
      otherwise
        [T,R,~]=lattice(a,3,2); % R(:,n) is the starting work set that
                                % will expanded if needed using 
                                % the "provision" of vectors in array T
    end
    
    %% Optimization of the work set R
    %  
    %  If a vector T(:,m) is not colinear 
    %  with any of the vectors of the current work set R, 
    %  this means that the direction T(:,m) is unexplored by the
    %  work set => update work set R by including T(:,m).
    %
    %  If a vector T(:,m) is colinear with some vector R(:,n)
    %  of the current work set R and pointing in the same
    %  direction as R(:,n), there is no need :
    %  [-] to update the work set R;
    %  [-] to consider the T(:,j) that are beyond (as viewed from
    %      the origin) the plane normal to R(:,n) and containing R(:,n).
    
    nodes=columns(R); newmax=columns(R); exclus=[];
    for m=nodes+1:columns(T)     % Search beyond first set of nodes
        if (~ismember(m,exclus)) % Check if excluded by a previous step
            colinear=false;
            for n=1:newmax
                sp=T(1,m)*R(1,n)+T(2,m)*R(2,n)+T(3,m)*R(3,n);
                cosinus=sp/(T(4,m)*R(4,n));
                if (abs(cosinus-1)<tol)
                    colinear=true; % T(:,m) and R(:n) are colinear and
                    break;         % pointing in the same direction
                end
            end
 
            if (colinear)
                for j=m+1:columns(T)
                    % prp = position of T(:,j) relatively to plane normal
                    % to R(:,n) and including R(:,n)
                    prp=R(1,n)*T(1,j)+R(2,n)*T(2,j)+R(3,n)*T(3,j)-R(5,n);
                    if (prp>0) % T(:,j) beyond plane (as viewed from the origin)
                        exclus = [exclus j]; % Excluding T(;,j)
                    end
                end
            else
                newmax=columns(R)+1;% If T(:,m) not in the same direction
                R(:,newmax)=T(:,m); % as any vector of the current work
                                    % set, include T(:,m) in work set
            end
        end
    end
    if(verbose)
        fprintf(['Function WignerSeitz using %d nodes after searching ' ...
                 'colinear vectors from the origin.\n'],columns(R))
    end
    
    %% VERTICES OF 3D WIGNER-SEITZ CELL
    
    c=R/2; % The planes of the Wigner-Seitz unit cell pass through
           % a subset of the points c=R/2
    
    for i=1:columns(R) 
        d(i) = R(1,i)*c(1,i) + R(2,i)*c(2,i)+ R(3,i)*c(3,i);
    end
    
    % The equation of a plane orthogonal to R(1:3,i) and containing
    % c(1:3,i) = R(1:3,i)/2 reads:
    %
    % R(1,i)*x(1) + R(2,i)*x(2) + R(3,i)*x(3) = d(i) 
    % where
    % d(i) = R(1,i)*c(1,i) + R(2,i)*c(2,i)+ R(3,i)*c(3,i)
    
    % The intersection between 3 planes i,j,k 
    % = a unit cell vertex = solution of the 3X3 linear system of equations:
    %
    % R(1,i)*x(1) + R(2,i)*x(2) + R(3,i)*x(3) = d(i) 
    % R(1,j)*x(1) + R(2,j)*x(2) + R(3,j)*x(3) = d(j) 
    % R(1,k)*x(1) + R(2,k)*x(2) + R(3,k)*x(3) = d(k) 
    
    %% Search possible vertices
    
    nodes=columns(R); nmax=(nodes-1)*(nodes-2)*(nodes-3);
    u=zeros(3,nmax);    
    
    n=0; 
    if(verbose)
        fprintf('Function WignerSeitz searching intersections ...');
    end
    for i=2:nodes
        for j=i+1:nodes 
            for k=j+1:nodes
                b = [d(i), d(j), d(k)]';
                S(1,1:3) = R(1:3,i);
                S(2,1:3) = R(1:3,j);
                S(3,1:3) = R(1:3,k);
                if( abs(det(S)) > tol ) 
                    x = S\b; % Solve 3X3 linear system of equations
                    n=n+1; u(1:3,n)=x(1:3);
                end
            end
        end
    end
    u(:,n+1:nmax)=[];
    if (verbose) 
        fprintf([' %d intersections as possible ' ...
                 'vertices.\n'],n);
    end
    
    %% Search intersections closest to the origin 
    %% than to any other lattice vector
    
    k=0; v=zeros(3,columns(u));
    if(verbose)
        fprintf(['Function WignerSeitz searching intersections closest ' ...
                 'to origin...\n']);  
    end
    
    for i=1:n
        flag=true; dorigin = u(1,i)^2 + u(2,i)^2 + u(3,i)^2;
        for j=2:columns(R)
            dr = (u(1,i)-R(1,j))^2 + (u(2,i)-R(2,j))^2 + (u(3,i)-R(3,j))^2;
            if (dr < dorigin-tol)
                flag=false; % Discard intersection closer to
                break;      % any lattice vector than to the origin
            end
        end
        if (flag)
            k=k+1; v(1:3,k)=u(1:3,i);  
        end
    end
    v(:,k+1:n)=[];
    if (verbose)  
        fprintf(['Function WignerSeitz: %d intersections closest ' ...
                 'to origin than to any other lattice vector.\n'],columns(v));
    end
    
    %% Discard eventual vertex duplicates
    
    k=0;
    for i=1:columns(v)
        check=true;
        for j=i+1:columns(v)
            if ( abs(v(1,i)-v(1,j)) < tol && abs(v(2,i)-v(2,j)) < tol && ...
                 abs(v(3,i)-v(3,j)) < tol)
                check=false; % Discard duplicate
            end
        end
        if (check) 
            k=k+1; vertex(1:3,k)=v(1:3,i); 
        end
    end
    if (verbose) 
        fprintf(['Function WignerSeitz: %d vertices after discarding ' ...
                 'duplicates.\n'], columns(vertex));
    end
    
    %% EDGES OF WIGNER-SEITZ CELL
    
    % The edges of the unit cell are found by considering all
    % pairs of vertices. If a pair of vertices share two planes of the
    % Wigner-Seitz cell boundaries, there is an edge between these two vertices.
    
    n=0;
    for i=1:columns(vertex)
        for j=i+1:columns(vertex)
            for k=2:columns(R) 
                % w(1) = 0 if vertex i in plane k
                w(1) = R(1,k)*vertex(1,i)+R(2,k)*vertex(2,i)+...
                       R(3,k)*vertex(3,i)-d(k);
                % w(2) = 0 if vertex j in plane k
                w(2) = R(1,k)*vertex(1,j)+R(2,k)*vertex(2,j)+...
                       R(3,k)*vertex(3,j)-d(k); 
                if (abs(w(1))<tol && abs(w(2))<tol) 
                    for l=k+1:columns(R) 
                        % w(3) = 0 if vertex i in plane l
                        w(3) = R(1,l)*vertex(1,i)+R(2,l)*vertex(2,i)+...
                               R(3,l)*vertex(3,i)-d(l); 
                        % w(4) = 0 if vertex j in plane l
                        w(4) = R(1,l)*vertex(1,j)+R(2,l)*vertex(2,j)+...
                               R(3,l)*vertex(3,j)-d(l);
                        if(abs(w(3))<tol && abs(w(4))<tol)% Edge connecting
                            n=n+1;edg(1,n)=i;edg(2,n)=j;  % vertices i & j     
                        end
                    end
                end
            end
        end   
    end
    
    %% Discard eventual edge duplicates
    
    k=0;
    for i=1:columns(edg)
        check=true;
        for j=i+1:columns(edg)
            if ( abs(edg(1,i)-edg(1,j)) < tol && abs(edg(2,i)-edg(2,j)) < tol)
                check=false; % Discard duplicate
            end
        end
        if (check) 
            k=k+1; edge(1:2,k)=edg(1:2,i); 
        end
    end
    
    if (verbose)
        fprintf('Function WignerSeitz: %d edges have been found.\n',...       
                 columns(edge))
    end
        
    %% FACES OF WIGNER-SEITZ CELL
    
    % The number of faces is known by Euler's characteristic 
    % formula for convex polyhedrons
    
    Euler=2+columns(edge)-columns(vertex);
    
    %% The faces of the unit cell are found
    %% if a plane contains more than two vertices.

    nmax=0; kf=0; face=zeros(columns(R),columns(vertex));
    for k=2:columns(R) % Avoid starting with R=[0,0,0]
        for i=1:columns(vertex)
            n=0; buf=[];
            for j=1:columns(vertex)
                % w(1) = 0 if vertex i in plane k
                w(1) = R(1,k)*vertex(1,i)+R(2,k)*vertex(2,i)+...
                       R(3,k)*vertex(3,i)-d(k);
                % w(2) = 0 if vertex j in plane k
                w(2) = R(1,k)*vertex(1,j)+R(2,k)*vertex(2,j)+...
                       R(3,k)*vertex(3,j)-d(k);            
                if ( abs(w(1))<tol && abs(w(2))<tol )
                    n=n+1; buf(n)=j;  
                end
            end
            if(n>2) % More than 2 vertices in plane k
                    % Add to face list if not already found
                if (kf==0 || (kf>0 && ~ismember(k,vnind))) 
                    kf=kf+1;
                    face(kf,1:n)=buf(1:n);
                    vnind(kf)=k; % Store index of vector normal to face
                    if (n>nmax)
                        nmax=n;
                    end
                end
            end
        end
    end
      
    face(:,nmax+1:end)=[]; face(kf+1:end,:)=[]; % Minimize array shape 
    
    if(rows(face) ~= Euler)
        fprintf(['\nWARNING: inconsistency in function WignerSeitz: number of ' ...
                 'faces = %d not equal to Euler characteristics = %d\n'],...
                rows(face),Euler)
    else
        fprintf(['Function WignerSeitz in 3D %s space: %d vertices, ' ...
                 '%d edges, %d faces.\n'],...
                 space,columns(vertex),columns(edge),Euler);
    end

    for i=1:rows(face)
        kept(i)=columns(face);
        for j=columns(face):-1:1
            if(face(i,j)==0)
                face(i,j) = NaN; kept(i)=j-1; 
                % See Matlab/Octave manual: NaN allows to plot faces with
                % different numbers of vertices
            end
        end
    end    
 
    %% ORDERING VERTICES OF EACH FACE AS A CONVEX POLYGON FOR PLOTTING
    
    for i=1:rows(face) % Number of faces = number of vectors normal
                       % to faces
        vn(1:5,i)=R(1:5,vnind(i));    % Lattice vector normal to face
        theta= acos(vn(3,i)/vn(4,i)); % Spherical coordinate angles
        phi  = atan2(vn(2,i),vn(1,i));% of vector normal to face i
                                        
        ROTZ=RotationAxe(3,-phi);   % Rotation around axe 3
        ROTY=RotationAxe(2,-theta); % Rotation around axe 2
        
        ROTATION=ROTY*ROTZ; % Rotation bringing vn(:i) parallel to axe 3
       
        pt=[]; ptxy=[]; ppt=[];
        for j=1:kept(i)     % Rotation  applied to all vertices
            pt(1:3,j)=ROTATION*[vertex(1:3,face(i,j))];
        end
        % The resulting pt(:,j) have all the same third component
        % (the rotation moved them in a plane parallel to the plane
        % of Cartesian basis vectors 1 and 2).
        % This third component must be stored to allow a correct 
        % inverse rotation later. Therefore pt is duplicated 
        % as ptxy of which we keep only the components 
        % in the plane of axes 1 and 2 as arguments of the 
        % Jarvis' march algorithm.
        ptxy=pt; ptxy(3,:)=[];
       
        ppt=Jarvis(ptxy); 
        % Function Jarvis returns the convex hull of
        % a set of points in a plane, thereby providing the correct 
        % rotation order of the vertices defining a face 
        % that will be suitable for plotting by the function
        % "patch" (see Matlab/Octave documentation).
        
        for j=1:kept(i)
            pt(1,j)=ppt(1,j); pt(2,j)=ppt(2,j); 
            % pt is now in the correct order to define a convex
            % polygon for the current face. The third component
            % was not changed because it is identical for
            % all points after the rotation of axes.
        end
        
        ROTZ=RotationAxe(3,phi);   % Inverse rotation around axe 3
        ROTY=RotationAxe(2,theta); % Inverse rotation around axe 2
        
        ROTATION=ROTZ*ROTY; % Inverse rotation bringing axe 3 back in the
                            % direction of vn(:,i) 
        
        for j=1:kept(i) % Loop over each vertex of face i
            pto=ROTATION*[pt(1:3,j)]; % Re-ordered vertex after inverse rotation
            for k=1:columns(vertex)   % Search match of pto in array vertex(:,k)
                vk=[vertex(1:3,k)];   % Buffer
                if (norm(pto-vk)< tol)% Matching...
                    face(i,j)=k;      % ... modify vertex j of face i
                end
            end
        end
    end  
    if(verbose)
        fprintf(['Function WignerSeitz ordered vertices for plotting ' ...
                 'faces by function patch.\n']) 
        fprintf('\n***** End of function WignerSeitz ***************\n')
    end
end % End of 3-dimension
 

%% OPTIONAL PLOTS

if (plot_edges || plot_faces || plot_basis_vectors || plot_lattice_nodes ...
    || plot_cell_vertices || pdf_output || latex_output)
    
    origin=[0 0 0]';
    
    h = findall(0,'type','figure');
    if isempty(h)
        figure('NumberTitle', 'off','name',Title);
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
    
    if (default_color) 
        color_scheme = get(gca,'colororder');
    end
    
    if (plot_basis_vectors) 
        for i=1:3
            k=k+1; Legend{k}=bvec{i};
            quiver3(origin(1),origin(2),origin(3),a(1,i),a(2,i),a(3,i),...
                    'autoscale','off','color',color_scheme(i,1:3),...
                    'maxheadsize',0.1,...
                    'DisplayName',Legend{k}); hold on;
        end
    end
    
    if (plot_lattice_nodes) 
        k=k+1; Legend{k}='Lattice nodes used in computations';
        plot3(R(1,:),R(2,:),R(3,:),'o',...
              'Color',color_scheme(3,:),'DisplayName',Legend{k}); hold on;
        if (dimension == 2)
            k=k+1; Legend{k}='Lattice nodes normal to edges ';
            plot3(vne(1,:),vne(2,:),vne(3,:),'o',...
                  'MarkerFaceColor',color_scheme(4,:),...
                  'Color',color_scheme(3,:),'DisplayName',Legend{k}); hold on;
        end
        if (dimension == 3)
            k=k+1; Legend{k}='Lattice nodes normal to faces'; 
            plot3(vn(1,:),vn(2,:),vn(3,:),'o',...
                  'Color',color_scheme(2,:),'DisplayName',Legend{k}); hold on;
        end
    end
    
    if (plot_cell_vertices) 
        k=k+1; Legend{k}='Vertices of Wigner-Seitz cell';
        plot3(vertex(1,:),vertex(2,:),vertex(3,:),'o',...
              'color',color_scheme(2,:),'DisplayName',Legend{k}); hold on;
    end
  
    if (plot_edges)  
        if (~plot_faces)
            k=k+1;Legend{k}='Edges of Wigner-Seitz cell'; 
            for i=1:columns(edge)
                pts = [vertex(:,edge(1,i))'; vertex(:,edge(2,i))'];
                plot3(pts(:,1), pts(:,2), pts(:,3),...
                      '-','color',color_scheme(7,:),'DisplayName',Legend{k}); hold on; 
            end
        end
    end
    
    if (plot_faces) 
        k=k+1; Legend{k}='Wigner-Seitz cell'; 
        if (default_color_map) 
            cmap=colormap('default');
        else
            cmap=colormap(color_map);
        end
        
        switch dimension % Color and transparency of faces
          case 1          
            ty=0.2;
            color_face(1,1:3)=color_scheme(4,:);
            
          case 2
            ty=0.2;
            color_face(1,1:3)=color_scheme(4,:);
            
          case 3
            ty=0.5;      
            offset=rows(cmap)-10; % Avoid extreme colors when changing
                                  % face. Nice when using the
                                  % following corlormaps: viridis,
                                  % gray,bone,pink,autumn,ocean,winter
            for i=1:rows(face)
                color_face(i,1:3)=cmap(offset-i,1:3);
            end
        end
        
        patch ('Vertices', vertex', 'Faces', face, ...
               'LineStyle','-',...
               'LineWidth', 1,...
               'EdgeColor',color_scheme(7,:),... % Options: 'none', 'flat', 'interp', [R G B]
               'EdgeAlpha',1,...     % Edge transparency in range 0 (transparent) to 1 (opaque)
               'FaceColor','flat',...% Options: 'none', 'flat', 'interp', [R G B]
               'FaceAlpha',ty,...    % Face transparency in range 0 (transparent) to 1 (opaque)
               'FaceVertexCData',color_face,...
               'DisplayName',Legend{k}
        ); hold on;
    end
    
    %% Axes, Labels & Legend & View
    
    axis equal; grid on; legend(Legend,'Location','EastOutside');
    xlabel(coord{1}); ylabel(coord{2}); zlabel(coord{3});     
    
    if(~show_legend)
        legend(gca,'off');
    end
    
    % Point of view
    
    view([1,0.5,0.3]);

end

end % End of function WignerSeitz
    