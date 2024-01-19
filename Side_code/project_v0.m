clear all;  % Clear all variables
clc;        % Clear the command window
close all;  % Close all figures
source('mystartdefaults.m'); % Contains SI physical sonstants

tic

% Units for free electrons in vacuum
recipunit = 1.0E+10;
ekinscale = ((hbar * recipunit)^2 / (2.0 * elm))/qel

%% INPUT PARAMETERS

datafile='Electron-1D.dat'; % Output file name
nband=4;                    % Number of bands

a=4.08;                     % Lattice spacing [constant]

%% Reciprocal lattice vectors

ng=10; % Maximum number of reciprocal lattice vectors
ngx=floor(ng/2); % Number of reciprocal lattice vectors in x-direction
ng=2*ngx+1 % Make sure that the size of the Hamiltonian is odd

n=0;
for i=-ng:ng    % We will need a number of G vectors
% that is at least twice the size of the Hamiltonian
    n=n+1;
    G(n)=1; % In 1D, reciprocal lattice vectors
            % in 2*pi/a units are integers
end

%% Fourier components of the potential

fprintf(['\nReciprocal lattice vectors and associated Fourier'...
'coefficients of potential energy\n']);

fprintf(['n    G(2*pi/a units)    Re(V_G) [eV]     Im(V_G) [eV] \n']);

V_G=zeros(1,2*ng+1); % Fourier components of the potential
n=0
for i=-ng:ng
    n=n+1;
    if(abs(i)==1)
        V_G(n)=1;
    end
    if(abs(i)==2)
        V_G(n)=0.75;
    end
end

n=0;
for i=-ng:ng
    n=n+1;
    fprintf('%2d    %15.6G    %15.6G    %15.6G \n',n,G(n),real(V_G(n)),imag(V_G(n)));
end

%% Sampling of reciprocal space
% 2*pi/a units

nk=301; % Number of k-points
k_max=1.5
k_min=-k_max
dk=(k_max-k_min)/(nk-1);

%% Trick to draw BZ boundaries
% 2*pi/a units

ho=35 % [eV]

for i=1:nk
    k(i)=k_min+(i-1)*dk;
    zb(i)=0;
    for j=-ngx:ngx
        if(k(i)>(G(ng+1+j)+0.5))
            zb(i)=ho*(-1)^abs(j);
        end
        if(k(i)<(G(ng+1+j)+0.5))
            zb(i)=0;
        end
    end
end

%% Off-diagonal elements of the Hamiltonian are indep of k

H=zeros(ng,ng); % Ini Ham matrix
j=0;
for jx=-ngx:ngx
    j=j+1;
    i=0
    for ix=-ngx:ngx
        i=i+1;
        H(i,j)=V_G(ng+1+jx-ix);  % NB: G(ng+1)=0
    end
end

%% Main Loop

f1=fopen(datafile,'w'); % Open output file

fprintf('\nLOOP OVER %4d WAVEVECTORS\n',nk);
for m=1:nk
    %% Kinetic Energy
    n=0;
    for i=-ngx:ngx
        n=n+1;
        H(n,n)=ekinscale*(2*pi/a)^2*(k(m)-G(ng+1+i))^2 + V_G(ng+1);
    end

    %% Diagonalization
    tol=1.0E-10;
    if(~ishermitian(H,tol))
        error('\nERROR: H is not hermitian\n');
    end
    [v,ev]=eig(H); % Diagonalize H
    E=real(diag(ev)); % Extract eigenvalues
    [E,perm]=sort(E); % Sort eigenvalues
    v=v(:,perm); % Sort eigenvectors

    fprintf(f1,'%15.6G %15.6G',k(m), zb(m)); % Write k and z to file
    for i=1:nband
        fprintf(f1,' %15.6G',E(i)); % Write eigenvalues to file
    end
    fprintf(f1,'\n'); % Write newline to file
end
fclose(f1); % Close output file
fprintf('\nResults: %2d bands in file Electron-1D.dat\n',nband);
type(datafile); % Display output file