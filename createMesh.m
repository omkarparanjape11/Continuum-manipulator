function [xs, eConn,xConn] = createMesh(x0, L, Nnodes)
% Creates a mesh in 1D with 6DOF at each node
% Inputs:
%  x0: starting coordinate of the domain
%  L : length of the domain
%  Nnodes: number of nodes to be used for discretization
% Outputs:
%  xs: Nnodesx1 vector containing nodal coordinates in sequence
%  eConn: (Nnodes-1) x 12 matrix containing element connectivity
%  xConn: Connects eConn to cordinates correpondingly.

xs = linspace(x0, x0+L, Nnodes)';
eConn = zeros(Nnodes-1,12);
xConn = zeros(Nnodes-1,12);

%% For loop to fill in eConn and xConn
n = 1;
for i = 1:Nnodes-1
    eConn(i,:) = n:n+11;
    xConn(i,:) = [xs(i),xs(i),xs(i),xs(i),xs(i),xs(i),xs(i+1),xs(i+1),xs(i+1),xs(i+1),xs(i+1),xs(i+1)];
    n = n+6;
end



