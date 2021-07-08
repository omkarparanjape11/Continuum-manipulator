% Gives abscissas (qs) and weights (ws) for gauss quadrature
% Ngp: number of gauss points required
function [qs, ws] = quad_data(Ngp)

switch Ngp
   case 1
      qs = 0;
      ws = 2;
   case 2
      qs = [-0.5774; 0.5774];
      ws = [1; 1];
    case 3
       qs=[-0.7446;0;0.7446];
       ws=[0.555556;0.8888889;0.555556];
    case 4
        qs = [-0.8611;-0.3399;0.3399;0.8611];
        ws = [0.3478;0.6521;0.6521;0.3478];
   otherwise
      error('Quadrature order not defined')
end