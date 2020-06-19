function z_j = objective(x_j,t_j)
% The function "objective" return an objective vector z_j of x_j
% according  to a discrete stochastic nonlinear dynamical systems. 
% 
% z(t) = h(x(t)) for t = 0,...,T
% 
% Input: 
%   - x_j : state vector
% 
% Outputs : 
%   - z_j : objective vector 
%  
% Date : 10/06/20
% Author : Amaury Gouverneur & Antoine Aspeel
T_j = size(x_j,2)-1;
z_j = x_j(1,:).*sin((t_j:T_j+t_j).*x_j(2,:) ) + x_j(3,:);

end