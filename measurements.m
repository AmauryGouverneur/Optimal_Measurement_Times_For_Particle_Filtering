function y_j = measurements(x_j,T,t_j)
% The function "measurements" simulates a realisation of a measurement 
% vector of x_j according  to a discrete stochastic nonlinear dynamical
% systems. 
% 
% y(t) = g(x(t),v(t)) 
%
% Inputs: 
%   - x_j : state vector
%   - T : length of the time interval
% 
% Output: 
%   - y_j : observation vector 
%
% Date : 10/06/20
% Author : Amaury Gouverneur & Antoine Aspeel 
v_j =  randn(1,T+1) ;
y_j = x_j(1,:).*sin((t_j:T+t_j).*x_j(2,:) ) + x_j(3,:)+ v_j;

end
