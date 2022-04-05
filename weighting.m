function w = weighting(y_j,part,t,t_j)
% The function "weighting" return the probability of observing y_j(t) given
% part is the "true" state vector at time t
%
% y(t) | x_(t) ~ p(y(t) | x(t)) dy(t)  
%
% Input: 
%   - y_j : measurement vector
%   - part : particles 
%   - t : time 
% 
% Outputs : 
%   - w : weigth vector
% 
% Date : 10/06/20
% Author : Amaury Gouverneur & Antoine Aspeel
v = 1;
index_t = t+1;
w_l = log_normpdf(y_j(index_t),objective_part(part,t_j+t),v )';
w = exp(w_l - max(w_l)) ;
end
