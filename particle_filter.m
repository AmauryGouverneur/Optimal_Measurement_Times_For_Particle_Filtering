function [tau_j,part] = particle_filter(y_j,measurement_times,T,part,t_j) 
% The function "particle_filter" computes an estimate of the state vector
% x_j. This is performed using a SISR particle filter using measurement y_j
% given at a given measurements times. 
% 
% A discrete stochastic nonlinear dynamical system is modelled:
%   x(t+1) = f(x(t),w(t)) for t = 0,...,T ? 1,
%   y(t) = g(x(t),v(t)) for t in a given measurment set
%   z(t) = h(x(t)) for t = 0,...,T 
%   x(0) ~ F 
%   w(t) and v(t) are the process and measurement noise respectively witn
%   known probability density functions
%
% The following functions are implemented according to the discrete 
% stochastic nonlinear dynamical system : 
%   - initialization(n_part) : initializes the particles according to F 
%   - weighting(y_j,part,t) : return p(y_j|part) at time t
%   - mutation(part,t): mutates the particles according to the dynamical 
%   system
% 
% Inputs: 
%   - y_j : observation vector
%   - measurement_times : binary vector of size T+1, 1 indicating a 
%   measurment time, otherwise 0 
%   - n_part : number of particles
%   - T : length of the time interval
% 
% Outputs: 
%   - tau_j : estimate of x using the particle filter
%   - part : particles at the last time step 
%
% Date : 10/06/20
% Author : Amaury Gouverneur & Antoine Aspeel

n_part = size(part,2);
size_objective = size(objective(part(:,1),0),1);
tau_j = zeros(size_objective,T+1);

if measurement_times(0+1)
        w = weighting(y_j,part,0,t_j); % weighting
        tau_j(:,1) = (objective_part(part,t_j)*w)/sum(w); % estimation
        ind = randsample(n_part,n_part,true,w);
        part = part(:,ind);
else 
    tau_j(:,1) = mean(objective_part(part,t_j));
end

for t = 1:T
    %mutation 
    part = mutation(part,t_j+t-1);
    index_t = t+1; 
    if measurement_times(index_t)
        w = weighting(y_j,part,t,t_j);
        tau_j(:,index_t) = (objective_part(part,t_j+t)*w)/sum(w) ;
        %selection
        if norm(w)==0 
            tau_j(:,index_t) = mean(objective_part(part,t_j+t)); 
        else
            ind = randsample(n_part,n_part,true,w);
            part = part(:,ind);
        end
    else
        tau_j(:,index_t) = mean(objective_part(part,t_j+t));
    end 
end

end