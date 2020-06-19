function part = initialization(n_draw,y,meas_1_j)
% Function initializing the particles according to the model 
% X_0 ~ Xsi(x_0) dx_0
% 
% Input : 
%   - n_part : number of particles
% 
% Output : 
%   - part : particles
% 
% Date : 10/06/20
% Author : Amaury Gouverneur & Antoine Aspeel

a_minus = 8.8;
a_plus = 24;
ds = 0.25;
omega_minus = 0.13/ds;
omega_plus = 0.21/ds;

b_minus = -5.8;
b_plus = 5.8;

minus = [a_minus;omega_minus;b_minus];
plus = [a_plus;omega_plus;b_plus];

if nargin == 1
    part = [unifrnd(minus(1),plus(1),[1,n_draw]);unifrnd(minus(2),plus(2),[1,n_draw]);unifrnd(minus(3),plus(3),[1,n_draw])];
else
    if meas_1_j == 0 
        part = [unifrnd(minus(1),plus(1),[1,n_draw]);unifrnd(minus(2),plus(2),[1,n_draw]);unifrnd(minus(3),plus(3),[1,n_draw])];
    else 
        t_j = meas_1_j(end);
        index_t_j = t_j + 1;
        y_j = y(1:index_t_j);
        measurement_times = zeros(1,t_j+1); 
        measurement_times(meas_1_j+1) = 1;
        n_part = max(n_draw,10000);
        part = initialization(n_part);
        [~,part] = particle_filter(y_j,measurement_times,t_j,part,0)  ;
        part = part(:,1:n_draw);
    end
end

end