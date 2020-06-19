%% Plot function 
close all;

n_measurements = 11; %number of observations 
T = 30; %number of time steps 
measurements_spacing = 1; 
n_part_fine = 1000; 


n_part = 100; %number of particles 250
n_draw = 200; %number of draws for the MC MSE estimator 10

pop_size = 30; %population size of the GA algorithm 
max_gen = 15; %maximum number of generations

n_draw_comp = 500; %number of draws to realize 

try
    load(sprintf('online_computations_%d_particles_%d_draws',n_part,n_draw));
    n_draw_done = length(find(mse_reg));
    mse_reg = [mse_reg(1:n_draw_done),zeros(1,n_draw_comp)];
    mse_GA = [mse_GA(1:n_draw_done),zeros(1,n_draw_comp)];
    mse_GA_online = [mse_GA_online(1:n_draw_done),zeros(1,n_draw_comp)];

catch
    mse_reg =  zeros(1,n_draw_comp); 
    mse_GA =  zeros(1,n_draw_comp); 
    mse_GA_online = zeros(1,n_draw_comp); 
    n_draw_done = 0;
end

n_draw_expected = 10000;
mse_reg_apriori = zeros(1,n_draw_expected);
x0 = initialization(n_draw_expected,0,0);
for j = 1:n_draw_expected
        %1.Simulation

        %1.1. Random motion model : X_j
        x_j = model(T,x0(:,j),0);

        %1.2. Artificial data record Y_j
        y_j = measurements(x_j,T,0);

        %2. Filtering
        part = initialization(n_part_fine,0,0);
        tau_j_reg = particle_filter(y_j,measurements_reg,T,part,0);
        mse_reg_apriori(j) = mean((objective(x_j,0)-tau_j_reg).^2); 
end


%% 1. Comparison of the performances of the methods by running n_draw_comp draws
display(['Comparison of the performances of the methods by running draw number ',num2str(n_draw_done+1,'%.0f'),' to draw number ',num2str(n_draw_done + n_draw_comp,'%.0f')]);
t_start_comparison = tic;

meas_reg = round(linspace(0,T,n_measurements));

measurements_reg = zeros(1,T+1); 
measurements_reg(meas_reg+1) = 1;


for j = n_draw_done+1:n_draw_done+n_draw_comp
x0 = initialization(1,0,0);
x = model(T,x0,0); %True state vector we want to reconstruct using online OIMPF
y = measurements(x,T,0);
% Computation of the measurement times reg, GA, GA_online 

disp('GA online (+GA) computations started');
t_start_GA_online = tic;
[meas_GA_online,meas_GA] = online_optimob(y,n_measurements,T,pop_size,max_gen,n_part,n_draw,measurements_spacing);
t_elapsed_GA_online = toc(t_start_GA_online);

measurements_GA = zeros(1,T+1); 
measurements_GA(meas_GA+1) = 1;

measurements_GA_online = zeros(1,T+1); 
measurements_GA_online(meas_GA_online+1) = 1;

part = initialization(n_part_fine,0,0);
tau_j_reg = particle_filter(y,measurements_reg,T,part,0);
tau_j_GA = particle_filter(y,measurements_GA,T,part,0);
tau_j_GA_online = particle_filter(y,measurements_GA_online,T,part,0);

mse_reg(j) = mean((objective(x,0)-tau_j_reg).^2,2);
mse_GA(j) =  mean((objective(x,0)-tau_j_GA).^2,2);
mse_GA_online(j) =  mean((objective(x,0)-tau_j_GA_online ).^2,2);

display(['GA online computations completed, mean gain = ',num2str( (mean(mse_reg)-mean(mse_GA_online))/mean(mse_reg)*100,'%.1f%%'), ', gain apriori = ',num2str((mean(mse_reg)-mean(mse_GA))/mean(mse_reg)*100,'%.1f%%')  ,', time elapsed = ',num2str(t_elapsed_GA_online,'%.0f sec')]);

save(sprintf('online_computations_%d_particles_%d_draws',n_part,n_draw),'mse_reg','mse_GA','mse_GA_online')
end

t_elapsed_comparison = toc(t_start_comparison);

display(['Computations of ',num2str(n_draw_done, '%.0f'),'draws, time elapsed = ',num2str(t_elapsed_comparison,'%.0f sec')]);


%% 2. Results and save

n_draw_done = length(find(mse_reg));
mse_reg = mse_reg(1:n_draw_done);
mse_GA = mse_GA(1:n_draw_done);
mse_GA_online = mse_GA_online(1:n_draw_done);
B = 1000; 
alpha = 0.05;


mse_error_margin_reg = norminv(1-alpha/2)*std(mse_reg)/(sqrt(n_draw_done));
mse_error_margin_reg_apriori = norminv(1-alpha/2)*std(mse_reg_apriori)/(sqrt(n_draw_expected));

gains_GA = (mse_reg-mse_GA)./mse_reg;
gains_GA_online  = (mse_reg-mse_GA_online)./mse_reg; 

bootstrap_error_median_GA = sort(median(gains_GA)-bootstrp(B,@median,gains_GA));
bootstrap_error_pos_frac_gain_GA = sort(mean(gains_GA>=0)-bootstrp(B,@(x) mean(x>=0),gains_GA));

median_error_margin_GA = [bootstrap_error_median_GA(ceil(alpha/2*B)),bootstrap_error_median_GA(ceil((1-alpha/2)*B))];
pos_frac_error_margin_GA = [bootstrap_error_pos_frac_gain_GA(ceil(alpha/2*B)),bootstrap_error_pos_frac_gain_GA(ceil((1-alpha/2)*B))];

bootstrap_error_median_GA_online = sort(median(gains_GA_online)-bootstrp(B,@median,gains_GA_online));
bootstrap_error_pos_frac_gain_GA_online = sort(mean(gains_GA_online>=0)-bootstrp(B,@(x) mean(x>=0),gains_GA_online));

median_error_margin_GA_online = [bootstrap_error_median_GA_online(ceil(alpha/2*B)),bootstrap_error_median_GA_online(ceil((1-alpha/2)*B))];
pos_frac_error_margin_GA_online = [bootstrap_error_pos_frac_gain_GA_online(ceil(alpha/2*B)),bootstrap_error_pos_frac_gain_GA_online(ceil((1-alpha/2)*B))];


mean_gain_error_margin_GA = norminv(1-alpha/2)*std(gains_GA)/(sqrt(n_draw_done));
mse_error_margin_GA = norminv(1-alpha/2)*std(mse_GA)/(sqrt(n_draw_done)*mean(mse_GA));

mean_gain_error_margin_GA_online = norminv(1-alpha/2)*std(gains_GA_online)/(sqrt(n_draw_done));
mse_error_margin_GA_online = norminv(1-alpha/2)*std(mse_GA_online)/(sqrt(n_draw_done)*mean(mse_GA_online));


display(['GA : gain on expected mse = ',num2str((1-(mean(mse_GA))/mean(mse_reg_apriori))*100,'%.1f%%'),' (+/- ',num2str((mse_error_margin_GA)/(1-mse_error_margin_reg_apriori)*100,'%.1f%%)') , ', average gain = ' num2str(mean(gains_GA)*100,'%.1f%%'),' (+/- ',num2str(mean_gain_error_margin_GA*100,'%.1f%%)'), ', median gain = ' num2str(mean((median(gains_GA)+median_error_margin_GA)*100),'%.1f%%'),' (+/- ',num2str(mean(abs(median_error_margin_GA))*100,'%.1f%%)'), ', fraction of positive gain = ' num2str(mean((mean(gains_GA>=0)+pos_frac_error_margin_GA)*100),'%.1f%%'),' (+/- ',num2str(mean(abs(pos_frac_error_margin_GA))*100,'%.1f%%)')]);
display(['GA online : gain on expected mse = ',num2str((1-(mean(mse_GA_online))/mean(mse_reg_apriori))*100,'%.1f%%'),' (+/- ',num2str((mse_error_margin_GA_online)/(1-mse_error_margin_reg_apriori)*100,'%.1f%%)') , ', average gain = ' num2str(mean(gains_GA_online)*100,'%.1f%%'),' (+/- ',num2str(mean_gain_error_margin_GA_online*100,'%.1f%%)'), ', median gain = ' num2str(mean((median(gains_GA_online)+median_error_margin_GA_online)*100),'%.1f%%'),' (+/- ',num2str(mean(abs(median_error_margin_GA_online))*100,'%.1f%%)'), ', fraction of positive gain = ' num2str(mean((mean(gains_GA_online>=0)+pos_frac_error_margin_GA_online)*100),'%.1f%%'),' (+/- ',num2str(mean(abs(pos_frac_error_margin_GA_online))*100,'%.1f%%)')]);

%% 3. Display histogram 
if true 
figure;    
bins = linspace(-2,1,16);
histogram(gains_GA,bins,'Normalization','pdf')
xlim([-2, 1]);
ylim([0,2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GA algorithm')    
    
figure;    
histogram(gains_GA_online,bins,'Normalization','pdf')
xlim([-2, 1]);
ylim([0,2]);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GA online algorithm')
end