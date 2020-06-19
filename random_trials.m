function [meas_RT,cost_RT,avgCostHist_RT,minCostHist_RT] = random_trials(T,n_measurements, pop_size, n_part, n_draw,measurements_spacing,y,meas_1_j)

if nargin < 6
    measurements_spacing = 1;
end

if nargin < 8 
    online = false ;
    meas_1_j = 0;
    y = 0;
else
    online = true ;
end

len=n_measurements;                     % The length of the genomes  
popSize=pop_size;           % The size of the population (must be an even number)
pop = zeros(popSize,len);

if online 
    accessible_meas = measurements_spacing:measurements_spacing:T;
else 
    accessible_meas = 0:measurements_spacing:T;
end

for i=1:popSize
    pop(i,:) = sort(accessible_meas(randperm(length(accessible_meas),len)));
end


% To identify copies in population
pop = sortrows(pop);
fitnessVals = localFitnessFunction(pop);

all_costs_RT = -fitnessVals;
[cost_RT,ind_opt] = min(all_costs_RT);
meas_RT = pop(ind_opt,:);

if nargout > 2
    all_costs_RT = all_costs_RT(randperm(popSize));
end

n_eval = length(all_costs_RT);

avgCostHist_RT = zeros(1,n_eval);
minCostHist_RT = zeros(1,n_eval);
for i=1:n_eval
    avgCostHist_RT(i) = mean(all_costs_RT(1:i));
    minCostHist_RT(i) = min(all_costs_RT(1:i));
end

    function fitness = localFitnessFunction(pop)
       % function to MAXIMIZE
       % is designed to compute only once the fitness in cases of copies of
       % individuals
        [popLocSize,~] = size(pop);
        
        % first time that an individual appears
        indFirstCopy = find(sum( (pop(1:end-1,:)-pop(2:end,:)).^2 ,2)~=0)'+1;
        indFirstCopy = [1 indFirstCopy];
        
        % measurements of these first individuals
        firstMeas = pop(indFirstCopy,:);
        
        firstFitnesses = zeros(1,length(indFirstCopy)-1);
        %parfor j = 1:length(indFirstCopy)
        for j = 1:length(indFirstCopy)
            meas = firstMeas(j,:);
             if online
                firstFitnesses(j)  = - MC_MSE_estimator(meas,T,n_draw,n_part,y,meas_1_j);
            else 
                firstFitnesses(j)  = - MC_MSE_estimator(meas,T,n_draw,n_part);
             end            
        end
        
        % copy the fitnesses for similar individuals
        fitness = repelem(firstFitnesses,diff([indFirstCopy popLocSize+1]));
    end
end
