function [meas_GF,cost_GF,avgCostHist,minCostHist] = greedy_forward_algo(n_measurements,T,n_part,n_draw,measurements_spacing,y,meas_1_j)
%%Simulated Annealing algorithm

if nargin < 5
    measurements_spacing = 1;
end
if nargin < 6
    online = false ;
    meas_1_j = 0;
    y = 0;
else
    online = true ;
end

visualizationFlag=0;       % 0 => don't visualize bit frequencies
                           % 1 => visualize bit frequencies

verboseFlag=0;             % 1 => display details of each generation
                           % 0 => run quietly

% preallocate vectors for recording the average and maximum fitness in each
% generation
avgFitnessHist=zeros(1,n_measurements);
maxFitnessHist=zeros(1,n_measurements);


% To identify copies in population
elite = [] ; 

if online 
    initial_candidates = measurements_spacing:measurements_spacing:T;
else 
    initial_candidates = 0:measurements_spacing:T;
end


candidates = initial_candidates;
pop = candidates'; 
while length(elite)<n_measurements  
    
    gen = length(elite)+1;
    fitnessVals=localFitnessFunction(pop);
    [maxFitnessHist(1,gen),maxIndex]=max(fitnessVals);
    avgFitnessHist(1,gen)=mean(fitnessVals,'omitnan');

    elite = pop(maxIndex,:);
    candidates_measurement_times = ones(1,length(initial_candidates)); 
    candidates_measurement_times(ismember(initial_candidates,elite)) = 0;
    candidates = initial_candidates(ones(1,length(initial_candidates)) == candidates_measurement_times);

    if verboseFlag
        display(['gen=' num2str(gen,'%.3d') '   avgFitness=' ...
            num2str(avgFitnessHist(1,gen),'%3.3f') '   maxFitness=' ...
            num2str(maxFitnessHist(1,gen),'%3.3f') ]);
    end
    if visualizationFlag
        popSize = length(candidates);
        figure(1)
        set (gcf, 'color', 'w');
        hold off
        if online
            histogram(pop,1:T,'Normalization','countdensity'); hold on;
            plot(pop(maxIndex,:)+0.5,0*pop(maxIndex,:)+popSize,'.','Markersize',25);
            axis([1 T 0 popSize]);
        else 
            histogram(pop,0:T,'Normalization','countdensity'); hold on;
            plot(pop(maxIndex,:)+0.5,0*pop(maxIndex,:)+popSize,'.','Markersize',25);
            axis([0 T 0 popSize]);
        end
        title(['Generation = ' num2str(gen) ', Average Fitness = ' sprintf('%0.3f', avgFitnessHist(1,gen))]);
        ylabel('Frequency of measure in t');
        xlabel('time t');
        drawnow;
    end    
    
    pop = [ones(length(candidates),length(elite)).*elite,candidates'];
    pop = sort(pop,2);
end

avgCostHist = -avgFitnessHist;
minCostHist = -maxFitnessHist;

meas_GF = elite;
cost_GF = -maxFitnessHist(end);

function fitness = localFitnessFunction(pop)
   % function to MAXIMIZE
   fitness = zeros(size(pop,1),1);
    %parfor j = 1:size(pop,1)
    for j = 1:size(pop,1)
        meas = pop(j,:);
        if online
            fitness(j)  = - MC_MSE_estimator(meas,T,n_draw,n_part,y,meas_1_j);
        else 
            fitness(j)  = - MC_MSE_estimator(meas,T,n_draw,n_part);
        end
    end
end
end
