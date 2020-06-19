function [meas_GB,cost_GB,avgCostHist,minCostHist] = greedy_backward_algo(n_measurements,T,n_part,n_draw,measurements_spacing,y,meas_1_j)
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
convergenceFlag=0;         % 1 => plot convergence curve
                           % 0 => does not


% pre-generate two ?repositories? of random binary digits from which the  
% the masks used in mutation and uniform crossover will be picked. 
% maskReposFactor determines the size of these repositories.




% To identify copies in population

if online 
    elite = measurements_spacing:measurements_spacing:T;
else 
    elite = 0:measurements_spacing:T;
end

% preallocate vectors for recording the average and maximum fitness in each
% generation

avGBitnessHist=zeros(1,length(elite)-n_measurements);
maxFitnessHist=zeros(1,length(elite)-n_measurements);

gen = 1;

while length(elite)>n_measurements  
    pop = [ones(length(elite)).*elite];
    pop = reshape(pop(eye(length(elite))~=1),length(elite),length(elite)-1);
    
    % evaluate the fitness of the population. The vector of fitness values 
    % returned  must be of dimensions 1 x popSize.
    popSize = size(pop,1);
    fitnessVals=localFitnessFunction(pop);
    [maxFitnessHist(1,gen),maxIndex]=max(fitnessVals);
    avGBitnessHist(1,gen)=mean(fitnessVals,'omitnan');

    elite = pop(maxIndex,:);
     
    % display the generation number, the average Fitness of the population,
    % and the maximum fitness of any individual in the population
    % Conditionally perform bit-frequency visualization
    if verboseFlag
        display(['gen=' num2str(gen,'%.3d') '   avGBitness=' ...
            num2str(avGBitnessHist(1,gen),'%3.3f') '   maxFitness=' ...
            num2str(maxFitnessHist(1,gen),'%3.3f') ]);
    end
    if visualizationFlag
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
        title(['Generation = ' num2str(gen) ', Average Fitness = ' sprintf('%0.3f', avGBitnessHist(1,gen))]);
        ylabel('Frequency of measure in t');
        xlabel('time t');
        drawnow;
    end    
    gen = gen+1;
end

avgCostHist = -avGBitnessHist;
minCostHist = -maxFitnessHist;

meas_GB = elite;
cost_GB = MC_MSE_estimator(elite,T,500,200);

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
