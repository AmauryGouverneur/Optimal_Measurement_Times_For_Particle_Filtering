function [meas_SA,cost_SA,avgCostHist,minCostHist] = SA_algo(n_measurements,T,pop_size,max_gen,n_part,n_draw,measurements_spacing,y,meas_1_j,meas_j_T)
%Simulated Annealing algorithm
if nargin < 7
    measurements_spacing = 1;
end

if nargin < 8 
    online = false ;
    meas_1_j = 0;
    y = 0;
else
    online = true ;
end

len=n_measurements;        % The length of the genomes  
popSize=pop_size;          % The size of the population (must be an even number)
maxGens=max_gen;                % The maximum number of generations allowed in a run
probMutation=ceil(0.25*n_measurements)/n_measurements;        % The mutation probability (per bit)
probChanges = 1;
visualizationFlag=0;       % 0 => don't visualize bit frequencies
                           % 1 => visualize bit frequencies

verboseFlag=0;             % 1 => display details of each generation
                           % 0 => run quietly
convergenceFlag=0;         % 1 => plot convergence curve
                           % 0 => does not


temperature = 10;
% preallocate vectors for recording the average and maximum fitness in each
% generation
avgFitnessHist=zeros(1,maxGens+1);
maxFitnessHist=zeros(1,maxGens+1);
rateAcceptanceHist=zeros(1,maxGens+1);
rateAcceptanceHist_Lowering=zeros(1,maxGens+1);

% the population is a popSize by len matrix of randomly generated boolean
% values

pop = zeros(popSize,len);

if online 
    accessible_meas = measurements_spacing:measurements_spacing:T;
    pop = ones(popSize,1)*meas_j_T;
    
    masks=rand(popSize, len)<probChanges;    
    
    if n_measurements == 1
            pop = sort((1-masks).*pop + masks.*(accessible_meas(unidrnd(length(accessible_meas),popSize,len)))',2);
    else
        % masks(i,j)==1 iff pop(i,j) has to be mutated (0 elsewhere) 
        pop = sort((1-masks).*pop + masks.*(accessible_meas(unidrnd(length(accessible_meas),popSize,len))),2);
        % Replace duplicates measurements and sort
        %parfor i=1:popSize
        for i=1:popSize
            pop(i,:) = replace_duplicates(pop(i,:),accessible_meas);
        end
    end
else 
    accessible_meas = 0:measurements_spacing:T;
    for i=1:popSize
        pop(i,:) = sort(accessible_meas(randperm(length(accessible_meas),len)));
    end
end




% To identify copies in population
pop = sortrows(pop);
gen = 0 ; 

while gen<maxGens  
    % evaluate the fitness of the population. The vector of fitness values 
    % returned  must be of dimensions 1 x popSize.
    counter = 0; 
    counter_Lowering =0;
    fitnessVals=localFitnessFunction(pop);
    
    if gen>0 && maxFitnessHist(1,gen) > max(fitnessVals)
        maxFitnessHist(1,gen+1)= maxFitnessHist(1,gen);
    else 
        [maxFitnessHist(1,gen+1),maxIndex]=max(fitnessVals);
        meas_SA = sort(pop(maxIndex,:));
    end
    avgFitnessHist(1,gen+1)=mean(fitnessVals,'omitnan');
     
    % display the generation number, the average Fitness of the population,
    % and the maximum fitness of any individual in the population
    % Conditionally perform bit-frequency visualization
    
    if visualizationFlag
        figure(1)
        set (gcf, 'color', 'w');
        hold off
        if online
            histogram(pop,accessible_meas,'Normalization','countdensity'); hold on;
            plot(pop(maxIndex,:)+0.5,0*pop(maxIndex,:)+popSize,'.','Markersize',25);
            axis([measurements_spacing T 0 popSize]);
        else 
            histogram(pop,accessible_meas,'Normalization','countdensity'); hold on;
            plot(pop(maxIndex,:)+0.5,0*pop(maxIndex,:)+popSize,'.','Markersize',25);
            axis([0 T 0 popSize]);
        end
        title(['Generation = ' num2str(gen) ', Average Fitness = ' sprintf('%0.3f', avgFitnessHist(1,gen+1))]);
        ylabel('Frequency of measure in t');
        xlabel('time t');
        drawnow;
    end
    
    
    masks=rand(popSize, len)< probMutation*ones(popSize,len);
    % masks(i,j)==1 iff pop(i,j) has to be mutated (0 elsewhere) 
    if n_measurements == 1
        pop_mutation = sort((1-masks).*pop + masks.*(accessible_meas(unidrnd(length(accessible_meas),popSize,len)))',2);
    else
        pop_mutation = sort((1-masks).*pop + masks.*(accessible_meas(unidrnd(length(accessible_meas),popSize,len))),2);
        % Replace duplicates measurements and sort
        %parfor i=1:popSize
        for i=1:popSize
            pop_mutation(i,:) = replace_duplicates(pop_mutation(i,:),accessible_meas);
        end
    end
    
    pop_mutation = sort(pop_mutation,2);
    pop_mutation = sortrows(pop_mutation);
    
    fitnessVals_mutation =localFitnessFunction(pop_mutation);
    
    for i=1:popSize
        if fitnessVals_mutation(i)>fitnessVals(i)
            pop(i,:) = pop_mutation(i,:);
            counter = counter+1;
        else
            if rand < exp(-(fitnessVals(i)-fitnessVals_mutation(i))/temperature)
                if (pop(i,:)~=pop_mutation(i,:))
                    pop(i,:) = pop_mutation(i,:);
                end
                counter = counter+1;
                counter_Lowering = counter_Lowering+1;
            end
        end
    end
    
    
    pop = sort(pop,2);
    pop = sortrows(pop); % to identify copies in population
    
    temperature = temperature*0.8;
    rateAcceptanceHist(gen+1)=counter/pop_size;
    rateAcceptanceHist_Lowering(gen+1)=counter_Lowering/pop_size;
    if verboseFlag
        display(['gen=' num2str(gen,'%.3d') '   avgFitness=' ...
            num2str(avgFitnessHist(1,gen+1),'%3.3f') '   maxFitness=' ...
            num2str(maxFitnessHist(1,gen+1),'%3.3f') '   rateAcc=' ...
            num2str(rateAcceptanceHist(1,gen+1),'%3.3f') '   rateAcc=' num2str(rateAcceptanceHist_Lowering(1,gen+1),'%3.3f') '   temperature=' num2str(temperature,'%3.3f')]);
    end
    gen = gen+1; 
end

display(['fraction of changes= ' num2str(sum(rateAcceptanceHist)/maxGens,'%3.3f')]);
display(['fraction of lowering changes= ' num2str(sum(rateAcceptanceHist_Lowering)/maxGens,'%3.3f')]);

avgCostHist = -avgFitnessHist;
minCostHist = -maxFitnessHist;

cost_SA = -maxFitnessHist(end);


% plot and print
if convergenceFlag
    figure
    set(gcf,'Color','w');
    hold off
    plot(0:maxGens,-avgFitnessHist,'k-');
    hold on
    plot(0:maxGens,-maxFitnessHist,'c-');
    title('Minimum and Average Cost');
    xlabel('Generation');
    ylabel('Cost');
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
