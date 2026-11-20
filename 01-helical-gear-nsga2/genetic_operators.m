%% Genetic operators: selection, crossover, mutation
% Binary tournament -> simulated binary crossover (SBX) -> polynomial mutation,
% the standard real-coded NSGA-II operators. z1 is rounded and every gene is
% clipped to its bounds at the end.
function offspring = genetic_operators(pop, rank, crowdDist, pc, pm, etaC, etaM, lb, ub)
    popSize = size(pop, 1);
    numVar  = size(pop, 2);

    % Binary tournament on the crowded-comparison operator: lower rank wins,
    % ties go to the less crowded individual.
    matingPool = zeros(popSize, numVar);
    for i = 1:popSize
        idx1 = randi(popSize);
        idx2 = randi(popSize);

        if rank(idx1) < rank(idx2) || (rank(idx1) == rank(idx2) && crowdDist(idx1) > crowdDist(idx2))
            matingPool(i,:) = pop(idx1,:);
        else
            matingPool(i,:) = pop(idx2,:);
        end
    end

    % SBX (Deb & Agrawal, 1995): the spread factor beta is drawn so that the
    % children's spread mimics single-point crossover on binary strings. etaC
    % controls how far children land from the parents. min(i+1, popSize) pairs
    % the last parent with itself when popSize is odd.
    offspring = zeros(popSize, numVar);
    for i = 1:2:popSize
        if rand() < pc
            parent1 = matingPool(i,:);
            parent2 = matingPool(min(i+1, popSize),:);

            beta = zeros(1, numVar);   % SBX spread factor (-), unrelated to the helix angle
            mu   = rand(1, numVar);

            for j = 1:numVar
                if mu(j) <= 0.5
                    beta(j) = (2 * mu(j))^(1 / (etaC + 1));
                else
                    beta(j) = (1 / (2 * (1 - mu(j))))^(1 / (etaC + 1));
                end

                % The two children are symmetric about the parents' mean.
                offspring(i,j) = 0.5 * ((1 + beta(j)) * parent1(j) + (1 - beta(j)) * parent2(j));
                offspring(min(i+1, popSize),j) = 0.5 * ((1 - beta(j)) * parent1(j) + (1 + beta(j)) * parent2(j));
            end
        else
            offspring(i,:) = matingPool(i,:);
            offspring(min(i+1, popSize),:) = matingPool(min(i+1, popSize),:);
        end
    end

    % Polynomial mutation (basic, unbounded form): perturbation delta in [-1, 1]
    % scaled by the variable's range; large etaM keeps most steps small.
    % Out-of-range results are clipped below rather than reflected.
    for i = 1:popSize
        for j = 1:numVar
            if rand() < pm
                y = offspring(i,j);
                delta = zeros(1, numVar);
                r = rand();

                if r < 0.5
                    delta(j) = (2 * r)^(1 / (etaM + 1)) - 1;
                else
                    delta(j) = 1 - (2 * (1 - r))^(1 / (etaM + 1));
                end

                offspring(i,j) = y + delta(j) * (ub(j) - lb(j));
            end
        end
    end

    % Tooth count must be an integer.
    offspring(:,1) = round(offspring(:,1));

    % Clip to the design box.
    for i = 1:popSize
        for j = 1:numVar
            offspring(i,j) = max(lb(j), min(ub(j), offspring(i,j)));
        end
    end
end
