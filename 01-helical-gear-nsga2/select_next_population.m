%% Environmental selection
% Sort by (rank ascending, crowding distance descending) and keep the first
% popSize. Equivalent to NSGA-II's front-by-front fill with the last front
% truncated by crowding distance.
function nextIdx = select_next_population(rank, crowdDist, popSize)
    [~, sortedIdx] = sortrows([rank, -crowdDist]);
    nextIdx = sortedIdx(1:popSize);
end
