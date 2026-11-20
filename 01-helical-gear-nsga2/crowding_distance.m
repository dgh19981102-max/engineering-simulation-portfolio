%% Crowding distance
% Sum over objectives of the normalised gap between each point's two neighbours
% on the same front. Large values mark sparse regions; preferring them keeps the
% front evenly spread. Boundary points get Inf so the extremes are never lost.
function crowdDist = crowding_distance(obj, rank)
    n = size(obj, 1);
    crowdDist = zeros(n, 1);

    uniqueRanks = unique(rank);

    for i = 1:length(uniqueRanks)
        currentFront = find(rank == uniqueRanks(i));

        if length(currentFront) > 2
            for m = 1:size(obj, 2)
                [sortedObj, sortIdx] = sort(obj(currentFront, m));

                crowdDist(currentFront(sortIdx(1)))   = Inf;
                crowdDist(currentFront(sortIdx(end))) = Inf;

                % Normalise by the objective's range on this front so MPa-scale
                % f1 and smaller f2 contribute comparably. A flat objective
                % (zero range) is skipped to avoid 0/0.
                for j = 2:length(currentFront)-1
                    if sortedObj(end) - sortedObj(1) ~= 0
                        crowdDist(currentFront(sortIdx(j))) = crowdDist(currentFront(sortIdx(j))) + ...
                            (sortedObj(j+1) - sortedObj(j-1)) / (sortedObj(end) - sortedObj(1));
                    end
                end
            end
        else
            % With one or two members every point is a boundary point.
            crowdDist(currentFront) = Inf;
        end
    end
end
