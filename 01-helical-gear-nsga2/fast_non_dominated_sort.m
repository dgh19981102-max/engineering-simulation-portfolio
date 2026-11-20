%% Fast non-dominated sort (constrained domination)
% Deb's constrained-domination rule: a smaller total violation always dominates;
% only when violations are equal (normally both 0, i.e. both feasible) are the
% objectives compared. Feasible designs therefore always outrank infeasible ones
% without any penalty weight to tune. O(M*N^2) for N individuals.
% rank(p) = front index of p; fronts{k} = indices in front k.
function [rank, fronts] = fast_non_dominated_sort(obj, constrViolation)
    n    = size(obj, 1);
    S    = cell(n, 1);    % S{p}: individuals dominated by p
    n_p  = zeros(n, 1);   % n_p(p): number of individuals dominating p
    rank = zeros(n, 1);

    fronts    = cell(1);
    fronts{1} = [];

    for p = 1:n
        S{p}   = [];
        n_p(p) = 0;

        for q = 1:n
            if p ~= q
                if (constrViolation(p) < constrViolation(q)) || ...
                   (constrViolation(p) == constrViolation(q) && ...
                    ((obj(p,1) < obj(q,1) && obj(p,2) <= obj(q,2)) || ...
                     (obj(p,1) <= obj(q,1) && obj(p,2) < obj(q,2))))
                    % p dominates q
                    S{p} = [S{p}, q];
                elseif (constrViolation(q) < constrViolation(p)) || ...
                       (constrViolation(p) == constrViolation(q) && ...
                        ((obj(q,1) < obj(p,1) && obj(q,2) <= obj(p,2)) || ...
                         (obj(q,1) <= obj(p,1) && obj(q,2) < obj(p,2))))
                    % q dominates p
                    n_p(p) = n_p(p) + 1;
                end
            end
        end

        if n_p(p) == 0
            rank(p)   = 1;
            fronts{1} = [fronts{1}, p];
        end
    end

    % Peel fronts: removing front i may leave members of S with no dominator left.
    i = 1;
    while ~isempty(fronts{i})
        Q = [];
        for p = fronts{i}
            for q = S{p}
                n_p(q) = n_p(q) - 1;
                if n_p(q) == 0
                    rank(q) = i + 1;
                    Q = [Q, q];
                end
            end
        end
        i = i + 1;
        fronts{i} = Q;
    end

    % The loop always ends on an empty front; drop it.
    fronts = fronts(1:end-1);
end
