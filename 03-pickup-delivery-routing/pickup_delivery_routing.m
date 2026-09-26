% PICKUP_DELIVERY_ROUTING  Pickup-and-delivery tour with part-type constraints,
% solved by four methods on the same instance.
%
% Purpose
%   A single carrier (e.g. a pick-and-place arm) starts at a depot, visits
%   n_pairs pickup points and n_pairs place points, carries one part at a time,
%   and returns to the depot. Each part has a type; a part may only be placed on
%   a place point of the same type. Tour length is minimised. This is a TSP with
%   a transition constraint, solved four ways:
%     1. Hybrid : ACO for up to 30% of the budget, then tabu search from its best
%     2. ACO    : plain Ant System
%     3. TS     : tabu search with 2-opt and swap neighbourhoods
%     4. Greedy : nearest feasible neighbour (baseline)
%
% Usage
%   >> pickup_delivery_routing
%   rng(6) at the top fixes the instance and every stochastic run.
%
% Inputs
%   None. Problem size and algorithm settings are in "Parameters"; node layout
%   comes from generate_array_coordinates at the bottom of the file.
%
% Outputs
%   Figure 1 : convergence of the best tour length for all four methods.
%   Figure 2 : hybrid and greedy routes side by side (drawn by
%              plot_route_comparison, which opens its own window).
%   Nothing is written to disk.
%
% Dependencies
%   pdist2 (Statistics and Machine Learning Toolbox); see the NOTE in section 4
%   for a toolbox-free replacement. Everything else is base MATLAB.
%
% Known limitations
%   - The type constraint lives in type_matrix: depot -> pickup,
%     pickup -> place of the same type, place -> any pickup. Pairing is by type,
%     not by a specific pickup-place pair.
%   - The closing edge back to the depot is always added to the length and is
%     never checked against type_matrix.
%   - Construction heuristics fall back to "any unvisited node" when no
%     feasible node is left. That fallback is never reached here, because every
%     type has as many place points as pickups, but it would produce infeasible
%     tours on other instances.
%   - run_Greedy is deterministic: all max_iter repetitions build the same tour,
%     so its curve is flat.
%   - Tabu entries are stored per route position pair (i, j), not per node pair.
%   - The subplot reserved for the route in Figure 1 stays empty, because
%     plot_route_comparison creates a new figure.

clear;
clc;
rng(6);   % fixed seed: same instance and same stochastic runs every time

%% 1. Parameters
n_total  = 37;    % nodes including the depot (-); must equal 2*n_pairs + 1
n_pairs  = 18;    % pickup-place pairs (-)
n_types  = 4;     % part types (-)
n_ants   = 50;    % ants per ACO iteration (-); about the node count is a common choice
max_iter = 100;   % iteration budget shared by all four methods (-)
alpha    = 1.5;   % pheromone exponent (-); 1-2 typical, higher = more exploitation of past tours
beta     = 2;     % heuristic (1/distance) exponent (-); 2-5 typical, higher = greedier ants
rho      = 0.15;  % pheromone evaporation rate per iteration (-); 0.1-0.5 typical.
                  %   run_Hybrid replaces it with adaptive_rho
Q        = 1;     % pheromone deposit constant (-); each ant lays Q/L on its edges, so only
                  %   Q relative to the initial level tau0 = 0.1 and the tour length L matters
prize    = 2;     % elitist bonus (-): prize/L extra pheromone on a new best tour (hybrid only)
punish   = 1;     % penalty coefficient (-); passed to run_Hybrid but not used

%% 2. Part types
% Draw how many pairs belong to each type, summing to n_pairs. Each of the
% first n_types-1 types gets 1..ceil(remaining/types_left); the last type takes
% the remainder, so it is usually the largest.
type_counts = zeros(1, n_types);
remaining = n_pairs;
for i = 1:n_types-1
    if remaining > 0
        type_counts(i) = randi([1, ceil(remaining/(n_types-i+1))]);
        remaining = remaining - type_counts(i);
    end
end
type_counts(n_types) = remaining;

%% 3. Node indices, types and coordinates
% Node numbering: 1 = depot, 2..n_pairs+1 = pickups, n_pairs+2..n_total = places.
% Pickups are grouped by type, and place point k has the same type as pickup k.
coordinates = zeros(n_total, 2);
coordinates(1,:) = [0, 0];   % depot (mm)

% NOTE the random coordinates drawn here are overwritten by
% generate_array_coordinates below. The loop still matters: it builds
% pickup_types, and its rand calls advance the RNG stream, so removing it would
% change every later random draw and every result.
pickup_points = 2:(n_pairs+1);
pickup_types = [];
current_idx = 1;
for t = 1:n_types
    for c = 1:type_counts(t)
        coordinates(current_idx + 1,:) = rand(1,2) * 50;
        pickup_types(current_idx) = t;
        current_idx = current_idx + 1;
    end
end

place_points = (n_pairs+2):n_total;
place_types = pickup_types;   % place point k takes the type of pickup k
% Earlier random placement of the place points, superseded by
% generate_array_coordinates:
% current_idx = n_pairs + 1;
% for t = 1:n_types
%     for c = 1:type_counts(t)
%         coordinates(current_idx + 1,:) = rand(1,2) * 50 + 100;
%         current_idx = current_idx + 1;
%     end
% end
coordinates = generate_array_coordinates(n_total, n_pairs, n_types, type_counts);

%% 4. Distance matrix
% Euclidean distance between every pair of nodes (mm).
% NOTE pdist2 belongs to the Statistics and Machine Learning Toolbox. A
% toolbox-free equivalent (R2016b+ implicit expansion) that gives the same values:
%     dist_matrix = sqrt(sum((permute(coordinates, [1 3 2]) - ...
%                             permute(coordinates, [3 1 2])).^2, 3));
% The Gram-matrix form sqrt(max(s + s.' - 2*G, 0)), with G = coordinates*coordinates.'
% and s = diag(G), is faster for large N but does not guarantee exact zeros on
% the diagonal (round-off in s + s.' - 2*G). A tiny non-zero value there makes
% eta finite, eta(isinf(eta)) = 0 no longer clears it, and self-loops get a huge
% heuristic weight, so add dist_matrix(1:N+1:end) = 0 if you use it.
dist_matrix = pdist2(coordinates, coordinates);

%% 5. Type-constraint matrix
% type_matrix(a, b) = 1 if the carrier may travel directly from node a to b.
% The carrier holds one part at a time, so the tour must alternate
% pickup -> matching-type place -> pickup -> ...
type_matrix = zeros(n_total);
% Depot -> any pickup (the carrier starts empty).
type_matrix(1, pickup_points) = 1;
% Pickup -> place points of the same type only.
for i = 1:length(pickup_points)
    p_idx = pickup_points(i);
    for j = 1:length(place_points)
        pl_idx = place_points(j);
        if pickup_types(i) == place_types(j)
            type_matrix(p_idx, pl_idx) = 1;
        end
    end
end
% Place -> any pickup (the carrier is empty again).
type_matrix(place_points, pickup_points) = 1;

%% 6. Run the four methods on the same instance
[hybrid_best_route, hybrid_min_lengths, hybrid_avg_lengths] = run_Hybrid(dist_matrix, type_matrix, pickup_types, place_types, n_ants, max_iter, alpha, beta, rho, Q, prize, punish);

[aco_best_route, aco_min_lengths, aco_avg_lengths] = run_ACO(dist_matrix, type_matrix, pickup_types, place_types, n_ants, max_iter, alpha, beta, rho, Q);

[ts_best_route, ts_min_lengths, ts_avg_lengths] = run_TS(dist_matrix, type_matrix, pickup_types, place_types, max_iter);

[greedy_best_route, greedy_min_lengths] = run_Greedy(dist_matrix, type_matrix, pickup_types, place_types, max_iter);

%% 7. Plots
figure('Position', [100 100 1200 600]);

% Route comparison. plot_route_comparison opens its own figure, so this subplot
% area stays empty.
subplot('Position', [0.05 0.1 0.55 0.85]);
plot_route_comparison(coordinates, hybrid_best_route, greedy_best_route, pickup_types);

% Best-so-far tour length (mm) against iteration for all four methods.
subplot('Position', [0.65 0.1 0.3 0.85]);
iterations = 1:length(hybrid_min_lengths);
plot(iterations, hybrid_min_lengths, 'b-', ...
     iterations, aco_min_lengths, 'r--', ...
     iterations, ts_min_lengths, 'k:', ...
     iterations, greedy_min_lengths, 'g-.');
legend('混合优化算法', '蚁群算法', '禁忌搜索算法', '贪心算法');
title('最短路径长度对比');
xlabel('迭代次数');
ylabel('最短路径长度');
grid on;

%% Greedy baseline
% Nearest feasible neighbour from the depot. Repeated max_iter times only so
% its curve has the same length as the others.
% NOTE nothing in the loop is random, so every repetition builds the same tour.
function [greedy_best_route, greedy_min_lengths] = run_Greedy(dist_matrix, type_matrix, pickup_types, place_types, max_iter)
    n_cities    = size(dist_matrix, 1);
    best_length = inf;
    best_route  = [];
    min_lengths = zeros(max_iter, 1);

    for iter = 1:max_iter
        route = ones(1, n_cities);
        visited = zeros(1, n_cities);
        visited(1) = 1;
        current_pos = 1;

        for i = 2:n_cities
            current = route(current_pos);
            allowed = find(~visited & type_matrix(current,:));

            % Fallback: take any unvisited node (see file header).
            if isempty(allowed)
                allowed = find(~visited);
            end

            [~, min_idx] = min(dist_matrix(current, allowed));
            next_city = allowed(min_idx);

            route(i) = next_city;
            visited(next_city) = 1;
            current_pos = i;
        end

        current_length = calculate_path_length(route, dist_matrix);

        if current_length < best_length
            best_length = current_length;
            best_route  = route;
        end

        min_lengths(iter) = best_length;
    end

    greedy_best_route  = best_route;
    greedy_min_lengths = min_lengths;
end

%% Route comparison figure (hybrid vs greedy)
function plot_route_comparison(coordinates, hybrid_route, greedy_route, pickup_types)
    figure('Position', [100 100 1500 600]);

    subplot(1, 2, 1);
    plot_single_route(coordinates, hybrid_route, pickup_types, '混合优化算法路径图');

    subplot(1, 2, 2);
    plot_single_route(coordinates, greedy_route, pickup_types, '贪心算法路径图');
end

%% Draw one route
% Depot = red diamond; pickups = filled circles; places = hollow circles; colour
% = part type. Arrows at edge midpoints show the direction of travel.
function plot_single_route(coordinates, route, pickup_types, title_text)
    hold on;

    % One colour per type (MATLAB default palette entries 2, 4, 5, 6).
    colors = {[0.8500 0.3250 0.0980],  % orange-red
              [0.4940 0.1840 0.5560],  % purple
              [0.4660 0.6740 0.1880],  % green
              [0.3010 0.7450 0.9330]}; % sky blue

    scatter(coordinates(1,1), coordinates(1,2), 250, [0.8 0 0], 'filled', 'diamond', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    text(coordinates(1,1), coordinates(1,2), '1', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'FontSize', 12);

    % Place point k is node k + 1 + n_pairs, the partner slot of pickup k.
    n_pairs = length(pickup_types);
    for i = 1:n_pairs
        pickup_idx = i + 1;
        type_idx   = pickup_types(i);
        color_idx  = mod(type_idx-1, length(colors)) + 1;   % wrap if there are more types than colours

        scatter(coordinates(pickup_idx,1), coordinates(pickup_idx,2), 150, ...
            colors{color_idx}, 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        text(coordinates(pickup_idx,1), coordinates(pickup_idx,2), ...
            num2str(pickup_idx), 'HorizontalAlignment', 'right', 'FontWeight', 'bold');

        place_idx = pickup_idx + n_pairs;
        scatter(coordinates(place_idx,1), coordinates(place_idx,2), 120, ...
            colors{color_idx}, 'o', 'LineWidth', 2.5);
        text(coordinates(place_idx,1), coordinates(place_idx,2), ...
            num2str(place_idx), 'HorizontalAlignment', 'right');
    end

    % Edges; arrow length is 1/10 of the edge, placed at its midpoint.
    for i = 1:length(route)-1
        p1 = coordinates(route(i),:);
        p2 = coordinates(route(i+1),:);
        plot([p1(1),p2(1)], [p1(2),p2(2)], 'k-', 'LineWidth', 1);
        arrow_pos = (p1 + p2) / 2;
        quiver(arrow_pos(1), arrow_pos(2), ...
            (p2(1)-p1(1))/10, (p2(2)-p1(2))/10, 0, ...
            'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);
    end

    % Closing edge back to the depot (the tour is a cycle).
    p1 = coordinates(route(end),:);
    p2 = coordinates(route(1),:);
    plot([p1(1),p2(1)], [p1(2),p2(2)], 'k-', 'LineWidth', 1);
    arrow_pos = (p1 + p2) / 2;
    quiver(arrow_pos(1), arrow_pos(2), ...
        (p2(1)-p1(1))/10, (p2(2)-p1(2))/10, 0, ...
        'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);

    % Legend built from invisible NaN markers, one per type, so each type
    % appears once however many points it has.
    h = zeros(max(pickup_types) + 1, 1);
    legend_entries = cell(max(pickup_types) + 1, 1);

    h(1) = scatter(NaN, NaN, 250, [0.8 0 0], 'filled', 'diamond', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    legend_entries{1} = '起始点';

    for i = 1:max(pickup_types)
        color_idx = mod(i-1, length(colors)) + 1;
        h(i+1) = scatter(NaN, NaN, 150, colors{color_idx}, 'filled', 'o', 'MarkerEdgeColor', 'k');
        legend_entries{i+1} = sprintf('型号%d (实心:抓取点 空心:放置点)', i);
    end

    legend(h, legend_entries, 'Location', 'best');
    title(title_text);
    xlabel('X坐标(mm)');
    ylabel('Y坐标(mm)');
    grid on;
    axis equal;
    hold off;
end

%% Hybrid: ACO followed by tabu search
% ACO explores globally and quickly finds a good basin; tabu search then
% refines that tour locally. ACO gets at most 30% of max_iter and stops early
% after max_stagnant iterations without improvement; TS uses the rest, so both
% phases together use exactly max_iter iterations, like the other methods.
% Two additions over plain ACO: an elitist bonus on each new best tour, and a
% 2-opt local search every 10 iterations.
function [best_route, min_lengths, avg_lengths] = run_Hybrid(dist_matrix, type_matrix, pickup_types, place_types, n_ants, max_iter, alpha, beta, rho, Q, prize, punish)
    n_cities     = size(dist_matrix, 1);
    aco_max_iter = floor(max_iter * 0.3);   % ACO phase cap (-)

    %% Phase 1: ant colony
    tau = ones(n_cities) * 0.1;   % initial pheromone tau0 (-); uniform, so the first ants follow distance only
    eta = 1 ./ dist_matrix;       % heuristic desirability 1/d (1/mm): shorter edges more attractive
    eta(isinf(eta)) = 0;          % diagonal d = 0 gives Inf; set to 0 so no self-loops

    best_route  = [];
    min_length  = inf;
    min_lengths = zeros(max_iter, 1);
    avg_lengths = zeros(max_iter, 1);

    % Early stop: 20 iterations without a new best means the colony has
    % converged and further ACO iterations are better spent in TS.
    stagnant_count   = 0;
    max_stagnant     = 20;
    last_best_length = inf;   % assigned but not read

    actual_aco_iter = 0;
    for iter = 1:aco_max_iter
        actual_aco_iter = iter;
        routes  = zeros(n_ants, n_cities);
        lengths = zeros(n_ants, 1);

        for ant = 1:n_ants
            route = build_route_with_types(tau, eta, alpha, beta, dist_matrix, type_matrix, pickup_types, place_types);
            routes(ant,:) = route;
            lengths(ant)  = calculate_path_length(route, dist_matrix);
        end

        [iter_min_length, best_ant] = min(lengths);
        if iter_min_length < min_length
            min_length = iter_min_length;
            best_route = routes(best_ant,:);
            stagnant_count = 0;

            % Elitist bonus prize/L on the new best tour's edges pulls the colony
            % towards it; shorter tours get a larger bonus.
            delta_tau_prize = prize / min_length;
            for i = 1:n_cities-1
                tau(best_route(i), best_route(i+1)) = tau(best_route(i), best_route(i+1)) + delta_tau_prize;
            end
        else
            stagnant_count = stagnant_count + 1;
        end

        % 2-opt polish of the best tour every 10 iterations. Doing it every
        % iteration would cost O(n^2) route evaluations each time.
        if mod(iter, 10) == 0
            [best_route, min_length] = local_search(best_route, min_length, dist_matrix, type_matrix);
        end

        min_lengths(iter) = min_length;
        avg_lengths(iter) = mean(lengths);

        % Ant System update: evaporate, then every ant deposits Q/L on its
        % edges. The evaporation rate grows over the phase (adaptive_rho),
        % overriding the rho argument.
        rho = adaptive_rho(iter, aco_max_iter);
        delta_tau = zeros(n_cities);
        for ant = 1:n_ants
            route = routes(ant,:);
            for i = 1:n_cities-1
                delta_tau(route(i),route(i+1)) = delta_tau(route(i),route(i+1)) + Q/lengths(ant);
            end
        end
        tau = (1-rho) * tau + delta_tau;

        if stagnant_count >= max_stagnant
            break;
        end

        last_best_length = min_length;
    end

    %% Phase 2: tabu search from the ACO best
    current_route  = best_route;
    current_length = min_length;

    % Tenure 7: a reversed move stays forbidden for 7 iterations; 5-10 is the
    % usual range for problems of this size.
    tabu_tenure         = 7;
    tabu_matrix         = zeros(n_cities);
    ts_stagnation_count = 0;
    ts_max_stagnation   = 20;   % iterations without improvement before diversifying

    for iter = actual_aco_iter+1:max_iter
        [next_route, next_length, move_i, move_j, neighbor_lengths] = find_best_neighbor(current_route, dist_matrix, type_matrix, tabu_matrix, iter, min_length);

        if next_length < min_length
            min_length = next_length;
            best_route = next_route;
            ts_stagnation_count = 0;
        else
            ts_stagnation_count = ts_stagnation_count + 1;
        end

        % Age all tabu entries by one, then forbid the move just made.
        % move_i = 0 means no admissible neighbour was found this iteration.
        tabu_matrix = max(0, tabu_matrix - 1);
        if move_i > 0 && move_j > 0
            tabu_matrix(move_i, move_j) = tabu_tenure;
            tabu_matrix(move_j, move_i) = tabu_tenure;
        end

        % Diversification: after a long stall, perturb the current tour and
        % clear the tabu list so the search can leave the local optimum.
        if ts_stagnation_count >= ts_max_stagnation
            current_route  = diversify_solution(current_route, type_matrix);
            current_length = calculate_path_length(current_route, dist_matrix);
            ts_stagnation_count = 0;
            tabu_matrix = zeros(n_cities);
        else
            % Always move to the best admissible neighbour, even if it is
            % worse; that is what lets TS climb out of local minima.
            current_route  = next_route;
            current_length = next_length;
        end

        min_lengths(iter) = min_length;
        avg_lengths(iter) = mean(neighbor_lengths);
    end
end

%% Plain ant colony (Ant System)
% Same construction and deposit rule as the hybrid's first phase, with a fixed
% evaporation rate, no elitist bonus and no local search. Serves as the
% reference for what those additions buy.
function [best_route, min_lengths, avg_lengths] = run_ACO(dist_matrix, type_matrix, pickup_types, place_types, n_ants, max_iter, alpha, beta, rho, Q)
    n_cities = size(dist_matrix, 1);
    tau = ones(n_cities) * 0.1;   % initial pheromone tau0 (-)
    eta = 1 ./ dist_matrix;       % heuristic desirability 1/d (1/mm)
    eta(isinf(eta)) = 0;

    best_route  = [];
    min_length  = inf;
    min_lengths = zeros(max_iter, 1);
    avg_lengths = zeros(max_iter, 1);

    for iter = 1:max_iter
        routes  = zeros(n_ants, n_cities);
        lengths = zeros(n_ants, 1);

        for ant = 1:n_ants
            route = build_route_with_types(tau, eta, alpha, beta, dist_matrix, type_matrix, pickup_types, place_types);
            routes(ant,:) = route;
            lengths(ant)  = calculate_path_length(route, dist_matrix);
        end

        [iter_min_length, best_ant] = min(lengths);
        if iter_min_length < min_length
            min_length = iter_min_length;
            best_route = routes(best_ant,:);
        end

        min_lengths(iter) = min_length;
        avg_lengths(iter) = mean(lengths);

        % Evaporate, then deposit Q/L per ant: short tours reinforce their
        % edges more, and evaporation forgets old choices.
        delta_tau = zeros(n_cities);
        for ant = 1:n_ants
            route = routes(ant,:);
            for i = 1:n_cities-1
                delta_tau(route(i),route(i+1)) = delta_tau(route(i),route(i+1)) + Q/lengths(ant);
            end
        end
        tau = (1-rho) * tau + delta_tau;
    end
end

%% Plain tabu search
% Starts from a random feasible tour. Each iteration moves to the best
% admissible neighbour (2-opt or swap), even if it is worse.
function [best_route, min_lengths, avg_lengths] = run_TS(dist_matrix, type_matrix, pickup_types, place_types, max_iter)
    n_cities = size(dist_matrix, 1);

    current_route  = generate_initial_solution(type_matrix, pickup_types, place_types);
    current_length = calculate_path_length(current_route, dist_matrix);

    best_route  = current_route;
    best_length = current_length;
    min_lengths = zeros(max_iter, 1);
    avg_lengths = zeros(max_iter, 1);

    % Tabu entries are position pairs (i, j) of the last moves, not whole
    % tours: storing tours would cost memory and almost never repeat exactly.
    tabu_tenure = 7;                % iterations a move stays forbidden (-); 5-10 usual
    tabu_matrix = zeros(n_cities);  % remaining tabu time per position pair

    stagnation_count = 0;
    max_stagnation   = 20;   % iterations without improvement before diversifying

    for iter = 1:max_iter
        [next_route, next_length, move_i, move_j, neighbor_lengths] = find_best_neighbor(current_route, dist_matrix, type_matrix, tabu_matrix, iter,best_length);

        if next_length < best_length
            best_length = next_length;
            best_route  = next_route;
            stagnation_count = 0;
        else
            stagnation_count = stagnation_count + 1;
        end

        % Age all entries, then forbid the move just made.
        % move_i = 0 means no admissible neighbour was found this iteration.
        tabu_matrix = max(0, tabu_matrix - 1);
        if move_i > 0 && move_j > 0
            tabu_matrix(move_i, move_j) = tabu_tenure;
            tabu_matrix(move_j, move_i) = tabu_tenure;
        end

        % Diversify after a stall. Unlike the hybrid, the tabu list is kept.
        if stagnation_count >= max_stagnation
            current_route  = diversify_solution(current_route, type_matrix);
            current_length = calculate_path_length(current_route, dist_matrix);
            stagnation_count = 0;
        else
            current_route  = next_route;
            current_length = next_length;
        end

        min_lengths(iter) = best_length;
        avg_lengths(iter) = mean(neighbor_lengths);
    end
end

%% Best admissible neighbour (2-opt + swap)
% Scans every position pair (i, j), i >= 2 so the depot stays first. Two
% neighbourhoods are used because 2-opt (reverse a segment) repairs crossing
% edges while swap (exchange two nodes) moves nodes between distant parts of the
% tour. Neighbours that break the type constraint are skipped.
% Returns the best admissible tour, the move (i, j) and all feasible neighbour
% lengths (used for the average-length curve).
function [best_route, best_length, move_i, move_j, neighbor_lengths] = find_best_neighbor(current_route, dist_matrix,  type_matrix, tabu_matrix, iter,global_best_length)
    n_cities         = length(current_route);
    best_route       = current_route;
    best_length      = inf;
    move_i           = 0;
    move_j           = 0;
    neighbor_lengths = [];

    for i = 2:n_cities-1
        for j = i+1:n_cities
            new_route_2opt = current_route;
            new_route_2opt(i:j) = current_route(j:-1:i);

            new_route_swap = current_route;
            new_route_swap([i,j]) = new_route_swap([j,i]);

            for route_idx = 1:2
                if route_idx == 1
                    new_route = new_route_2opt;
                else
                    new_route = new_route_swap;
                end

                valid = check_type_constraints(new_route, type_matrix);

                if valid
                    new_length = calculate_path_length(new_route, dist_matrix);
                    neighbor_lengths = [neighbor_lengths; new_length];

                    % Aspiration criterion: a tabu move is still allowed if it
                    % beats the best tour found so far, since it cannot be
                    % cycling back to a visited solution.
                    if new_length < global_best_length || tabu_matrix(i,j) == 0
                        if new_length < best_length
                            best_length = new_length;
                            best_route  = new_route;
                            move_i      = i;
                            move_j      = j;
                        end
                    end
                end
            end
        end
    end

    if isempty(neighbor_lengths)
        neighbor_lengths = best_length;
    end
end

%% Type-constraint check
% True if every consecutive edge of the route is allowed by type_matrix.
% The closing edge back to the depot is not checked.
function valid = check_type_constraints(route, type_matrix)
    valid = true;
    for k = 1:length(route)-1
        if ~type_matrix(route(k), route(k+1))
            valid = false;
            break;
        end
    end
end

%% Diversification
% 2-4 random swaps, each kept only if the tour stays feasible. A small kick:
% large enough to leave the current basin, small enough to keep most of the
% tour's structure.
function new_route = diversify_solution(route, type_matrix)
    n = length(route);
    new_route = route;
    num_swaps = randi([2,4]);

    for i = 1:num_swaps
        % Positions from 2 so the depot stays first.
        idx1 = randi([2,n]);
        idx2 = randi([2,n]);

        temp_route = new_route;
        temp_route([idx1,idx2]) = temp_route([idx2,idx1]);

        if check_type_constraints(temp_route, type_matrix)
            new_route = temp_route;
        end
    end
end

%% Ant tour construction
% From the current node, choose among feasible unvisited nodes with probability
% proportional to tau^alpha * eta^beta (random-proportional rule of Ant System).
% Restricting the choice to type_matrix builds the constraint into the
% encoding, so no penalty term is needed.
function route = build_route_with_types(tau, eta, alpha, beta, dist_matrix, type_matrix, pickup_types, place_types)
    n_cities = size(tau, 1);
    route = ones(1, n_cities);
    visited = zeros(1, n_cities);
    visited(1) = 1;
    current_pos = 1;

    for i = 2:n_cities
        current = route(current_pos);
        allowed = find(~visited & type_matrix(current,:));

        % Fallback: take any unvisited node (see file header).
        if isempty(allowed)
            allowed = find(~visited);
        end

        P = zeros(1, length(allowed));
        for k = 1:length(allowed)
            next = allowed(k);
            P(k) = (tau(current,next)^alpha) * (eta(current,next)^beta);
        end
        P = P / sum(P);

        ncity = allowed(roulettewheel(P));
        route(i) = ncity;
        visited(ncity) = 1;
        current_pos = i;
    end
end

%% Roulette-wheel selection
% Returns index k with probability P(k): the first cumulative sum >= a uniform draw.
function selected = roulettewheel(P)
    r = rand();
    C = cumsum(P);
    selected = find(r <= C, 1);
end

%% Closed tour length (mm)
% Sum of consecutive edges plus the closing edge back to the start.
function total_dist = calculate_path_length(route, dist_matrix)
    total_dist = 0;
    for i = 1:length(route)-1
        total_dist = total_dist + dist_matrix(route(i), route(i+1));
    end
    total_dist = total_dist + dist_matrix(route(end), route(1));
end

%% Adaptive evaporation rate
% rho = 1 - exp(-lambda*iter/max_iter), clipped to [0.1, 0.7]. Low evaporation
% early keeps the colony exploring; high evaporation late forgets poor edges
% fast and concentrates the search. lambda = 2 (-) sets how fast rho ramps up.
function rho = adaptive_rho(iter, max_iter)
    lambda = 2;
    rho = 1 - exp(-lambda * iter / max_iter);
    rho = min(max(rho, 0.1), 0.7);
end

%% 2-opt local search
% Tries every segment reversal of the tour and keeps the single best
% improvement (best-improvement, one pass). Only feasible reversals count.
function [improved_route, improved_length] = local_search(route, current_length, dist_matrix, type_matrix)
    improved_route  = route;
    improved_length = current_length;

    for i = 2:length(route)-1
        for j = i+1:length(route)
            new_route = route;
            new_route(i:j) = route(j:-1:i);

            % Same test as check_type_constraints, written inline.
            valid = true;
            for k = 1:length(new_route)-1
                if ~type_matrix(new_route(k), new_route(k+1))
                    valid = false;
                    break;
                end
            end

            if valid
                new_length = calculate_path_length(new_route, dist_matrix);
                if new_length < improved_length
                    improved_route  = new_route;
                    improved_length = new_length;
                end
            end
        end
    end
end

%% Random feasible initial tour for tabu search
% Like the greedy construction but picks uniformly among feasible nodes, so TS
% starts from an unbiased point.
function initial_route = generate_initial_solution(type_matrix, pickup_types, place_types)
    n_cities = size(type_matrix, 1);
    initial_route = ones(1, n_cities);
    visited = zeros(1, n_cities);
    visited(1) = 1;
    current_pos = 1;

    for i = 2:n_cities
        current = initial_route(current_pos);
        allowed = find(~visited & type_matrix(current,:));

        if isempty(allowed)
            allowed = find(~visited);
        end

        next_city = allowed(randi(length(allowed)));
        initial_route(i) = next_city;
        visited(next_city) = 1;
        current_pos = i;
    end
end

%% Single-route plot (not called)
% Older version of plot_single_route that draws into the current axes; kept for
% reuse, not called by the script.
function plot_route(coordinates, route, pickup_types)
    hold on;
    cla;   % clear the current axes first

    colors = {[0.8500 0.3250 0.0980],  % orange-red
              [0.4940 0.1840 0.5560],  % purple
              [0.4660 0.6740 0.1880],  % green
              [0.3010 0.7450 0.9330]}; % sky blue

    scatter(coordinates(1,1), coordinates(1,2), 250, [0.8 0 0], 'filled', 'diamond', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    text(coordinates(1,1), coordinates(1,2), '1', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'FontSize', 12);

    n_pairs = length(pickup_types);
    for i = 1:n_pairs
        pickup_idx = i + 1;
        type_idx   = pickup_types(i);
        color_idx  = mod(type_idx-1, length(colors)) + 1;

        % Pickup: large filled circle
        scatter(coordinates(pickup_idx,1), coordinates(pickup_idx,2), 150, ...
            colors{color_idx}, 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        text(coordinates(pickup_idx,1), coordinates(pickup_idx,2), ...
            num2str(pickup_idx), 'HorizontalAlignment', 'right', 'FontWeight', 'bold');

        % Place: smaller hollow circle with a thick edge
        place_idx = pickup_idx + n_pairs;
        scatter(coordinates(place_idx,1), coordinates(place_idx,2), 120, ...
            colors{color_idx}, 'o', 'LineWidth', 2.5);
        text(coordinates(place_idx,1), coordinates(place_idx,2), ...
            num2str(place_idx), 'HorizontalAlignment', 'right');
    end

    for i = 1:length(route)-1
        p1 = coordinates(route(i),:);
        p2 = coordinates(route(i+1),:);
        plot([p1(1),p2(1)], [p1(2),p2(2)], 'k-', 'LineWidth', 1);
        arrow_pos = (p1 + p2) / 2;
        quiver(arrow_pos(1), arrow_pos(2), ...
            (p2(1)-p1(1))/10, (p2(2)-p1(2))/10, 0, ...
            'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);
    end

    p1 = coordinates(route(end),:);
    p2 = coordinates(route(1),:);
    plot([p1(1),p2(1)], [p1(2),p2(2)], 'k-', 'LineWidth', 1);
    arrow_pos = (p1 + p2) / 2;
    quiver(arrow_pos(1), arrow_pos(2), ...
        (p2(1)-p1(1))/10, (p2(2)-p1(2))/10, 0, ...
        'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);

    h = zeros(max(pickup_types) + 1, 1);
    legend_entries = cell(max(pickup_types) + 1, 1);

    h(1) = scatter(NaN, NaN, 250, [0.8 0 0], 'filled', 'diamond', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    legend_entries{1} = '起始点';

    for i = 1:max(pickup_types)
        color_idx = mod(i-1, length(colors)) + 1;
        h(i+1) = scatter(NaN, NaN, 150, colors{color_idx}, 'filled', 'o', 'MarkerEdgeColor', 'k');
        legend_entries{i+1} = sprintf('型号%d (实心:抓取点 空心:放置点)', i);
    end

    legend(h, legend_entries, 'Location', 'best');

    title('路径规划结果（不同颜色和标记表示不同型号）');
    xlabel('X坐标(mm)');
    ylabel('Y坐标(mm)');
    grid on;
    axis equal;
    hold off;
end

%% Node layout
% Depot at (0, 0) mm. Pickups sit on a grid, 5 per row with 20 mm pitch,
% ordered by type and shifted up-right by 0.4*spacing. Place points are
% scattered at random in a band below-left of the grid, at least 0.8*spacing
% apart (rejection sampling), then shifted down-left by 0.4*spacing.
% pickup_points, pickup_types and place_points are rebuilt here but not returned.
function coordinates = generate_array_coordinates(n_total, n_pairs, n_types, type_counts)
    coordinates = zeros(n_total, 2);
    spacing  = 20;   % grid pitch (mm)
    row_size = 5;    % pickups per grid row (-)

    coordinates(1,:) = [0, 0];   % depot (mm)

    % Grid slots for the pickups, rows going downwards (-y).
    num_rows  = ceil(n_pairs / row_size);
    base_grid = zeros(n_pairs, 2);
    for i = 1:n_pairs
        row = floor((i-1) / row_size);
        col = mod(i-1, row_size);
        base_grid(i,:) = [col * spacing, -row * spacing];
    end

    % Pickups: grid slot plus an up-right offset of 0.4*spacing (8 mm).
    pickup_points = 2:(n_pairs+1);
    pickup_types  = [];
    current_idx   = 1;

    for t = 1:n_types
        for c = 1:type_counts(t)
            base_pos = base_grid(current_idx,:);
            coordinates(current_idx + 1,:) = base_pos + [spacing*0.4, spacing*0.4];
            pickup_types(current_idx) = t;
            current_idx = current_idx + 1;
        end
    end

    % Place-point band (mm), sized from the grid height so it scales with n_pairs.
    min_x = -spacing * (num_rows + 2);
    max_x = 0;
    min_y = -spacing*(num_rows +4);
    max_y = -spacing*(num_rows-1 );

    place_points   = (n_pairs+2):n_total;
    current_idx    = n_pairs + 1;
    used_positions = [];

    for t = 1:n_types
        for c = 1:type_counts(t)
            % Rejection sampling: redraw until the point is at least min_dist
            % from every earlier place point, so markers do not overlap and no
            % two place points coincide (a zero off-diagonal distance would be
            % zeroed out of eta by eta(isinf(eta)) = 0 and never chosen).
            while true
                x = min_x + rand() * (max_x - min_x);
                y = min_y + rand() * (max_y - min_y);

                min_dist = spacing * 0.8;   % minimum spacing (mm), 16 mm
                if isempty(used_positions) || all(sqrt(sum((used_positions - [x,y]).^2, 2)) > min_dist)
                    break;
                end
            end

            % Down-left offset of 0.4*spacing (8 mm).
            coordinates(current_idx + 1,:) = [x - spacing*0.4, y - spacing*0.4];
            used_positions = [used_positions; x, y];
            current_idx = current_idx + 1;
        end
    end

    return
end
