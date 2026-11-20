% HELICAL_GEAR_NSGA2  Multi-objective geometry optimisation of a helical gear pair.
%
% Purpose
%   For a fixed centre distance, face width and gear ratio, search the geometry
%   (z1, mn, beta, xn1, xn2) that (f1) minimises flank contact stress and
%   (f2) minimises the gap between pinion and wheel root bending stress, so that
%   both gears are equally utilised in bending. The two goals conflict, so the
%   result is a Pareto front. The problem set-up follows the paper
%   "Multi-objective optimisation of helical gear geometry based on MOPSO";
%   NSGA-II (Deb et al., 2002) is used here in place of MOPSO.
%
% Usage
%   Run the script from its own folder:  >> helical_gear_nsga2
%   The helper functions live in separate .m files in the same folder
%   (evaluate_population, tooth_factors, calculate_transverse_contact_ratio,
%   fast_non_dominated_sort, crowding_distance, genetic_operators,
%   select_next_population, solve_inv). The script sets rng(42) before
%   initialising the population, so runs are repeatable; change the seed in
%   the "Reproducibility" section below for a different run.
%
% Inputs
%   None. All operating data and algorithm settings are hard-coded in the
%   "Problem parameters" and "NSGA-II settings" sections below.
%
% Outputs
%   Console : every feasible non-dominated design, a feasibility check of the
%             paper's baseline design, and the balanced trade-off design.
%   Figure  : Pareto front with representatives A / B / C and the baseline.
%   File    : pareto_front.png (200 dpi) in the current folder.
%
% Dependencies
%   Base MATLAB only. exportgraphics needs R2020a or later.
%
% Known limitations
%   - Load factors K1 = K2 = 1: no application, dynamic or load-distribution
%     factors (KA, KV, KHbeta, ...). Stresses are nominal, not ISO 6336 ratings.
%   - Form factor YF, root thickness sFn and bending arm hFe (tooth_factors) come
%     from simplified empirical fits, not the ISO 6336-3 method B construction.
%   - Transverse contact ratio and tip thickness use spur-gear relations with the
%     normal module and normal pressure angle (see the NOTE in each function).
%   - Helix factor Ybeta = 1 - beta/120deg is the ISO form with the overlap ratio
%     capped at 1; it is only exact when eps_beta >= 1.
%   - The wheel bending stress sigmaF2 divides the pinion torque by z2 instead of
%     z1 (see the NOTE at objective 2), which understates it by a factor of u.
%   - The wheel tip-thickness check always uses 0.4*mn and ignores saCoef
%     (see the NOTE at constraints 7-8).
%   - z1 is a continuous gene rounded on evaluation, so the search is not a true
%     mixed-integer method.

clear; clc; close all;

%% Reproducibility
% Fixed seed, set before the population is initialised, so the Pareto front
% is repeatable run to run (the original single-file version set no seed).
% rng is supported by both MATLAB and GNU Octave.
rng(42);

%% Problem parameters
% Operating point and material data. Values follow the reference paper.
T        = 286.47;   % pinion torque (N*m); 2000*T below gives 2*T in N*mm
K1       = 1;        % pinion load factor (-); 1 = no KA/KV/KH factors, i.e. nominal stress
K2       = 1;        % wheel load factor (-); same assumption as K1
ZE       = 189.8;    % elasticity factor (sqrt(MPa)); 189.8 = steel on steel, E = 206 GPa, nu = 0.3
ZH       = 2.5;      % zone factor (-); 2.5 = 20 deg pressure angle, unshifted spur pair
                     %   (unshifted pairs: about 2.5 at beta = 0 down to 2.4 at beta = 20 deg,
                     %   lower with positive shift); held constant here
ha       = 1.25;     % addendum coefficient (-); NOTE ISO 53 basic rack uses 1.0, and 1.25
                     %   is normally the dedendum coefficient (ha* + c*). 1.25 is kept as in
                     %   the reference; it enlarges da, which thins the tip land and raises
                     %   the contact ratio.
cn       = 0.25;     % bottom clearance coefficient (-); ISO 53 standard value 0.25
pf       = 0.3;      % root fillet radius coefficient rho_f/mn (-); ISO 53 rack 0.38, 0.25-0.4
                     %   common. Not used: tooth_factors hard-codes rhoF = 0.3
sigmaFP1 = 450;      % pinion permissible root bending stress (MPa); value from the reference
sigmaFP2 = 450;      % wheel permissible root bending stress (MPa)
sigmaHP1 = 1300;     % pinion permissible contact stress (MPa); case-hardened steel ~1200-1500
sigmaHP2 = 1300;     % wheel permissible contact stress (MPa)
a        = 180;      % required centre distance (mm)
b        = 60;       % face width (mm)
u        = 3.72;     % nominal gear ratio z2/z1 (-); z2 is rounded, so the actual ratio varies

%% NSGA-II settings
% NOTE saCoef is the minimum tip-land thickness as a multiple of mn. Use 0.25*mn
% for through-hardened (quenched and tempered) gears and 0.4*mn for case-hardened
% (carburised, nitrided) gears, whose hard case would chip on a thin tip. The value
% directly decides which designs are feasible: raising it to 0.4 removes the
% strongly positive-shifted pinions (thin tips) from the feasible set, so the
% Pareto front itself changes, not just its labelling.
saCoef  = 0.25;   % tip-thickness criterion coefficient (-): 0.25 through-hardened, 0.4 case-hardened
popSize = 200;    % population size (-); 100-200 is typical for 2 objectives, 5 variables
maxGen  = 800;    % number of generations (-); budget set by the reference comparison
pc      = 0.9;    % SBX crossover probability (-); Deb's recommended value
pm      = 0.1;    % per-gene mutation probability (-); the 1/numVar rule would give 0.2
etaC    = 20;     % SBX distribution index (-); 10-20 common, larger = children closer to parents
etaM    = 20;     % polynomial mutation distribution index (-); 20 is the usual default

%% Design variables and bounds
% x = [z1 (-), mn (mm), beta (rad), xn1 (-), xn2 (-)]
% z1 >= 17 keeps an unshifted 20 deg pinion free of undercut; beta 8-20 deg is the
% usual helix range (below 8 deg little overlap gain, above 20 deg high axial load).
lb     = [17, 2, 8*pi/180, 0.2, -0.5];    % lower bounds
ub     = [30, 4, 20*pi/180, 0.5, 0.5];    % upper bounds
numVar = length(lb);                      % number of design variables

%% Initial population
% Uniform sampling inside the box. z1 is rounded straight away so every
% individual represents a manufacturable tooth count.
pop = zeros(popSize, numVar);
for i = 1:popSize
    for j = 1:numVar
        pop(i,j) = lb(j) + rand()*(ub(j) - lb(j));
    end
    pop(i,1) = round(pop(i,1));
end

%% NSGA-II main loop
% Each generation: evaluate parents, create offspring, then keep the best popSize
% of parents + offspring (elitist mu+lambda survival), ranked by constraint
% violation, then Pareto rank, then crowding distance.
for gen = 1:maxGen
    [obj, constr] = evaluate_population(pop, T, K1, K2, ZH, ZE, u, b, sigmaFP1, sigmaFP2, sigmaHP1, sigmaHP2, ha, cn, a, saCoef);

    % Total violation = sum of the positive parts of all g(x) <= 0 constraints.
    % The constraints are in mixed units (MPa, mm, -); they are not normalised,
    % so the stress constraints dominate the sum when they are violated.
    constrViolation = zeros(size(pop, 1), 1);
    for i = 1:size(pop, 1)
        constrViolation(i) = sum(max(0, constr(i,:)));
    end

    % Rank and crowding of the parents drive the mating tournament.
    [rank, ~] = fast_non_dominated_sort(obj, constrViolation);
    crowdDist = crowding_distance(obj, rank);
    offspring = genetic_operators(pop, rank, crowdDist, pc, pm, etaC, etaM, lb, ub);

    [objOffspring, constrOffspring] = evaluate_population(offspring, T, K1, K2, ZH, ZE, u, b, sigmaFP1, sigmaFP2, sigmaHP1, sigmaHP2, ha, cn, a, saCoef);
    constrViolationOffspring = zeros(size(offspring, 1), 1);
    for i = 1:size(offspring, 1)
        constrViolationOffspring(i) = sum(max(0, constrOffspring(i,:)));
    end

    % Elitism: parents compete with their children, so a good design is never lost.
    combinedPop             = [pop; offspring];
    combinedObj             = [obj; objOffspring];
    combinedConstrViolation = [constrViolation; constrViolationOffspring];

    [rankCombined, frontsCombined] = fast_non_dominated_sort(combinedObj, combinedConstrViolation);
    crowdDistCombined = crowding_distance(combinedObj, rankCombined);

    nextIdx = select_next_population(rankCombined, crowdDistCombined, popSize);
    pop = combinedPop(nextIdx, :);

    if mod(gen, 100) == 0
        fprintf('第 %d 代迭代完成\n', gen);
    end
end

%% Evaluate the final population
[obj, constr] = evaluate_population(pop, T, K1, K2, ZH, ZE, u, b, sigmaFP1, sigmaFP2, sigmaHP1, sigmaHP2, ha, cn, a, saCoef);
constrViolation = zeros(size(pop, 1), 1);
for i = 1:size(pop, 1)
    constrViolation(i) = sum(max(0, constr(i,:)));
end

%% Extract the Pareto front
% Rank 1 alone is not enough: if no individual were feasible, rank 1 would hold
% the least-infeasible ones. The violation filter keeps only feasible designs;
% 1e-6 absorbs floating-point round-off in the constraint sums.
[rank, fronts] = fast_non_dominated_sort(obj, constrViolation);
paretoOptimal = pop(rank == 1 & constrViolation < 1e-6, :);
paretoFront   = obj(rank == 1 & constrViolation < 1e-6, :);

%% Print the Pareto set
% Columns: z1 (-), mn (mm), beta (deg), xn1 (-), xn2 (-), f1 (MPa), f2 (MPa)
fprintf('\nPareto最优解:\n');
fprintf('z1\t\tmn\t\tbeta\t\txn1\t\txn2\t\tf1\t\tf2\n');
for i = 1:size(paretoOptimal, 1)
    fprintf('%d\t\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', round(paretoOptimal(i,1)), ...
        paretoOptimal(i,2), paretoOptimal(i,3)*180/pi, paretoOptimal(i,4), ...
        paretoOptimal(i,5), paretoFront(i,1), paretoFront(i,2));
end

%% Feasibility check of the baseline design
% The paper's design is evaluated with the same model before it is plotted.
% An earlier version drew it on the front without checking, which suggested it
% was a valid reference point when under this model it may violate constraints.
original_design = [18, 4, 17.27*pi/180, 0.2, 0.3115];   % [z1, mn (mm), beta (rad), xn1, xn2] from the paper
[objOrig, constrOrig] = evaluate_population(original_design, T, K1, K2, ZH, ZE, ...
    u, b, sigmaFP1, sigmaFP2, sigmaHP1, sigmaHP2, ha, cn, a, saCoef);
violOrig     = sum(max(0, constrOrig));
origFeasible = violOrig < 1e-6;

fprintf('\n--- 原始设计校验 ---\n');
fprintf('f1 = %.4f MPa,  f2 = %.4f MPa\n', objOrig(1), objOrig(2));
if origFeasible
    fprintf('原始设计可行\n');
else
    fprintf('原始设计不满足约束，违反度 = %.4f\n', violOrig);
    badIdx = find(constrOrig > 0);
    fprintf('违反的约束编号: %s\n', mat2str(badIdx));
end

%% Plot the Pareto front
% x-axis f2 (bending-stress gap, MPa), y-axis f1 (contact stress, MPa).
figure('Color','w','Position',[100 100 900 650]);
hold on; box on; grid on;

scatter(paretoFront(:,2), paretoFront(:,1), 45, [0.20 0.45 0.70], 'filled', ...
    'MarkerFaceAlpha', 0.75, 'DisplayName', 'Pareto前沿解');

% Representatives: the two extremes, plus the point closest to the ideal
% (min f1, min f2) after scaling each objective to [0, 1]. Scaling matters
% because f1 is ~1000 MPa and f2 is ~10-100 MPa; unscaled, f1 would decide alone.
[~, idxA] = min(paretoFront(:,2));   % A: smallest bending-stress gap
[~, idxC] = min(paretoFront(:,1));   % C: smallest contact stress
dn = (paretoFront - min(paretoFront)) ./ max(max(paretoFront) - min(paretoFront), eps);
[~, idxB] = min(sum(dn.^2, 2));      % B: normalised compromise (eps guards a one-point front)

scatter(paretoFront(idxA,2), paretoFront(idxA,1), 130, 'd', 'filled', ...
    'MarkerFaceColor', [0.85 0.33 0.10], 'DisplayName', 'A  弯曲应力差最小');
scatter(paretoFront(idxB,2), paretoFront(idxB,1), 130, 's', 'filled', ...
    'MarkerFaceColor', [0.47 0.67 0.19], 'DisplayName', 'B  综合折中');
scatter(paretoFront(idxC,2), paretoFront(idxC,1), 130, 'o', 'filled', ...
    'MarkerFaceColor', [0.00 0.30 0.55], 'DisplayName', 'C  接触应力最小');

text(paretoFront(idxA,2), paretoFront(idxA,1), '  A', 'FontSize', 13, 'FontWeight', 'bold');
text(paretoFront(idxB,2), paretoFront(idxB,1), '  B', 'FontSize', 13, 'FontWeight', 'bold');
text(paretoFront(idxC,2), paretoFront(idxC,1), '  C', 'FontSize', 13, 'FontWeight', 'bold');

% An infeasible baseline is drawn hollow and grey so it is not read as a
% competitor to the front.
if origFeasible
    scatter(objOrig(2), objOrig(1), 160, 'p', 'filled', ...
        'MarkerFaceColor', [0.80 0.00 0.00], 'DisplayName', '原始设计');
    text(objOrig(2), objOrig(1), '  原始设计', 'FontSize', 12);
else
    scatter(objOrig(2), objOrig(1), 160, 'p', 'MarkerEdgeColor', [0.55 0.55 0.55], ...
        'LineWidth', 1.4, 'DisplayName', '原始设计（不满足约束）');
    text(objOrig(2), objOrig(1), '  原始设计（不可行）', 'FontSize', 12, ...
        'Color', [0.45 0.45 0.45]);
end

xlabel('f_2  齿根最大弯曲应力差值 (MPa)', 'FontSize', 13);
ylabel('f_1  齿面接触应力 (MPa)', 'FontSize', 13);
title(sprintf('斜齿轮传动几何参数多目标优化 Pareto前沿（NSGA-II，%d个非支配解）', ...
    size(paretoFront,1)), 'FontSize', 15);
legend('Location', 'northeast', 'FontSize', 11, 'Box', 'on');
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;

exportgraphics(gcf, 'pareto_front.png', 'Resolution', 200);

%% Print the balanced design
% NOTE this uses the raw (unscaled) distance to the ideal point, so f1 (MPa in
% the thousands) is weighted more heavily than f2, and the design printed here
% can differ from the normalised point B drawn in the figure.
[~, minIdx] = min(sqrt((paretoFront(:,1) - min(paretoFront(:,1))).^2 + (paretoFront(:,2) - min(paretoFront(:,2))).^2));
optimalSolution   = paretoOptimal(minIdx, :);
optimalObjectives = paretoFront(minIdx, :);

fprintf('\n平衡两个目标的最优解:\n');
fprintf('z1 = %d\n', round(optimalSolution(1)));
fprintf('mn = %.4f mm\n', optimalSolution(2));
fprintf('beta = %.4f 度\n', optimalSolution(3)*180/pi);
fprintf('xn1 = %.4f\n', optimalSolution(4));
fprintf('xn2 = %.4f\n', optimalSolution(5));
fprintf('接触应力 = %.4f MPa\n', optimalObjectives(1));
fprintf('弯曲应力差值 = %.4f MPa\n', optimalObjectives(2));
