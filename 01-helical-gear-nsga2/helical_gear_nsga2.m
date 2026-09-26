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
%   The population is initialised with rand() and no seed is set, so the front
%   changes slightly between runs. Call rng(<seed>) first for a repeatable run.
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

%% Objective and constraint evaluation
% Returns obj (N x 2, MPa) and constr (N x 10). Every constraint is written as
% g(x) <= 0, so a positive entry is the amount of violation in that row's unit.
function [obj, constr] = evaluate_population(pop, T, K1, K2, ZH, ZE, u, b, sigmaFP1, sigmaFP2, sigmaHP1, sigmaHP2, ha, cn, a, saCoef)
    popSize = size(pop, 1);
    obj     = zeros(popSize, 2);
    constr  = zeros(popSize, 10);

    for i = 1:popSize
        z1   = round(pop(i,1));
        mn   = pop(i,2);    % normal module (mm)
        beta = pop(i,3);    % helix angle at the reference circle (rad)
        xn1  = pop(i,4);    % pinion profile shift coefficient (-)
        xn2  = pop(i,5);    % wheel profile shift coefficient (-)

        % Tooth counts must be integers, so the real ratio z2/z1 deviates slightly from u.
        z2     = round(u * z1);
        alphan = 20 * pi/180;                  % normal pressure angle (rad); 20 deg ISO standard
        alphat = atan(tan(alphan)/cos(beta));  % transverse pressure angle (rad): tooth profile
                                               %   seen in the plane of rotation

        [YF1, YS1] = tooth_factors(z1, beta, alphan, xn1);
        [YF2, YS2] = tooth_factors(z2, beta, alphan, xn2);

        % Helix factor (-): inclined contact lines spread the root load, lowering
        % bending stress. ISO 6336-3 Ybeta = 1 - eps_beta*beta/120deg with
        % eps_beta capped at 1; here eps_beta = 1 is assumed.
        Ybeta = 1 - beta/(120*pi/180);

        % Objective 1: flank contact stress (MPa), Hertzian line contact.
        % 2000*T/(mn*z1/cos(beta)) is the tangential force in N (T in N*m,
        % reference diameter d1 = mn*z1/cos(beta) in mm). The cos(beta)^3 term
        % collects cos(beta)^2 from d1^2 and the helix factor Zbeta^2 = cos(beta).
        % NOTE the extra 1/u under the first root is not in the ISO 6336-2 form
        % sqrt(Ft/(b*d1) * (u+1)/u); it is reproduced from the reference and lowers
        % sigmah by sqrt(u). Check it before using the absolute values.
        sigmah = ZH * ZE * sqrt((2000 * T * cos(beta)^3) / (mn^2 * z1^2 * b * u)) * sqrt((u+1)/u) * K1;
        obj(i,1) = sigmah;

        % Objective 2: gap between pinion and wheel root bending stress (MPa).
        % Lewis-type root stress sigmaF = Ft*YF*YS*Ybeta/(b*mn), with
        % Ft = 2000*T*cos(beta)/(mn*z1) in N. Minimising the gap balances the
        % bending life of the two gears, which is what profile shift is used for.
        % NOTE sigmaF2 divides by z2, but both gears carry the same tangential
        % force, which is set by the pinion (z1). As written, sigmaF2 is too low by
        % a factor of u = 3.72, so the wheel always looks lightly loaded in bending.
        sigmaF1 = (2000 * T * K1 * YF1 * YS1 * Ybeta * cos(beta)) / (mn^2 * z1 * b);
        sigmaF2 = (2000 * T * K2 * YF2 * YS2 * Ybeta * cos(beta)) / (mn^2 * z2 * b);
        obj(i,2) = abs(sigmaF1 - sigmaF2);

        % Constraints 1-2: contact fatigue (MPa). Both gears see the same Hertzian
        % stress, so these only differ if the two permissible stresses differ.
        constr(i,1) = sigmah - sigmaHP1;
        constr(i,2) = sigmah - sigmaHP2;

        % Constraints 3-4: root bending fatigue (MPa).
        constr(i,3) = sigmaF1 - sigmaFP1;
        constr(i,4) = sigmaF2 - sigmaFP2;

        % Constraint 5: centre distance (mm).
        % NOTE the centre distance must be solved through the working pressure
        % angle alphat' (atp). Shifting the profiles changes the tooth thickness
        % on the reference circle, so the pair only meshes without backlash at the
        % pressure angle given by
        %     inv(alphat') = inv(alphat) + 2*(xn1 + xn2)*tan(alphan)/(z1 + z2)
        % and the centre distance is then a' = a_std*cos(alphat)/cos(alphat').
        % If alphat is used in place of alphat', the ratio becomes 1 and
        % xn1 + xn2 drops out of the equation. The check then treats every shifted
        % pair as if it sat at the unshifted distance, so no shifted pair can
        % satisfy the constraint through its shift and the search piles up
        % against the helix-angle bound instead.
        a_std    = mn * (z1 + z2) / (2 * cos(beta));   % reference centre distance, no shift (mm)
        inv_at   = tan(alphat) - alphat;               % involute function inv(x) = tan(x) - x (-)
        inv_atp  = inv_at + 2*(xn1 + xn2)*tan(alphan)/(z1 + z2);
        atp      = solve_inv(inv_atp, alphat);         % working transverse pressure angle (rad)
        a_actual = a_std * cos(alphat) / cos(atp);     % working centre distance (mm)
        % Tolerance 0.5 mm (was 0.01 mm). With a continuous mn and beta the
        % working distance can never be hit exactly; 0.5 mm is the band a small
        % housing-bore or backlash allowance can take up.
        constr(i,5) = abs(a_actual - a) - 0.5;

        % Constraint 6: total contact ratio eps_gamma >= 2 (-), i.e. at least two
        % tooth pairs share the load at all times, for smooth running and lower
        % noise. eps_beta is the overlap ratio contributed by the helix.
        ep_alpha = calculate_transverse_contact_ratio(z1, z2, alphan, beta, mn, xn1, xn2);
        ep_beta  = (b * tan(beta)) / (pi * mn);
        ep_gamma = ep_alpha + ep_beta;
        constr(i,6) = 2 - ep_gamma;

        % Constraints 7-8: tip-land thickness (mm). Positive shift thickens the
        % root but sharpens the tip; a pointed tip chips or through-hardens.
        % s_a = d_a*[(pi + 4*x*tan(alphan))/(2*z) + inv(alpha) - inv(alpha_a)]
        % is the tooth thickness transferred from the reference circle to the tip
        % circle. NOTE spur-gear form: diameters use mn*z (not mn*z/cos(beta)) and
        % alpha_a uses the normal pressure angle, so this approximates the
        % thickness in the normal section.
        da1          = mn * (z1 + 2 * ha + 2 * xn1);   % tip diameter (mm)
        da2          = mn * (z2 + 2 * ha + 2 * xn2);
        alpha_a1     = acos(z1 * cos(alphan) / (z1 + 2 * ha + 2 * xn1));   % pressure angle at the tip (rad)
        alpha_a2     = acos(z2 * cos(alphan) / (z2 + 2 * ha + 2 * xn2));
        inv_alphat   = tan(alphat) - alphat;
        inv_alpha_a1 = tan(alpha_a1) - alpha_a1;
        inv_alpha_a2 = tan(alpha_a2) - alpha_a2;

        s_a1 = da1 * ((pi + 4 * xn1 * tan(alphan)) / (2 * z1) + inv_alphat - inv_alpha_a1);
        s_a2 = da2 * ((pi + 4 * xn2 * tan(alphan)) / (2 * z2) + inv_alphat - inv_alpha_a2);

        % NOTE saCoef (0.25*mn through-hardened, 0.4*mn case-hardened) only acts on
        % the pinion. The wheel row is hard-coded to 0.4*mn, so the wheel is always
        % checked against the case-hardened limit whatever saCoef is set to.
        constr(i,7) = saCoef * mn - s_a1;
        constr(i,8) = 0.4 * mn - s_a2;

        % Constraints 9-10: minimum profile shift to avoid undercut (-).
        % x_min = ha - z*sin(alphat)^2/(2*cos(beta)); below this the generating
        % rack tip cuts into the involute near the base circle.
        constr(i,9)  = (ha - z1 * sin(alphat)^2 / (2 * cos(beta))) - xn1;
        constr(i,10) = (ha - z2 * sin(alphat)^2 / (2 * cos(beta))) - xn2;
    end
end

%% Form factor and stress-correction factor
% Returns YF (-) and YS (-) for one gear. Lengths are in units of mn
% (dimensionless), so the module cancels out.
% NOTE hFe and sFn are empirical fits from the reference, not the ISO 6336-3
% 30-degree tangent construction; YF is only indicative. zv and alphat are
% computed but not used (the fit takes z, not the virtual tooth count).
function [YF, YS] = tooth_factors(z, beta, alphan, xn)
    alphat = atan(tan(alphan)/cos(beta));
    zv     = z / cos(beta)^3;   % virtual spur tooth count (-) of the normal section

    % Bending moment arm hFe/mn (-): shrinks with more teeth (stubbier profile)
    % and grows with positive shift (load applied further from the root).
    hFe = 0.8 * z^(-0.5) + 0.5 * xn + 0.15;

    % Root chord sFn/mn (-): pi/2 is half the reference pitch; positive shift
    % adds 2*x*tan(alphan) of thickness, which is why shift strengthens the root.
    sFn = (pi/2 + 2 * xn * tan(alphan)) * cos(alphan);

    % Load angle at the tip (rad); taken equal to alphan as a simplification.
    alphaFen = alphan;

    % Form factor (-): cantilever bending stress 6*M/(b*s^2) at the root chord.
    YF = (6 * hFe / cos(alphaFen)) * (sFn / cos(alphan))^(-2);

    % Root fillet radius rho_F/mn (-). 0.3 is assumed; the ISO 53 rack gives 0.38.
    rhoF = 0.3;

    % Notch parameter qs = sFn/(2*rhoF) (-): smaller fillet = sharper notch.
    qs = sFn / (2 * rhoF);

    % Ratio of root chord to moment arm (-).
    L = sFn / hFe;

    % Stress-correction factor (-), ISO 6336-3 eq. for YS; valid for 1 <= qs < 8.
    YS = (1.2 + 0.13 * L) * qs^(1 / (1.21 + 2.3/L));
end

%% Transverse contact ratio
% eps_alpha (-) = length of the path of contact / base pitch.
% NOTE spur-gear relations with the normal module: d = z*mn instead of
% z*mn/cos(beta), base circle and base pitch from alphan instead of alphat, and the
% reference (not working) centre distance. The addendum is re-declared locally as
% 1.25 and so ignores the ha passed to the caller.
function ep_alpha = calculate_transverse_contact_ratio(z1, z2, alphan, beta, mn, xn1, xn2)
    alphat = atan(tan(alphan)/cos(beta));

    ha  = 1.25;                          % addendum coefficient (-), same value as the script
    da1 = mn * (z1 + 2 * ha + 2 * xn1);  % tip diameters (mm)
    da2 = mn * (z2 + 2 * ha + 2 * xn2);

    db1 = z1 * mn * cos(alphan);         % base diameters (mm)
    db2 = z2 * mn * cos(alphan);

    d1 = z1 * mn;                        % reference diameters (mm)
    d2 = z2 * mn;

    a = (d1 + d2) / 2;                   % reference centre distance (mm)

    ra1 = da1 / 2;                       % radii (mm)
    ra2 = da2 / 2;
    rb1 = db1 / 2;
    rb2 = db2 / 2;

    % Path of contact g_alpha (mm): the part of the line of action between the
    % two tip circles, minus the centre-distance projection a*sin(alpha).
    ga = sqrt(ra1^2 - rb1^2) + sqrt(ra2^2 - rb2^2) - a * sin(alphat);

    % Divide by the base pitch pi*mn*cos(alphan) (mm).
    ep_alpha = ga / (pi * mn * cos(alphan));
end

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

%% Environmental selection
% Sort by (rank ascending, crowding distance descending) and keep the first
% popSize. Equivalent to NSGA-II's front-by-front fill with the last front
% truncated by crowding distance.
function nextIdx = select_next_population(rank, crowdDist, popSize)
    [~, sortedIdx] = sortrows([rank, -crowdDist]);
    nextIdx = sortedIdx(1:popSize);
end

%% Inverse involute function
% Solves tan(ang) - ang = invval for the pressure angle ang (rad) by Newton's
% method. There is no closed-form inverse. d/dx(tan x - x) = tan(x)^2, which
% vanishes at 0, so the start is kept away from 0 and each step is clamped to
% [1e-4, 1.5] rad (1.5 rad stays below the pole of tan at pi/2).
% x0: starting guess (rad); the reference transverse pressure angle is close.
function ang = solve_inv(invval, x0)
    ang = max(x0, 0.01);
    for k = 1:60
        f  = tan(ang) - ang - invval;
        df = tan(ang)^2;
        if abs(df) < 1e-12, break; end
        step = f / df;
        ang  = ang - step;
        ang  = min(max(ang, 1e-4), 1.5);
        if abs(step) < 1e-12, break; end
    end
end
