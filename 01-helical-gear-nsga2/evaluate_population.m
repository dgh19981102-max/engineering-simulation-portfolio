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
