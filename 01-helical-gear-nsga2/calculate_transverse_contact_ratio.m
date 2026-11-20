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
