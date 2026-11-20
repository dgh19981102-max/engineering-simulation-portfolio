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
