%%Switch Phase Space
clear; 
clc;

%% Constants
Kd1       = 115.6;
KdT       = 2.1;
Kd2_fixed = 105;
Lvec      = logspace(log10(0.1), log10(1e5), 80).';

opts_fsolve = optimoptions('fsolve','Display','off','TolFun',1e-10,'TolX',1e-10);
opts_lsq    = optimoptions('lsqcurvefit','Display','off');

%% Panel A 
R0  = 960;
T0  = 800;
P0  = 8;

Kd2_vec  = logspace(log10(0.021), log10(210), 25).';
EC50_vec = nan(size(Kd2_vec));
Hill_vec = nan(size(Kd2_vec));
Dyn_vec  = nan(size(Kd2_vec));

for j = 1:numel(Kd2_vec)
    params.R0   = R0;
    params.T0   = T0;
    params.P0   = P0;
    params.Kd1  = Kd1;
    params.Kd2  = Kd2_vec(j);
    params.KdT  = KdT;
    params.Cmax = min(R0, P0);
    [Dyn_vec(j), EC50_vec(j), Hill_vec(j)] = compute_metrics(params, Lvec, opts_fsolve, opts_lsq);
end

Kd2_norm = Kd2_vec ./ KdT;

figure;
loglog(Kd2_norm, EC50_vec, '-o','LineWidth',2);
xlabel('K_{d2} / K_{dT}');
ylabel('EC_{50}');
title('EC_{50} vs K_{d2}/K_{dT}');
grid on;

figure;
semilogx(Kd2_norm, Hill_vec, '-o','LineWidth',2);
xlabel('K_{d2} / K_{dT}');
ylabel('Hill coefficient');
title('Hill vs K_{d2}/K_{dT}');
grid on;

figure;
semilogx(Kd2_norm, Dyn_vec, '-o','LineWidth',2);
xlabel('K_{d2} / K_{dT}');
ylabel('Dynamic range (max C_{eq} / P_0)');
title('Dynamic range vs K_{d2}/K_{dT}');
grid on;

Results1 = table(Kd2_norm, EC50_vec, Hill_vec, Dyn_vec, ...
    'VariableNames', {'Kd2_over_KdT','EC50','HillCoeff','DynamicRange'});
writetable(Results1,'ultrasensitivity_sweep.csv');

%% Panel B
P0  = 8;
Kd2 = Kd2_fixed;

nR   = 25;
nT   = 25;
R_vec = logspace(log10(P0/100), log10(100*P0), nR).';
T_vec = logspace(log10(P0/100), log10(100*P0), nT).';

EC50_RT = nan(nT, nR);
Hill_RT = nan(nT, nR);
Dyn_RT  = nan(nT, nR);

for iR = 1:nR
    R0 = R_vec(iR);
    for iT = 1:nT
        T0 = T_vec(iT);
        params.R0   = R0;
        params.T0   = T0;
        params.P0   = P0;
        params.Kd1  = Kd1;
        params.Kd2  = Kd2;
        params.KdT  = KdT;
        params.Cmax = min(R0, P0);
        [Dyn_RT(iT,iR), EC50_RT(iT,iR), Hill_RT(iT,iR)] = compute_metrics(params, Lvec, opts_fsolve, opts_lsq);
    end
end

[R_grid, T_grid] = meshgrid(R_vec/P0, T_vec/P0);

figure;
contourf(R_grid, T_grid, Dyn_RT, 20,'LineStyle','none');
set(gca,'XScale','log','YScale','log');
colorbar;
xlabel('[R]_0 / [P]_0');
ylabel('[T]_0 / [P]_0');
title('Dynamic range (R,T)');

figure;
contourf(R_grid, T_grid, EC50_RT, 20,'LineStyle','none');
set(gca,'XScale','log','YScale','log');
colorbar;
xlabel('[R]_0 / [P]_0');
ylabel('[T]_0 / [P]_0');
title('EC_{50} (R,T)');

figure;
contourf(R_grid, T_grid, Hill_RT, 20,'LineStyle','none');
set(gca,'XScale','log','YScale','log');
colorbar;
xlabel('[R]_0 / [P]_0');
ylabel('[T]_0 / [P]_0');
title('Hill (R,T)');

R_flat    = R_grid(:);
T_flat    = T_grid(:);
Dyn_flat  = Dyn_RT(:);
EC50_flat = EC50_RT(:);
Hill_flat = Hill_RT(:);

Results2 = table(R_flat, T_flat, Dyn_flat, EC50_flat, Hill_flat, ...
    'VariableNames', {'R_over_P','T_over_P','DynamicRange','EC50','HillCoeff'});
writetable(Results2,'RT_contour_Kd2_105nM.csv');

Rnorm = R_vec / P0;
Tnorm = T_vec / P0;

DynTable_RT = array2table(Dyn_RT, ...
    'VariableNames', compose('R_over_P=%g', Rnorm), ...
    'RowNames',      compose('T_over_P=%g', Tnorm));
writetable(DynTable_RT,'DynamicRange_RT_matrix.csv','WriteRowNames',true);

EC50Table_RT = array2table(EC50_RT, ...
    'VariableNames', compose('R_over_P=%g', Rnorm), ...
    'RowNames',      compose('T_over_P=%g', Tnorm));
writetable(EC50Table_RT,'EC50_RT_matrix.csv','WriteRowNames',true);

HillTable_RT = array2table(Hill_RT, ...
    'VariableNames', compose('R_over_P=%g', Rnorm), ...
    'RowNames',      compose('T_over_P=%g', Tnorm));
writetable(HillTable_RT,'Hill_RT_matrix.csv','WriteRowNames',true);

%% Panel C 
T0  = 800;
nP  = 25;
nR2 = 25;

P_vec = logspace(log10(T0/100), log10(100*T0), nP).';    % partner
R_vec2 = logspace(log10(T0/100), log10(100*T0), nR2).';  % receptor

EC50_PR = nan(nP, nR2);
Hill_PR = nan(nP, nR2);
Dyn_PR  = nan(nP, nR2);

for iP = 1:nP
    P0 = P_vec(iP);
    for iR = 1:nR2
        R0 = R_vec2(iR);
        params.R0   = R0;
        params.T0   = T0;
        params.P0   = P0;
        params.Kd1  = Kd1;
        params.Kd2  = Kd2;
        params.KdT  = KdT;
        params.Cmax = min(R0, P0);
        [Dyn_PR(iP,iR), EC50_PR(iP,iR), Hill_PR(iP,iR)] = compute_metrics(params, Lvec, opts_fsolve, opts_lsq);
    end
end

Rnorm2 = R_vec2 / T0;
Pnorm2 = P_vec  / T0;
[R_grid2, P_grid2] = meshgrid(Rnorm2, Pnorm2);

figure;
contourf(R_grid2, P_grid2, Dyn_PR, 15,'LineStyle','none');
set(gca,'XScale','log','YScale','log');
colorbar;
xlabel('[R]_0 / [T]_0');
ylabel('[P]_0 / [T]_0');
title('Dynamic range (P,R)');

figure;
contourf(R_grid2, P_grid2, EC50_PR, 15,'LineStyle','none');
set(gca,'XScale','log','YScale','log');
colorbar;
xlabel('[R]_0 / [T]_0');
ylabel('[P]_0 / [T]_0');
title('EC_{50} (P,R)');

figure;
contourf(R_grid2, P_grid2, Hill_PR, 15,'LineStyle','none');
set(gca,'XScale','log','YScale','log');
colorbar;
xlabel('[R]_0 / [T]_0');
ylabel('[P]_0 / [T]_0');
title('Hill (P,R)');

DynTable_PR = array2table(Dyn_PR, ...
    'VariableNames', compose('R_over_T=%.3g', Rnorm2), ...
    'RowNames',      compose('P_over_T=%.3g', Pnorm2));
writetable(DynTable_PR,'DynamicRange_PR_matrix.csv','WriteRowNames',true);

EC50Table_PR = array2table(EC50_PR, ...
    'VariableNames', compose('R_over_T=%.3g', Rnorm2), ...
    'RowNames',      compose('P_over_T=%.3g', Pnorm2));
writetable(EC50Table_PR,'EC50_PR_matrix.csv','WriteRowNames',true);

HillTable_PR = array2table(Hill_PR, ...
    'VariableNames', compose('R_over_T=%.3g', Rnorm2), ...
    'RowNames',      compose('P_over_T=%.3g', Pnorm2));
writetable(HillTable_PR,'Hill_PR_matrix.csv','WriteRowNames',true);

%% Panel D 
R0  = 960;
nP3 = 25;
nT2 = 25;

P_vec3 = logspace(log10(R0/100), log10(100*R0), nP3).';  % partner
T_vec2 = logspace(log10(R0/100), log10(100*R0), nT2).';  % titrant

EC50_PT = nan(nP3, nT2);
Hill_PT = nan(nP3, nT2);
Dyn_PT  = nan(nP3, nT2);

for iP = 1:nP3
    P0 = P_vec3(iP);
    for iT = 1:nT2
        T0 = T_vec2(iT);
        params.R0   = R0;
        params.T0   = T0;
        params.P0   = P0;
        params.Kd1  = Kd1;
        params.Kd2  = Kd2;
        params.KdT  = KdT;
        params.Cmax = min(R0, P0);
        [Dyn_PT(iP,iT), EC50_PT(iP,iT), Hill_PT(iP,iT)] = compute_metrics(params, Lvec, opts_fsolve, opts_lsq);
    end
end

Tnorm3 = T_vec2 / R0;
Pnorm3 = P_vec3 / R0;
[T_grid3, P_grid3] = meshgrid(Tnorm3, Pnorm3);

figure;
contourf(T_grid3, P_grid3, Dyn_PT, 15,'LineStyle','none');
set(gca,'XScale','log','YScale','log');
colorbar;
xlabel('[T]_0 / [R]_0');
ylabel('[P]_0 / [R]_0');
title('Dynamic range (P,T)');

figure;
contourf(T_grid3, P_grid3, EC50_PT, 15,'LineStyle','none');
set(gca,'XScale','log','YScale','log');
colorbar;
xlabel('[T]_0 / [R]_0');
ylabel('[P]_0 / [R]_0');
title('EC_{50} (P,T)');

figure;
contourf(T_grid3, P_grid3, Hill_PT, 15,'LineStyle','none');
set(gca,'XScale','log','YScale','log');
colorbar;
xlabel('[T]_0 / [R]_0');
ylabel('[P]_0 / [R]_0');
title('Hill (P,T)');

DynTable_PT = array2table(Dyn_PT, ...
    'VariableNames', compose('T_over_R=%.3g', Tnorm3), ...
    'RowNames',      compose('P_over_R=%.3g', Pnorm3));
writetable(DynTable_PT,'DynamicRange_PT_matrix.csv','WriteRowNames',true);

EC50Table_PT = array2table(EC50_PT, ...
    'VariableNames', compose('T_over_R=%.3g', Tnorm3), ...
    'RowNames',      compose('P_over_R=%.3g', Pnorm3));
writetable(EC50Table_PT,'EC50_PT_matrix.csv','WriteRowNames',true);

HillTable_PT = array2table(Hill_PT, ...
    'VariableNames', compose('T_over_R=%.3g', Tnorm3), ...
    'RowNames',      compose('P_over_R=%.3g', Pnorm3));
writetable(HillTable_PT,'Hill_PT_matrix.csv','WriteRowNames',true);

function [Dyn, EC50, Hill] = compute_metrics(p, Lvec, opts_fsolve, opts_lsq)
    nL  = numel(Lvec);
    Ceq = nan(nL,1);
    for k = 1:nL
        L0 = Lvec(k);
        if k == 1 || isnan(Ceq(k-1))
            Cguess = 0.01 * p.Cmax;
        else
            Cguess = Ceq(k-1);
        end
        fun = @(C) ligand_balance_eq21(C, L0, p);
        [Csol,~,exitflag] = fsolve(fun, Cguess, opts_fsolve);
        if exitflag > 0 && Csol > 0 && Csol < p.Cmax
            Ceq(k) = Csol;
        else
            Ceq(k) = NaN;
        end
    end
    Dyn = max(Ceq) / p.P0;
    valid = ~isnan(Ceq);
    xdata = Lvec(valid);
    ydata = Ceq(valid);
    if numel(xdata) < 10
        EC50 = NaN; Hill = NaN; return;
    end
    bottom0 = min(ydata);
    top0    = max(ydata);
    EC50_0  = 10^(mean(log10(xdata)));
    nH_0    = 2;
    p0 = [bottom0, top0, EC50_0, nH_0];
    lb = [0, 0, min(xdata), 0.1];
    ub = [Inf, Inf, max(xdata), 10];
    pfit = lsqcurvefit(@hill_fun, p0, xdata, ydata, lb, ub, opts_lsq);
    EC50 = pfit(3);
    Hill = pfit(4);
end

function F = hill_fun(p, x)
    F = p(1) + (p(2)-p(1)) ./ (1 + (p(3)./x).^p(4));
end

function F = ligand_balance_eq21(C, L0, p)
    if C <= 0 || C >= p.Cmax
        F = 1e6; 
        return;
    end
    alpha = p.Kd2 .* C ./ (p.P0 - C);                   % α = Kd2·C/(P0−C)
    denom = p.R0 - C - alpha - (p.T0 .* alpha) ./ (p.KdT + alpha);
    if denom <= 0
        F = 1e6;
        return;
    end
    L_expr = (p.Kd1 .* alpha) ./ denom ...              % [L]_0 from Eq. 21
             + alpha ...
             + (p.T0 .* alpha) ./ (p.KdT + alpha) ...
             + C;
    F = L_expr - L0;
end
