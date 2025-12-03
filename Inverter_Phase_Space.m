%%Inverter Phase Space
clear; clc;

opts = optimoptions('fsolve','Display','off', ...
    'FunctionTolerance',1e-12, 'StepTolerance',1e-12);

%% Panel B 
Kd1  = 115.6;
Kd2  = 2.0;
R0   = 200;   % Receptor
A0   = 20;    % Activator
P0   = 2;     % Partner

Lvec1 = logspace(-3, 6, 10).';

Kd_ratio_vec = logspace(-2, 2, 40).';
Cap_vs_ratio = nan(numel(Kd_ratio_vec),1);

for j = 1:numel(Kd_ratio_vec)
    Kd_off_j = Kd_ratio_vec(j) * Kd2;
    Cap_vs_ratio(j) = sequestration_cap(Kd1, Kd2, Kd_off_j, R0, A0, P0, Lvec1, opts);
end

figure;
semilogx(Kd_ratio_vec, Cap_vs_ratio, 'o-', 'LineWidth', 1.5);
xlabel('K_{d,off} / K_{d2}');
ylabel('Fractional Inhibition');
ylim([0 1.05]);
grid on;
title('Fractional Inhibition');


%% Panel C (left)

Kd1    = 115.6;
Kd2    = 2.1;
Kd_off = 200;
P0     = 2;

Lvec2 = logspace(-4, 6, 10).';

R_over_P = logspace(-2, 2, 35);
P_over_A = logspace(-2, 2, 35);

R_vals_raw = R_over_P .* P0;
A_vals_raw = P0 ./ P_over_A;

[nY, nX] = deal(numel(P_over_A), numel(R_over_P));
[R_grid_norm, PoverA_grid] = meshgrid(R_over_P, P_over_A);
Cap_RP = nan(nY, nX);

for iy = 1:nY
    A0 = A_vals_raw(iy);
    for ix = 1:nX
        R0 = R_vals_raw(ix);
        Cap_RP(iy, ix) = sequestration_cap(Kd1, Kd2, Kd_off, R0, A0, P0, Lvec2, opts);
    end
end

Cap_RP_filled = fill_nans(Cap_RP, R_vals_raw, A_vals_raw);
Cap_RP_filled = max(0, min(1, Cap_RP_filled));

figure;
contourf(R_grid_norm, PoverA_grid, Cap_RP_filled, 0:0.1:1, 'LineColor','none');
set(gca, 'XScale','log', 'YScale','log');
xlabel('[Receptor] / [Partner]');
ylabel('[Partner] / [Activator]');
title('Fractional Inhibition');
colorbar; caxis([0 1]);
hold on; contour(R_grid_norm, PoverA_grid, Cap_RP_filled, 0:0.1:1, 'k'); hold off;

ExportMat2 = cell(nY+1, nX+1);
ExportMat2{1,1}         = '';
ExportMat2(1,2:end)     = num2cell(R_over_P);
ExportMat2(2:end,1)     = num2cell(P_over_A');
ExportMat2(2:end,2:end) = num2cell(Cap_RP_filled);
writecell(ExportMat2, 'sequestration_phase_space_R_overP_vs_P_overA.csv');


%% Panel C (right)

Kd1    = 115.6;
Kd2    = 2.1;
Kd_off = 200;
A0     = 20;

Lvec3 = Lvec2;

R_over_A  = logspace(-2, 2, 35);
P_over_A2 = logspace(-2, 2, 35);

R_vals_raw2 = A0 .* R_over_A;
P_vals_raw2 = A0 .* P_over_A2;

[nP2, nR2] = deal(numel(P_over_A2), numel(R_over_A));
[R_grid_norm2, PoverA_grid2] = meshgrid(R_over_A, P_over_A2);
Cap_RA = nan(nP2, nR2);

for ip = 1:nP2
    P0 = P_vals_raw2(ip);
    for ir = 1:nR2
        R0 = R_vals_raw2(ir);
        Cap_RA(ip, ir) = sequestration_cap(Kd1, Kd2, Kd_off, R0, A0, P0, Lvec3, opts);
    end
end

Cap_RA_filled = fill_nans(Cap_RA, R_vals_raw2, P_vals_raw2);
Cap_RA_f
