function stats = tp1_repere_stator()
% TP1_REPERE_STATOR Manipulation I - Partie I :
% Demarrage direct a vide du moteur asynchrone dans le repere lie au stator (wa = 0).

    script_dir = fileparts(mfilename('fullpath'));
    fig_dir = fullfile(script_dir, '..', 'figures');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

    C = set_pub_style();
    p = parametres_mas();
    h = 5e-4;
    tf = 3.0;

    res = simuler_mas_rk4('stator', p, h, tf);

    % Extraction des grandeurs caracteristiques
    idx_95 = find(res.N_trmin >= 0.95 * res.N_trmin(end), 1);
    stats.t_dem_95   = res.t(idx_95);
    stats.Ce_max     = max(res.Ce);
    stats.Ce_min     = min(res.Ce);
    stats.Ce_fin     = res.Ce(end);
    stats.N_fin      = res.N_trmin(end);
    stats.Omega_fin  = res.Omega(end);
    stats.s_fin_pct  = 100 * (p.Ns - stats.N_fin) / p.Ns;
    stats.Is_max_pic = max(max(abs([res.isa; res.isb; res.isc])));
    stats.isa_max    = max(abs(res.isa));
    stats.Is_eff_max = max(res.Is_eff);
    stats.Is_eff_fin = res.Is_eff(end);
    stats.Ir_eff_max = max(res.Ir_eff);
    stats.Ir_eff_fin = res.Ir_eff(end);
    stats.ira_max    = max(abs(res.ira));
    stats.res        = res;

    % --- Figure 1 : 4 grandeurs principales dans le repere lie au stator ---
    fig1 = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 23 15.5]);
    tl = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, '\textbf{Manipulation I -- D\''emarrage direct \`a vide (Rep\`ere li\''e au stator, $\omega_a = 0$)}', ...
          'FontSize', 12, 'Interpreter', 'latex');

    % (a) Couple electromagnetique
    nexttile(tl);
    plot(res.t, res.Ce, 'Color', C.navy, 'LineWidth', 1.3); hold on;
    yline(stats.Ce_fin, '--', 'Color', C.crimson, 'LineWidth', 1.0);
    [cmax, imax] = max(res.Ce);
    plot(res.t(imax), cmax, 'o', 'Color', C.crimson, 'MarkerFaceColor', C.crimson, 'MarkerSize', 5);
    text(res.t(imax)+0.12, cmax-1.5, sprintf('$C_{e,\\max} = %.2f\\,\\mathrm{N\\cdot m}$', cmax), ...
         'Color', C.crimson, 'FontSize', 9);
    text(1.75, stats.Ce_fin+2.6, sprintf('$C_{e,\\infty} = %.2f\\,\\mathrm{N\\cdot m}$', stats.Ce_fin), ...
         'Color', C.crimson, 'FontSize', 9);
    xlabel('Temps $t$ (s)');
    ylabel('Couple $C_e$ ($\mathrm{N\cdot m}$)');
    title('(a) Couple \''electromagn\''etique $C_e(t)$', 'FontSize', 10);
    xlim([0 tf]); ylim([-5 34]);

    % (b) Vitesse mecanique
    nexttile(tl);
    yyaxis left;
    plot(res.t, res.Omega, 'Color', C.navy, 'LineWidth', 1.5); hold on;
    yline(p.Ws_mec, '--', 'Color', C.slate, 'LineWidth', 1.0);
    xline(stats.t_dem_95, ':', 'Color', C.emerald, 'LineWidth', 1.2);
    ylabel('Vitesse $\Omega$ (rad/s)', 'Color', C.navy);
    ylim([0 175]);
    ax = gca; ax.YAxis(1).Color = C.navy;

    yyaxis right;
    plot(res.t, res.N_trmin, 'Color', C.navy, 'LineStyle', 'none');
    ylabel('Vitesse $N$ (tr/min)', 'Color', C.crimson);
    ylim([0 175 * 60 / (2*pi)]);
    ax.YAxis(2).Color = C.crimson;

    text(stats.t_dem_95+0.08, 450, sprintf('$t_{95\\%%} = %.2f\\,\\mathrm{s}$', stats.t_dem_95), ...
         'Color', C.emerald, 'FontSize', 9);
    text(1.45, 1250, sprintf('$N_{\\infty} = %.1f\\,\\mathrm{tr/min}$ ($s=%.2f\\%%$)', ...
         stats.N_fin, stats.s_fin_pct), 'Color', C.navy, 'FontSize', 9);
    xlabel('Temps $t$ (s)');
    title('(b) Vitesse m\''ecanique $\Omega(t)$ et $N(t)$', 'FontSize', 10);
    xlim([0 tf]);

    % (c) Courant statorique
    nexttile(tl);
    plot(res.t, res.isd, 'Color', [0.043 0.235 0.365 0.45], 'LineWidth', 0.9); hold on;
    plot(res.t, res.isa, 'Color', C.crimson, 'LineWidth', 1.0);
    plot(res.t, sqrt(2)*res.Is_eff, 'Color', C.black, 'LineWidth', 1.3, 'LineStyle', '--');
    xlabel('Temps $t$ (s)');
    ylabel('Courants statoriques (A)');
    title('(c) Courants statoriques $i_{sd}(t)$, $i_{sa}(t)$ et enveloppe', 'FontSize', 10);
    legend({'$i_{sd}$ (Park)', '$i_{sa}$ (phase $a$)', '$\sqrt{2}\,I_{s,\mathrm{eff}}$'}, ...
           'Location', 'northeast', 'FontSize', 8);
    xlim([0 tf]); ylim([-25 25]);

    % (d) Courant rotorique
    nexttile(tl);
    plot(res.t, res.ird, 'Color', [0.118 0.518 0.286 0.40], 'LineWidth', 0.9); hold on;
    plot(res.t, res.ira, 'Color', C.amber, 'LineWidth', 1.15);
    plot(res.t, sqrt(2)*res.Ir_eff, 'Color', C.black, 'LineWidth', 1.2, 'LineStyle', '--');
    xlabel('Temps $t$ (s)');
    ylabel('Courants rotoriques (A)');
    title('(d) Courants rotoriques $i_{rd}(t)$ ($50\,\mathrm{Hz}$) et $i_{ra}(t)$ ($s f_s$)', 'FontSize', 10);
    legend({'$i_{rd}$ (rep\`ere stator)', '$i_{ra}$ (phase rotor)', '$\sqrt{2}\,I_{r,\mathrm{eff}}$'}, ...
           'Location', 'northeast', 'FontSize', 8);
    xlim([0 tf]); ylim([-58 58]);

    % Zoom inset place APRES tous les nexttile pour ne pas perturber tiledlayout
    ax_in1 = axes(fig1, 'Position', [0.27 0.11 0.165 0.115]);
    box on;
    idx_z = res.t >= 2.50 & res.t <= 2.56;
    plot(ax_in1, res.t(idx_z), res.isa(idx_z), 'Color', C.crimson, 'LineWidth', 1.1); hold(ax_in1, 'on');
    plot(ax_in1, res.t(idx_z), res.isd(idx_z), 'Color', C.navy, 'LineWidth', 0.9, 'LineStyle', '--');
    xlim(ax_in1, [2.50 2.56]); ylim(ax_in1, [-6.5 6.5]);
    title(ax_in1, 'Zoom 50 Hz', 'FontSize', 7.5);
    set(ax_in1, 'FontSize', 7);

    exportgraphics(fig1, fullfile(fig_dir, 'fig_tp1_stator_general.png'), 'Resolution', 300);
    close(fig1);

    % --- Figure 2 : Caracteristique dynamique Ce(N) vs statique + courants efficaces ---
    N_grid = linspace(0, p.Ns - 0.1, 600);
    s_grid = (p.Ns - N_grid) / p.Ns;
    Zs = p.Rs + 1i * p.ws * p.Ls;
    Xm = p.ws * p.M;
    Ce_stat = zeros(size(s_grid));
    for k = 1:length(s_grid)
        s = s_grid(k);
        Zr = p.Rr / s + 1i * p.ws * p.Lr;
        Is = p.V / (Zs + Xm^2 / Zr);
        Ir = -1i * Xm * Is / Zr;
        Ce_stat(k) = 3 * p.p * (p.Rr / (s * p.ws)) * abs(Ir)^2;
    end
    [Ce_dec, idec] = max(Ce_stat);
    stats.Ce_decrochage_stat = Ce_dec;
    stats.N_decrochage_stat  = N_grid(idec);
    stats.Ce_dem_stat        = Ce_stat(1);

    fig2 = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 22 8.8]);
    tl2 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl2, '\textbf{Caract\''eristique couple-vitesse et appel de courant au d\''emarrage direct}', ...
          'FontSize', 11.5, 'Interpreter', 'latex');

    nexttile;
    plot(res.N_trmin, res.Ce, 'Color', C.navy, 'LineWidth', 1.25); hold on;
    plot(N_grid, Ce_stat, '--', 'Color', C.crimson, 'LineWidth', 1.6);
    plot(N_grid, p.f * (N_grid*2*pi/60) + p.Cp, ':', 'Color', C.emerald, 'LineWidth', 1.4);
    plot(stats.N_decrochage_stat, stats.Ce_decrochage_stat, 's', 'Color', C.crimson, ...
         'MarkerFaceColor', C.crimson, 'MarkerSize', 6);
    text(stats.N_decrochage_stat - 460, stats.Ce_decrochage_stat + 1.8, ...
         sprintf('$C_{\\max}^{\\mathrm{stat}} = %.2f\\,\\mathrm{N\\cdot m}$ ($N=%.0f\\,\\mathrm{tr/min}$)', ...
         stats.Ce_decrochage_stat, stats.N_decrochage_stat), 'Color', C.crimson, 'FontSize', 8.5);
    xlabel('Vitesse de rotation $N$ (tr/min)');
    ylabel('Couple $C_e$ ($\mathrm{N\cdot m}$)');
    title('(a) Caract\''eristique dynamique $C_e(N)$ vs statique');
    legend({'Dynamique (RK4)', 'Statique th\''eorique', 'Pertes $f\Omega + C_p$'}, ...
           'Location', 'northwest', 'FontSize', 8.5);
    xlim([0 1550]); ylim([-4 33]);

    nexttile;
    plot(res.t, res.Is_eff, 'Color', C.navy, 'LineWidth', 1.5); hold on;
    plot(res.t, res.Ir_eff, 'Color', C.crimson, 'LineWidth', 1.5);
    yline(stats.Is_eff_fin, '--', 'Color', C.navy, 'LineWidth', 0.9);
    text(0.25, stats.Is_eff_max + 2.0, sprintf('$I_{s,\\mathrm{eff,max}} = %.2f\\,\\mathrm{A}$ ($%.1f\\times I_{s0}$)', ...
         stats.Is_eff_max, stats.Is_eff_max/stats.Is_eff_fin), 'Color', C.navy, 'FontSize', 8.5);
    text(0.45, stats.Ir_eff_max - 2.5, sprintf('$I_{r,\\mathrm{eff,max}} = %.2f\\,\\mathrm{A}$', ...
         stats.Ir_eff_max), 'Color', C.crimson, 'FontSize', 8.5);
    xlabel('Temps $t$ (s)');
    ylabel('Valeur efficace instantan\''ee (A)');
    title('(b) Enveloppes efficaces $I_{s,\mathrm{eff}}(t)$ et $I_{r,\mathrm{eff}}(t)$');
    legend({'$I_{s,\mathrm{eff}}(t)$ (stator)', '$I_{r,\mathrm{eff}}(t)$ (rotor)'}, ...
           'Location', 'northeast', 'FontSize', 8.5);
    xlim([0 tf]); ylim([0 44]);

    exportgraphics(fig2, fullfile(fig_dir, 'fig_tp1_stator_dynamique.png'), 'Resolution', 300);
    close(fig2);
end
