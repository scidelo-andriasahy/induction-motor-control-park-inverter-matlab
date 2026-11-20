function stats = tp1_repere_synchronisme()
% TP1_REPERE_SYNCHRONISME Manipulation I - Partie III (facultatif) :
% Demarrage direct a vide dans le repere lie au synchronisme (wa = ws).

    script_dir = fileparts(mfilename('fullpath'));
    fig_dir = fullfile(script_dir, '..', 'figures');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

    C = set_pub_style();
    p = parametres_mas();
    h = 5e-4;
    tf = 3.0;

    res = simuler_mas_rk4('synchronisme', p, h, tf);

    stats.Ce_max     = max(res.Ce);
    stats.Ce_fin     = res.Ce(end);
    stats.N_fin      = res.N_trmin(end);
    stats.isd_fin    = res.isd(end);
    stats.isq_fin    = res.isq(end);
    stats.ird_fin    = res.ird(end);
    stats.irq_fin    = res.irq(end);
    stats.Is_eff_fin = res.Is_eff(end);
    stats.Ir_eff_fin = res.Ir_eff(end);
    stats.res        = res;

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 22 15]);
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, '\textbf{Manipulation I -- D\''emarrage direct \`a vide (Rep\`ere li\''e au synchronisme, $\omega_a = \omega_s$)}', ...
          'FontSize', 12, 'Interpreter', 'latex');

    % (a) Couple
    nexttile;
    plot(res.t, res.Ce, 'Color', C.navy, 'LineWidth', 1.3); hold on;
    yline(stats.Ce_fin, '--', 'Color', C.crimson, 'LineWidth', 1.0);
    xlabel('Temps $t$ (s)');
    ylabel('Couple $C_e$ ($\mathrm{N\cdot m}$)');
    title('(a) Couple \''electromagn\''etique $C_e(t)$');
    xlim([0 tf]); ylim([-5 34]);

    % (b) Vitesse
    nexttile;
    plot(res.t, res.N_trmin, 'Color', C.navy, 'LineWidth', 1.5); hold on;
    yline(p.Ns, '--', 'Color', C.slate, 'LineWidth', 1.0);
    xlabel('Temps $t$ (s)');
    ylabel('Vitesse $N$ (tr/min)');
    title('(b) Vitesse de rotation $N(t)$');
    legend({'$N(t)$', '$N_s = 1500\,\mathrm{tr/min}$'}, 'Location', 'southeast', 'FontSize', 8.5);
    xlim([0 tf]); ylim([0 1650]);

    % (c) Courants statoriques dans le repere synchrone (grandeurs continues en regime etabli!)
    nexttile;
    plot(res.t, res.isd, 'Color', C.navy, 'LineWidth', 1.35); hold on;
    plot(res.t, res.isq, 'Color', C.crimson, 'LineWidth', 1.35);
    text(1.75, stats.isd_fin + 3.0, sprintf('$i_{sd,\\infty}^{(s)} = %.2f\\,\\mathrm{A}$', stats.isd_fin), ...
         'Color', C.navy, 'FontSize', 8.5);
    text(1.75, stats.isq_fin - 3.5, sprintf('$i_{sq,\\infty}^{(s)} = %.2f\\,\\mathrm{A}$', stats.isq_fin), ...
         'Color', C.crimson, 'FontSize', 8.5);
    xlabel('Temps $t$ (s)');
    ylabel('Courants statoriques (A)');
    title('(c) Composantes $i_{sd}^{(s)}, i_{sq}^{(s)}$ (constantes en r\''egime \''etabli)');
    legend({'$i_{sd}^{(s)}(t)$', '$i_{sq}^{(s)}(t)$'}, 'Location', 'northeast', 'FontSize', 8.5);
    xlim([0 tf]); ylim([-25 25]);

    % (d) Courants rotoriques dans le repere synchrone (grandeurs continues en regime etabli!)
    nexttile;
    plot(res.t, res.ird, 'Color', C.emerald, 'LineWidth', 1.35); hold on;
    plot(res.t, res.irq, 'Color', C.amber, 'LineWidth', 1.35);
    text(1.75, stats.ird_fin - 6.5, sprintf('$i_{rd,\\infty}^{(s)} = %.2f\\,\\mathrm{A}$', stats.ird_fin), ...
         'Color', C.emerald, 'FontSize', 8.5);
    text(1.75, stats.irq_fin + 6.5, sprintf('$i_{rq,\\infty}^{(s)} = %.2f\\,\\mathrm{A}$', stats.irq_fin), ...
         'Color', C.amber, 'FontSize', 8.5);
    xlabel('Temps $t$ (s)');
    ylabel('Courants rotoriques (A)');
    title('(d) Composantes $i_{rd}^{(s)}, i_{rq}^{(s)}$ (constantes en r\''egime \''etabli)');
    legend({'$i_{rd}^{(s)}(t)$', '$i_{rq}^{(s)}(t)$'}, 'Location', 'northeast', 'FontSize', 8.5);
    xlim([0 tf]); ylim([-58 58]);

    exportgraphics(fig, fullfile(fig_dir, 'fig_tp1_sync_general.png'), 'Resolution', 300);
    close(fig);
end
