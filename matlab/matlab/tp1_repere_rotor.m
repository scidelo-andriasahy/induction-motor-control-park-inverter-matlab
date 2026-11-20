function stats = tp1_repere_rotor()
% TP1_REPERE_ROTOR Manipulation I - Partie II :
% Demarrage direct a vide dans le repere lie au rotor (wa = w = p*Omega).

    script_dir = fileparts(mfilename('fullpath'));
    fig_dir = fullfile(script_dir, '..', 'figures');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

    C = set_pub_style();
    p = parametres_mas();
    h = 5e-4;
    tf = 3.0;

    res = simuler_mas_rk4('rotor', p, h, tf);

    stats.Ce_max     = max(res.Ce);
    stats.Ce_fin     = res.Ce(end);
    stats.N_fin      = res.N_trmin(end);
    stats.Is_eff_fin = res.Is_eff(end);
    stats.Ir_eff_fin = res.Ir_eff(end);
    stats.f_gliss    = p.fs * (p.Ns - stats.N_fin) / p.Ns;
    stats.res        = res;

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 22 15]);
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, '\textbf{Manipulation I -- D\''emarrage direct \`a vide (Rep\`ere li\''e au rotor, $\omega_a = \omega = p\Omega$)}', ...
          'FontSize', 12, 'Interpreter', 'latex');

    % (a) Couple electromagnetique
    nexttile;
    plot(res.t, res.Ce, 'Color', C.navy, 'LineWidth', 1.3); hold on;
    yline(stats.Ce_fin, '--', 'Color', C.crimson, 'LineWidth', 1.0);
    xlabel('Temps $t$ (s)');
    ylabel('Couple $C_e$ ($\mathrm{N\cdot m}$)');
    title('(a) Couple \''electromagn\''etique $C_e(t)$');
    xlim([0 tf]); ylim([-5 34]);

    % (b) Vitesse mecanique
    nexttile;
    plot(res.t, res.N_trmin, 'Color', C.navy, 'LineWidth', 1.5); hold on;
    yline(p.Ns, '--', 'Color', C.slate, 'LineWidth', 1.0);
    xlabel('Temps $t$ (s)');
    ylabel('Vitesse $N$ (tr/min)');
    title('(b) Vitesse de rotation $N(t)$');
    legend({'$N(t)$', '$N_s = 1500\,\mathrm{tr/min}$'}, 'Location', 'southeast', 'FontSize', 8.5);
    xlim([0 tf]); ylim([0 1650]);

    % (c) Courants statoriques dans le repere rotor (oscillation lente a f_r = s*f_s)
    nexttile;
    plot(res.t, res.isd, 'Color', C.navy, 'LineWidth', 1.25); hold on;
    plot(res.t, res.isq, 'Color', C.crimson, 'LineWidth', 1.15);
    xlabel('Temps $t$ (s)');
    ylabel('Courants statoriques (A)');
    title(sprintf('(c) Composantes $i_{sd}^{(r)}, i_{sq}^{(r)}$ (fr\\''equence $s f_s \\approx %.2f\\,\\mathrm{Hz}$)', stats.f_gliss));
    legend({'$i_{sd}^{(r)}(t)$', '$i_{sq}^{(r)}(t)$'}, 'Location', 'northeast', 'FontSize', 8.5);
    xlim([0 tf]); ylim([-25 25]);

    % (d) Courants rotoriques dans le repere rotor (proportionnels aux courants reels)
    nexttile;
    plot(res.t, res.ird, 'Color', C.emerald, 'LineWidth', 1.25); hold on;
    plot(res.t, res.irq, 'Color', C.amber, 'LineWidth', 1.15);
    xlabel('Temps $t$ (s)');
    ylabel('Courants rotoriques (A)');
    title('(d) Composantes $i_{rd}^{(r)}, i_{rq}^{(r)}$ (li\''ees directement au rotor)');
    legend({'$i_{rd}^{(r)}(t)$', '$i_{rq}^{(r)}(t)$'}, 'Location', 'northeast', 'FontSize', 8.5);
    xlim([0 tf]); ylim([-58 58]);

    exportgraphics(fig, fullfile(fig_dir, 'fig_tp1_rotor_general.png'), 'Resolution', 300);
    close(fig);
end
