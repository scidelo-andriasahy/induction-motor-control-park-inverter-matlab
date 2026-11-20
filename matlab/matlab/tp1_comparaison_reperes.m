function stats = tp1_comparaison_reperes(r_st, r_ro, r_sy)
% TP1_COMPARAISON_REPERES Manipulation I - Partie IV :
% Comparaison quantitative et qualitative des trois referentiels de Park
% (stator, rotor, synchronisme).

    script_dir = fileparts(mfilename('fullpath'));
    fig_dir = fullfile(script_dir, '..', 'figures');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

    C = set_pub_style();
    p = parametres_mas();
    if nargin < 3
        r_st = simuler_mas_rk4('stator',       p, 5e-4, 3.0);
        r_ro = simuler_mas_rk4('rotor',        p, 5e-4, 3.0);
        r_sy = simuler_mas_rk4('synchronisme', p, 5e-4, 3.0);
    end

    stats.err_Ce_st_ro    = max(abs(r_st.Ce - r_ro.Ce));
    stats.err_Ce_st_sy    = max(abs(r_st.Ce - r_sy.Ce));
    stats.err_N_st_ro     = max(abs(r_st.N_trmin - r_ro.N_trmin));
    stats.err_N_st_sy     = max(abs(r_st.N_trmin - r_sy.N_trmin));
    stats.err_isa_st_ro   = max(abs(r_st.isa - r_ro.isa));
    stats.err_isa_st_sy   = max(abs(r_st.isa - r_sy.isa));

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 23 15.2]);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, '\textbf{Manipulation I -- Comparaison des trois r\''ef\''erentiels de Park (Stator, Rotor, Synchronisme)}', ...
          'FontSize', 12, 'Interpreter', 'latex');

    % (a) Superposition des couples electromagnetiques
    nexttile(tl);
    plot(r_st.t, r_st.Ce, '-',  'Color', C.navy,    'LineWidth', 2.0); hold on;
    plot(r_ro.t, r_ro.Ce, '--', 'Color', C.crimson, 'LineWidth', 1.4);
    plot(r_sy.t, r_sy.Ce, ':',  'Color', C.emerald, 'LineWidth', 1.8);
    xlabel('Temps $t$ (s)');
    ylabel('Couple $C_e$ ($\mathrm{N\cdot m}$)');
    title('(a) Invariance du couple $C_e(t)$ (superposition)', 'FontSize', 10);
    legend({'Rep\`ere Stator', 'Rep\`ere Rotor', 'Rep\`ere Synchronisme'}, ...
           'Location', 'northeast', 'FontSize', 8);
    xlim([0 3.0]); ylim([-4 33]);

    % (b) Superposition des vitesses
    nexttile(tl);
    plot(r_st.t, r_st.N_trmin, '-',  'Color', C.navy,    'LineWidth', 2.0); hold on;
    plot(r_ro.t, r_ro.N_trmin, '--', 'Color', C.crimson, 'LineWidth', 1.4);
    plot(r_sy.t, r_sy.N_trmin, ':',  'Color', C.emerald, 'LineWidth', 1.8);
    xlabel('Temps $t$ (s)');
    ylabel('Vitesse $N$ (tr/min)');
    title('(b) Invariance de la vitesse $N(t)$ entre les 3 rep\`eres', 'FontSize', 10);
    legend({'Rep\`ere Stator', 'Rep\`ere Rotor', 'Rep\`ere Synchronisme'}, ...
           'Location', 'southeast', 'FontSize', 8);
    xlim([0 3.0]); ylim([0 1650]);

    % (c) Comparaison de i_sd(t) en regime etabli (signature frequentielle!)
    nexttile(tl);
    idx_z = r_st.t >= 2.0 & r_st.t <= 2.6;
    plot(r_st.t(idx_z), r_st.isd(idx_z), 'Color', [0.043 0.235 0.365 0.6], 'LineWidth', 0.9); hold on;
    plot(r_ro.t(idx_z), r_ro.isd(idx_z), 'Color', C.crimson, 'LineWidth', 1.6);
    plot(r_sy.t(idx_z), r_sy.isd(idx_z), 'Color', C.emerald, 'LineWidth', 2.0);
    xlabel('Temps $t$ (s)');
    ylabel('Composante $i_{sd}$ (A)');
    title('(c) Signature fr\''equentielle de $i_{sd}$ en r\''egime \''etabli', 'FontSize', 10);
    legend({'Stator ($50\,\mathrm{Hz}$)', 'Rotor ($s f_s \approx 0{,}6\,\mathrm{Hz}$)', ...
            'Synchronisme (continu DC)'}, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 7.8);
    xlim([2.0 2.6]); ylim([-7 7]);

    % (d) Invariance du courant reel de phase statorique i_sa(t) reconstruit
    nexttile(tl);
    idx_ph = r_st.t >= 2.40 & r_st.t <= 2.48;
    plot(r_st.t(idx_ph), r_st.isa(idx_ph), '-',  'Color', C.navy,    'LineWidth', 2.0); hold on;
    plot(r_ro.t(idx_ph), r_ro.isa(idx_ph), '--', 'Color', C.crimson, 'LineWidth', 1.4);
    plot(r_sy.t(idx_ph), r_sy.isa(idx_ph), ':',  'Color', C.emerald, 'LineWidth', 1.8);
    xlabel('Temps $t$ (s)');
    ylabel('Courant de phase $i_{sa}$ (A)');
    title('(d) Courant de phase $i_{sa}(t)$ reconstruit par $P^{-1}(\theta_s)$', 'FontSize', 10);
    legend({'Depuis Stator', 'Depuis Rotor', 'Depuis Synchronisme'}, ...
           'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 7.8);
    xlim([2.40 2.48]); ylim([-6 6]);

    % Inset ecart numerique place APRES tous les nexttile
    ax_e = axes(fig, 'Position', [0.73 0.72 0.185 0.115]);
    box on;
    semilogy(ax_e, r_st.t, abs(r_st.Ce - r_ro.Ce) + 1e-12, 'Color', C.crimson, 'LineWidth', 1.0); hold(ax_e, 'on');
    semilogy(ax_e, r_st.t, abs(r_st.Ce - r_sy.Ce) + 1e-12, 'Color', C.emerald, 'LineWidth', 1.0);
    xlim(ax_e, [0 3.0]); ylim(ax_e, [1e-7 1e-1]);
    title(ax_e, '$|\Delta C_e|$ ($\mathrm{N\cdot m}$)', 'FontSize', 7.5);
    set(ax_e, 'FontSize', 6.8);

    exportgraphics(fig, fullfile(fig_dir, 'fig_tp1_comparaison_reperes.png'), 'Resolution', 300);
    close(fig);
end
