function stats = tp2_gradateur_angle_phase()
% TP2_GRADATEUR_ANGLE_PHASE Manipulation II - Partie I :
% Moteur asynchrone alimente par gradateur triphase a thyristors
% en commande par decoupage de phase (angle d'amorcage alpha fixe et rampe progressive).

    script_dir = fileparts(mfilename('fullpath'));
    fig_dir = fullfile(script_dir, '..', 'figures');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

    C = set_pub_style();
    p = parametres_mas();
    h = 5e-5;
    tf = 5.0;

    alphas = [0, 80, 90, 95, 100];
    cols = {C.black, C.navy, C.emerald, C.amber, C.crimson};
    res_alp = cell(1, length(alphas));

    for k = 1:length(alphas)
        r = simuler_mas_gradateur('angle_phase', struct('alpha_deg', alphas(k)), p, h, tf);
        res_alp{k} = r;
        tail = r.t >= 4.0;
        stats.alphas(k)      = alphas(k);
        stats.N_fin(k)       = mean(r.N_trmin(tail));
        stats.isa_max(k)     = max(abs(r.isa));
        stats.Is_eff_max(k)  = max(r.Is_eff);
        stats.Is_eff_fin(k)  = mean(r.Is_eff(tail));
        stats.Ir_eff_fin(k)  = mean(r.Ir_eff(tail));
        stats.Ce_max(k)      = max(r.Ce);
        stats.Vsa_rms_fin(k) = rms(r.vsa(tail));
        idx95 = find(r.N_trmin >= 0.95 * stats.N_fin(k), 1);
        stats.t95(k) = r.t(idx95);
        [~, ~, thd_i] = calc_fft_thd(r.t, r.isa, 4.0, p.fs, 20);
        stats.thd_i(k) = thd_i;
    end

    % Demarrage progressif industriel (rampe d'angle alpha : 110 deg -> 70 deg sur 2.0 s)
    r_rampe = simuler_mas_gradateur('rampe_alpha', ...
        struct('alpha_init_deg', 110, 'alpha_fin_deg', 70, 't_rampe', 2.0), p, h, tf);
    tail = r_rampe.t >= 4.0;
    stats.rampe.N_fin      = mean(r_rampe.N_trmin(tail));
    stats.rampe.isa_max    = max(abs(r_rampe.isa));
    stats.rampe.Is_eff_max = max(r_rampe.Is_eff);
    stats.rampe.Ce_max     = max(r_rampe.Ce);
    idx95r = find(r_rampe.N_trmin >= 0.95 * stats.rampe.N_fin, 1);
    stats.rampe.t95        = r_rampe.t(idx95r);
    [~, ~, stats.rampe.thd_i] = calc_fft_thd(r_rampe.t, r_rampe.isa, 4.0, p.fs, 20);

    % Sauvegarde des structures utiles pour la comparaison globale
    stats.r_direct = res_alp{1};
    stats.r_90     = res_alp{3};
    stats.r_95     = res_alp{4};
    stats.r_rampe  = r_rampe;

    % --- Figure 1 : Grandeurs globales pour differents angles alpha + rampe ---
    fig1 = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 22 15.2]);
    tl1 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl1, '\textbf{Manipulation II -- Gradateur en commande par angle de phase ($\alpha$ fixe et rampe progressive)}', ...
          'FontSize', 11.5, 'Interpreter', 'latex');

    lbls = {'Direct ($\alpha = 0^\circ$)', '$\alpha = 80^\circ$', '$\alpha = 90^\circ$', ...
            '$\alpha = 95^\circ$', '$\alpha = 100^\circ$', 'Rampe $110^\circ \to 70^\circ$'};

    % Decimation legere pour tracage rapide et propre
    step_p = 10;

    % (a) Couple electromagnetique
    nexttile; hold on;
    for k = 1:length(alphas)
        plot(res_alp{k}.t(1:step_p:end), res_alp{k}.Ce(1:step_p:end), 'Color', cols{k}, 'LineWidth', 1.2);
    end
    plot(r_rampe.t(1:step_p:end), r_rampe.Ce(1:step_p:end), '--', 'Color', C.purple, 'LineWidth', 1.5);
    xlabel('Temps $t$ (s)'); ylabel('Couple $C_e$ ($\mathrm{N\cdot m}$)');
    title('(a) Couple \''electromagn\''etique $C_e(t)$');
    legend(lbls, 'Location', 'northeast', 'FontSize', 7.5);
    xlim([0 tf]); ylim([-3 33]);

    % (b) Vitesse mecanique
    nexttile; hold on;
    for k = 1:length(alphas)
        plot(res_alp{k}.t(1:step_p:end), res_alp{k}.N_trmin(1:step_p:end), 'Color', cols{k}, 'LineWidth', 1.35);
    end
    plot(r_rampe.t(1:step_p:end), r_rampe.N_trmin(1:step_p:end), '--', 'Color', C.purple, 'LineWidth', 1.6);
    xlabel('Temps $t$ (s)'); ylabel('Vitesse $N$ (tr/min)');
    title('(b) Vitesse de rotation $N(t)$');
    legend(lbls, 'Location', 'southeast', 'FontSize', 7.5);
    xlim([0 tf]); ylim([0 1620]);

    % (c) Courant efficace statorique
    nexttile; hold on;
    for k = 1:length(alphas)
        plot(res_alp{k}.t(1:step_p:end), res_alp{k}.Is_eff(1:step_p:end), 'Color', cols{k}, 'LineWidth', 1.35);
    end
    plot(r_rampe.t(1:step_p:end), r_rampe.Is_eff(1:step_p:end), '--', 'Color', C.purple, 'LineWidth', 1.6);
    xlabel('Temps $t$ (s)'); ylabel('$I_{s,\mathrm{eff}}$ (A)');
    title('(c) Courant efficace statorique $I_{s,\mathrm{eff}}(t)$');
    legend(lbls, 'Location', 'northeast', 'FontSize', 7.5);
    xlim([0 tf]); ylim([0 16]);

    % (d) Courant efficace rotorique
    nexttile; hold on;
    for k = 1:length(alphas)
        plot(res_alp{k}.t(1:step_p:end), res_alp{k}.Ir_eff(1:step_p:end), 'Color', cols{k}, 'LineWidth', 1.35);
    end
    plot(r_rampe.t(1:step_p:end), r_rampe.Ir_eff(1:step_p:end), '--', 'Color', C.purple, 'LineWidth', 1.6);
    xlabel('Temps $t$ (s)'); ylabel('$I_{r,\mathrm{eff}}$ (A)');
    title('(d) Courant efficace rotorique $I_{r,\mathrm{eff}}(t)$');
    legend(lbls, 'Location', 'northeast', 'FontSize', 7.5);
    xlim([0 tf]); ylim([0 44]);

    exportgraphics(fig1, fullfile(fig_dir, 'fig_tp2_angle_phase_global.png'), 'Resolution', 300);
    close(fig1);

    % --- Figure 2 : Zooms formes d'ondes electriques (decoupage thyristors) ---
    r90 = res_alp{3};
    r95 = res_alp{4};
    iz = r95.t >= 4.40 & r95.t <= 4.46; % 3 periodes reseau (60 ms)

    fig2 = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 22 14.5]);
    tl2 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl2, '\textbf{Formes d''ondes d\''etaill\''ees du gradateur \`a angle de phase en r\''egime \''etabli}', ...
          'FontSize', 12, 'Interpreter', 'latex');

    % (a) Tension reseau vs tension appliquee au stator (alpha = 95 deg)
    nexttile;
    plot(r95.t(iz)*1e3, r95.ea_net(iz), '--', 'Color', C.slate, 'LineWidth', 1.1); hold on;
    plot(r95.t(iz)*1e3, r95.vsa(iz),    '-',  'Color', C.navy,  'LineWidth', 1.45);
    xlabel('Temps $t$ (ms)'); ylabel('Tension (V)');
    title('(a) Tension r\''eseau $e_a$ et tension stator $v_{sa}$ ($\alpha = 95^\circ$)');
    legend({'R\''eseau $e_a(t)$', 'Stator $v_{sa}(t)$'}, 'Location', 'northeast', 'FontSize', 8);
    xlim([4400 4460]); ylim([-350 350]);

    % (b) Courant statorique decoupe (alpha = 90 deg vs alpha = 95 deg)
    nexttile;
    plot(r90.t(iz)*1e3, r90.isa(iz), 'Color', C.emerald, 'LineWidth', 1.4); hold on;
    plot(r95.t(iz)*1e3, r95.isa(iz), 'Color', C.crimson, 'LineWidth', 1.4);
    yline(0, ':', 'Color', C.black, 'LineWidth', 0.8);
    xlabel('Temps $t$ (ms)'); ylabel('Courant $i_{sa}$ (A)');
    title('(b) Courant de phase $i_{sa}(t)$ (paliers de blocage $i_{sa} = 0$)');
    legend({'$\alpha = 90^\circ$', '$\alpha = 95^\circ$'}, 'Location', 'northeast', 'FontSize', 8.5);
    xlim([4400 4460]); ylim([-7 7]);

    % (c) Courant rotorique i_ra(t) pour alpha = 90 deg
    nexttile;
    iz_r = r90.t >= 0 & r90.t <= 5.0;
    plot(r90.t(1:step_p:end), r90.ira(1:step_p:end), 'Color', C.emerald, 'LineWidth', 1.0); hold on;
    plot(r95.t(1:step_p:end), r95.ira(1:step_p:end), 'Color', [0.753 0.224 0.169 0.75], 'LineWidth', 0.9);
    xlabel('Temps $t$ (s)'); ylabel('Courant $i_{ra}$ (A)');
    title('(c) Courant rotorique de phase $i_{ra}(t)$ ($\alpha = 90^\circ$ et $95^\circ$)');
    legend({'$\alpha = 90^\circ$', '$\alpha = 95^\circ$'}, 'Location', 'northeast', 'FontSize', 8.5);
    xlim([0 5.0]); ylim([-30 30]);

    % (d) Nombre de phases conductrices (commutation naturelle 3 <-> 2 <-> 0)
    nexttile;
    stairs(r95.t(iz)*1e3, r95.nb_cond(iz), 'Color', C.amber, 'LineWidth', 1.5); hold on;
    stairs(r90.t(iz)*1e3, r90.nb_cond(iz) + 0.08, '--', 'Color', C.navy, 'LineWidth', 1.2);
    xlabel('Temps $t$ (ms)'); ylabel('Phases conductrices');
    title('(d) S\''equence de conduction naturelle des thyristors ($n_{\mathrm{cond}}$)');
    legend({'$\alpha = 95^\circ$', '$\alpha = 90^\circ$ (+0{,}08)'}, 'Location', 'southeast', 'FontSize', 8);
    xlim([4400 4460]); ylim([-0.3 3.6]); yticks([0 1 2 3]);

    exportgraphics(fig2, fullfile(fig_dir, 'fig_tp2_angle_phase_zoom.png'), 'Resolution', 300);
    close(fig2);
end
