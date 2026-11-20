function stats = tp2_onduleur()
% TP2_ONDULEUR Manipulation II - Travail facultatif :
% Moteur asynchrone alimente par un onduleur de tension triphase a deux niveaux :
%   1) Comparaison Direct sinusoidal vs Onduleur pleine onde 180 deg vs MLI SPWM (50 Hz)
%   2) Variation de vitesse par commande scalaire V/f = cte sous MLI (25 -> 37.5 -> 50 Hz)

    script_dir = fileparts(mfilename('fullpath'));
    fig_dir = fullfile(script_dir, '..', 'figures');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

    C = set_pub_style();
    p = parametres_mas();

    % 1. Comparaison a 50 Hz (meme tension fondamentale V1_eff = 220 V)
    h1 = 1e-5;   % 10 us -> 50 echantillons par periode porteuse (fp = 2000 Hz)
    tf1 = 2.5;

    r_sin = simuler_mas_rk4('stator', p, 5e-4, tf1);

    p_po.f_func = @(t) 50;
    p_po.Edc    = (pi / sqrt(2)) * p.V; % 488.72 V -> V1_eff = 220 V
    p_po.use_vf = false;
    r_po = simuler_mas_onduleur('pleine_onde', p_po, p, 2e-5, tf1);

    p_mli.f_func = @(t) 50;
    p_mli.fp     = 2000;                % m = fp / fs = 40
    p_mli.Edc    = 2 * sqrt(2) * p.V;   % 622.25 V -> V1_eff = 220 V pour r = 1
    r_mli = simuler_mas_onduleur('mli', p_mli, p, h1, tf1);

    % Extraction des statistiques en regime etabli (t >= 2.0 s)
    t_po  = r_po.t >= 2.0;
    t_mli = r_mli.t >= 2.0;

    stats.po.N_fin      = mean(r_po.N_trmin(t_po));
    stats.po.Ce_max     = max(r_po.Ce);
    stats.po.dCe_pp     = max(r_po.Ce(t_po)) - min(r_po.Ce(t_po));
    stats.po.std_Ce     = std(r_po.Ce(t_po));
    stats.po.isa_max    = max(abs(r_po.isa));
    stats.po.isa_rms    = rms(r_po.isa(t_po));
    [f_po, A_po, stats.po.thd_i, ~] = calc_fft_thd(r_po.t, r_po.isa, 2.0, p.fs, 20);

    stats.mli.N_fin     = mean(r_mli.N_trmin(t_mli));
    stats.mli.Ce_max    = max(r_mli.Ce);
    stats.mli.dCe_pp    = max(r_mli.Ce(t_mli)) - min(r_mli.Ce(t_mli));
    stats.mli.std_Ce    = std(r_mli.Ce(t_mli));
    stats.mli.isa_max   = max(abs(r_mli.isa));
    stats.mli.isa_rms   = rms(r_mli.isa(t_mli));
    [f_mli, A_mli, stats.mli.thd_i, ~] = calc_fft_thd(r_mli.t, r_mli.isa, 2.0, p.fs, 20);

    % --- Figure 1 : Comparaison Direct vs Pleine Onde 180 deg vs MLI SPWM ---
    fig1 = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 23.5 15.5]);
    tl1 = tiledlayout(fig1, 2, 2, 'TileSpacing', 'normal', 'Padding', 'compact');
    title(tl1, '\textbf{Manipulation II (Facultatif) -- Onduleur : Pleine onde $180^\circ$ vs MLI ($f_p = 2\,\mathrm{kHz}$)}', ...
          'FontSize', 11.5, 'Interpreter', 'latex');

    iz_po  = r_po.t  >= 2.20 & r_po.t  <= 2.22;
    iz_mli = r_mli.t >= 2.20 & r_mli.t <= 2.22;

    % (a) Couple electromagnetique
    nexttile(tl1); hold on;
    plot(r_sin.t,            r_sin.Ce,            'Color', C.black,   'LineWidth', 1.3);
    plot(r_po.t(1:20:end),   r_po.Ce(1:20:end),   'Color', C.crimson, 'LineWidth', 1.0);
    plot(r_mli.t(1:40:end),  r_mli.Ce(1:40:end),  'Color', C.navy,    'LineWidth', 1.1);
    xlabel('Temps $t$ (s)'); ylabel('Couple $C_e$ ($\mathrm{N\cdot m}$)');
    title('(a) Couple \''electromagn\''etique $C_e(t)$', 'FontSize', 9.8);
    legend({'Sinuso\"idal direct', 'Pleine onde $180^\circ$', 'MLI ($2\,\mathrm{kHz}$)'}, ...
           'Location', 'northeast', 'FontSize', 7.8);
    xlim([0 tf1]); ylim([-5 36]);

    % (b) Tension simple vsa(t) en regime etabli (1 periode = 20 ms)
    nexttile(tl1); hold on;
    plot(r_mli.t(iz_mli)*1e3, r_mli.vsa(iz_mli), 'Color', [0.043 0.235 0.365 0.55], 'LineWidth', 0.8);
    plot(r_po.t(iz_po)*1e3,   r_po.vsa(iz_po),   'Color', C.crimson, 'LineWidth', 1.6);
    xlabel('Temps $t$ (ms)'); ylabel('Tension $v_{sa}$ (V)');
    title('(b) Tension simple $v_{sa}(t)$ (6 marches vs MLI)', 'FontSize', 9.8);
    legend({'MLI ($E_{dc}=622\,\mathrm{V}$)', 'Pleine onde ($E_{dc}=489\,\mathrm{V}$)'}, ...
           'Location', 'northeast', 'FontSize', 7.8);
    xlim([2200 2220]); ylim([-460 460]);

    % (c) Courant de phase isa(t) en regime etabli (2 periodes = 40 ms)
    nexttile(tl1); hold on;
    iz2_po  = r_po.t  >= 2.20 & r_po.t  <= 2.24;
    iz2_mli = r_mli.t >= 2.20 & r_mli.t <= 2.24;
    plot(r_po.t(iz2_po)*1e3,   r_po.isa(iz2_po),   'Color', C.crimson, 'LineWidth', 1.4);
    plot(r_mli.t(iz2_mli)*1e3, r_mli.isa(iz2_mli), 'Color', C.navy,    'LineWidth', 1.25);
    xlabel('Temps $t$ (ms)'); ylabel('Courant $i_{sa}$ (A)');
    title('(c) Courant $i_{sa}(t)$ en r\''egime \''etabli (filtrage)', 'FontSize', 9.8);
    legend({'Pleine onde $180^\circ$ (d\''eform\''e)', 'MLI $2\,\mathrm{kHz}$ (quasi-sinuso\"idal)'}, ...
           'Location', 'northeast', 'FontSize', 7.8);
    xlim([2200 2240]); ylim([-6.5 6.5]);

    % (d) Spectre harmonique FFT compare (0 .. 2400 Hz)
    nexttile(tl1); hold on;
    m_po  = f_po  <= 2400;
    m_mli = f_mli <= 2400;
    plot(f_po(m_po),   A_po(m_po),   'Color', C.crimson, 'LineWidth', 1.35);
    plot(f_mli(m_mli), A_mli(m_mli), 'Color', C.navy,    'LineWidth', 1.15);
    xlabel('Fr\''equence $f$ (Hz)'); ylabel('Amplitude $|I_{sa}(f)|$ (A)');
    title('(d) Spectre harmonique de $i_{sa}(t)$ (rejet vers $f_p$)', 'FontSize', 9.8);
    legend({sprintf('Pleine onde (THD $= %.1f\\%%$, rangs $5,7,11,13$)', stats.po.thd_i), ...
            sprintf('MLI $2\\,\\mathrm{kHz}$ (THD $= %.1f\\%%$, rejet $f_p$)', stats.mli.thd_i)}, ...
           'Location', 'northeast', 'FontSize', 7.8);
    xlim([0 2400]); ylim([0 5.5]);

    % Inset zoom ondulation couple place APRES tous les nexttile
    ax_c = axes(fig1, 'Position', [0.27 0.61 0.175 0.115]);
    box on; hold(ax_c, 'on');
    plot(ax_c, r_po.t(iz_po)*1e3,   r_po.Ce(iz_po),   'Color', C.crimson, 'LineWidth', 1.1);
    plot(ax_c, r_mli.t(iz_mli)*1e3, r_mli.Ce(iz_mli), 'Color', C.navy,    'LineWidth', 0.9);
    xlim(ax_c, [2200 2220]); ylim(ax_c, [0.4 2.8]);
    title(ax_c, 'Ondulation $C_e$ (20 ms)', 'FontSize', 7.3);
    set(ax_c, 'FontSize', 6.8);

    exportgraphics(fig1, fullfile(fig_dir, 'fig_tp2_onduleur_comp.png'), 'Resolution', 300);
    close(fig1);

    % 2. Variation de vitesse par commande scalaire V/f = cte sous MLI
    % Profil : rampe douce 5->25 Hz (0..0.5s), palier 25 Hz (0.5..2.0s),
    %          rampe 25->37.5 Hz (2.0..2.5s), palier 37.5 Hz (2.5..3.8s),
    %          rampe 37.5->50 Hz (3.8..4.3s), palier 50 Hz (4.3..5.5s)
    f_prof = @(t) (5 + 40*t).*(t < 0.5) + ...
                  25.*(t >= 0.5 & t < 2.0) + ...
                  (25 + 25*(t - 2.0)).*(t >= 2.0 & t < 2.5) + ...
                  37.5.*(t >= 2.5 & t < 3.8) + ...
                  (37.5 + 25*(t - 3.8)).*(t >= 3.8 & t < 4.3) + ...
                  50.*(t >= 4.3);

    p_vf.f_func = f_prof;
    p_vf.fp     = 2000;
    p_vf.Edc    = 2 * sqrt(2) * p.V;
    h2 = 2e-5;
    tf2 = 5.5;

    r_vf = simuler_mas_onduleur('mli', p_vf, p, h2, tf2);

    stats.vf.N_25   = mean(r_vf.N_trmin(r_vf.t >= 1.7 & r_vf.t <= 1.95));
    stats.vf.N_375  = mean(r_vf.N_trmin(r_vf.t >= 3.5 & r_vf.t <= 3.75));
    stats.vf.N_50   = mean(r_vf.N_trmin(r_vf.t >= 5.2 & r_vf.t <= 5.5));
    stats.vf.Is_max = max(r_vf.Is_eff);
    stats.vf.Ce_max = max(r_vf.Ce);

    fig2 = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 23.5 15.0]);
    tl2 = tiledlayout(fig2, 2, 2, 'TileSpacing', 'normal', 'Padding', 'compact');
    title(tl2, '\textbf{Variation de vitesse par onduleur MLI \`a loi scalaire $V/f = \mathrm{cte}$ ($25 \to 37{,}5 \to 50\,\mathrm{Hz}$)}', ...
          'FontSize', 11.5, 'Interpreter', 'latex');

    sp = 20;
    % (a) Vitesse mecanique vs synchronisme instantane
    nexttile(tl2); hold on;
    plot(r_vf.t(1:sp:end), 60*r_vf.f_ref(1:sp:end)/p.p, '--', 'Color', C.slate, 'LineWidth', 1.3);
    plot(r_vf.t(1:sp:end), r_vf.N_trmin(1:sp:end),      '-',  'Color', C.navy,  'LineWidth', 1.5);
    xlabel('Temps $t$ (s)'); ylabel('Vitesse $N$ (tr/min)');
    title('(a) Poursuite de vitesse $N(t)$ et $N_s(t)$', 'FontSize', 9.8);
    legend({'$N_s(t) = 60 f_{\mathrm{ref}}/p$', 'Vitesse $N(t)$'}, 'Location', 'southeast', 'FontSize', 8.5);
    xlim([0 tf2]); ylim([0 1650]);

    % (b) Couple electromagnetique
    nexttile(tl2);
    plot(r_vf.t(1:sp:end), r_vf.Ce(1:sp:end), 'Color', C.emerald, 'LineWidth', 1.15);
    xlabel('Temps $t$ (s)'); ylabel('Couple $C_e$ ($\mathrm{N\cdot m}$)');
    title('(b) Couple $C_e(t)$ (d\''emarrage ma\^itris\''e)', 'FontSize', 9.8);
    xlim([0 tf2]); ylim([-4 19]);

    % (c) Courant efficace statorique (flux constant -> appel de courant maitrise!)
    nexttile(tl2);
    plot(r_vf.t(1:sp:end), r_vf.Is_eff(1:sp:end), 'Color', C.crimson, 'LineWidth', 1.3);
    xlabel('Temps $t$ (s)'); ylabel('$I_{s,\mathrm{eff}}$ (A)');
    title('(c) Courant statorique $I_{s,\mathrm{eff}}(t)$ (flux constant)', 'FontSize', 9.8);
    xlim([0 tf2]); ylim([0 9]);

    % (d) Courant rotorique ira(t) et frequence de commande f_ref(t)
    nexttile(tl2);
    yyaxis left;
    plot(r_vf.t(1:sp:end), r_vf.Ir_eff(1:sp:end), 'Color', C.amber, 'LineWidth', 1.3);
    ylabel('$I_{r,\mathrm{eff}}$ (A)', 'Color', C.amber);
    ylim([0 18]);
    ax = gca; ax.YAxis(1).Color = C.amber;

    yyaxis right;
    plot(r_vf.t(1:sp:end), r_vf.f_ref(1:sp:end), '--', 'Color', C.navy, 'LineWidth', 1.5);
    ylabel('Fr\''equence $f_{\mathrm{ref}}$ (Hz)', 'Color', C.navy);
    ylim([0 60]);
    ax.YAxis(2).Color = C.navy;

    xlabel('Temps $t$ (s)');
    title('(d) Courant rotorique $I_{r,\mathrm{eff}}(t)$ et consigne $f_{\mathrm{ref}}(t)$', 'FontSize', 9.8);
    xlim([0 tf2]);

    exportgraphics(fig2, fullfile(fig_dir, 'fig_tp2_onduleur_vf.png'), 'Resolution', 300);
    close(fig2);
end
