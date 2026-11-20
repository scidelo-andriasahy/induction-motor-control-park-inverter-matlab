function tp2_comparaison_gradateurs(st_ap, st_to)
% TP2_COMPARAISON_GRADATEURS Manipulation II - Partie III :
% Comparaison synthetique entre demarrage direct, gradateur a angle de phase
% (fixe et rampe progressive) et gradateur par train d'ondes pleines.

    script_dir = fileparts(mfilename('fullpath'));
    fig_dir = fullfile(script_dir, '..', 'figures');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

    C = set_pub_style();
    p = parametres_mas();
    if nargin < 2
        st_ap = tp2_gradateur_angle_phase();
        st_to = tp2_gradateur_train_ondes();
    end

    r_dir   = st_ap.r_direct;
    r_a90   = st_ap.r_90;
    r_rampe = st_ap.r_rampe;
    r_to68  = st_to.r_68;

    step_p = 10;

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 22 15.2]);
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, '\textbf{Manipulation II -- Synth\`ese comparative : Direct vs Angle de phase vs Train d''ondes}', ...
          'FontSize', 12, 'Interpreter', 'latex');

    % (a) Vitesse de rotation N(t)
    nexttile; hold on;
    plot(r_dir.t(1:step_p:end),   r_dir.N_trmin(1:step_p:end),   'Color', C.black,   'LineWidth', 1.4);
    plot(r_a90.t(1:step_p:end),   r_a90.N_trmin(1:step_p:end),   'Color', C.emerald, 'LineWidth', 1.4);
    plot(r_rampe.t(1:step_p:end), r_rampe.N_trmin(1:step_p:end), '--', 'Color', C.purple, 'LineWidth', 1.6);
    plot(r_to68.t(1:step_p:end),  r_to68.N_trmin(1:step_p:end),  'Color', C.crimson, 'LineWidth', 1.2);
    xlabel('Temps $t$ (s)'); ylabel('Vitesse $N$ (tr/min)');
    title('(a) \''Evolution de la vitesse $N(t)$');
    legend({'Direct ($\alpha = 0^\circ$)', 'Angle de phase ($\alpha = 90^\circ$)', ...
            'Rampe $\alpha(t)$ ($110^\circ \to 70^\circ$)', 'Train d''ondes ($6/8$)'}, ...
           'Location', 'southeast', 'FontSize', 7.8);
    xlim([0 5.0]); ylim([0 1620]);

    % (b) Couple electromagnetique Ce(t)
    nexttile; hold on;
    plot(r_dir.t(1:step_p:end),   r_dir.Ce(1:step_p:end),   'Color', [0.12 0.12 0.12 0.65], 'LineWidth', 1.1);
    plot(r_to68.t(1:step_p:end),  r_to68.Ce(1:step_p:end),  'Color', [0.753 0.224 0.169 0.75], 'LineWidth', 0.95);
    plot(r_rampe.t(1:step_p:end), r_rampe.Ce(1:step_p:end), 'Color', C.purple,  'LineWidth', 1.4);
    plot(r_a90.t(1:step_p:end),   r_a90.Ce(1:step_p:end),   'Color', C.emerald, 'LineWidth', 1.35);
    xlabel('Temps $t$ (s)'); ylabel('Couple $C_e$ ($\mathrm{N\cdot m}$)');
    title('(b) Comparaison des chocs de couple $C_e(t)$');
    legend({'Direct', 'Train d''ondes ($6/8$)', 'Rampe $\alpha(t)$', 'Angle fixe ($\alpha = 90^\circ$)'}, ...
           'Location', 'northeast', 'FontSize', 7.8);
    xlim([0 4.0]); ylim([-18 33]);

    % (c) Courant d'appel au demarrage (enveloppe Is_eff sur 0..3 s)
    nexttile; hold on;
    plot(r_dir.t(1:step_p:end),   r_dir.Is_eff(1:step_p:end),   'Color', C.black,   'LineWidth', 1.4);
    plot(r_to68.t(1:step_p:end),  r_to68.Is_eff(1:step_p:end),  'Color', C.crimson, 'LineWidth', 1.15);
    plot(r_rampe.t(1:step_p:end), r_rampe.Is_eff(1:step_p:end), '--', 'Color', C.purple, 'LineWidth', 1.6);
    plot(r_a90.t(1:step_p:end),   r_a90.Is_eff(1:step_p:end),   'Color', C.emerald, 'LineWidth', 1.4);
    xlabel('Temps $t$ (s)'); ylabel('$I_{s,\mathrm{eff}}$ (A)');
    title('(c) Ma\^itrise de l''appel de courant $I_{s,\mathrm{eff}}(t)$');
    legend({'Direct', 'Train d''ondes ($6/8$)', 'Rampe $\alpha(t)$', 'Angle fixe ($\alpha = 90^\circ$)'}, ...
           'Location', 'northeast', 'FontSize', 7.8);
    xlim([0 4.0]); ylim([0 18]);

    % (d) Spectre harmonique FFT du courant isa(t) en regime etabli
    nexttile; hold on;
    [f_a95, A_a95, thd_a95, ~] = calc_fft_thd(st_ap.r_95.t, st_ap.r_95.isa, 4.0,  p.fs, 20);
    [f_to,  A_to,  thd_to,  ~] = calc_fft_thd(r_to68.t,     r_to68.isa,     3.68, p.fs, 64);
    [f_a90, A_a90, thd_a90, ~] = calc_fft_thd(r_a90.t,      r_a90.isa,      4.0,  p.fs, 20);

    mask_f95 = f_a95 <= 550;
    mask_fto = f_to  <= 550;
    plot(f_to(mask_fto),  A_to(mask_fto),  'Color', C.crimson, 'LineWidth', 1.25);
    plot(f_a95(mask_f95), A_a95(mask_f95), 'Color', C.amber,   'LineWidth', 1.35);
    plot(f_a90(mask_f95), A_a90(mask_f95), 'Color', C.emerald, 'LineWidth', 1.35);
    xlabel('Fr\''equence $f$ (Hz)'); ylabel('Amplitude $|I_{sa}(f)|$ (A)');
    title('(d) Spectre harmonique de $i_{sa}(t)$ en r\''egime \''etabli');
    legend({sprintf('Train d''ondes $6/8$ (THD $= %.1f\\%%$)', thd_to), ...
            sprintf('Angle $\\alpha=95^\\circ$ (THD $= %.1f\\%%$)', thd_a95), ...
            sprintf('Angle $\\alpha=90^\\circ$ (THD $= %.1f\\%%$)', thd_a90)}, ...
           'Location', 'northeast', 'FontSize', 7.8);
    xlim([0 550]); ylim([0 7]);

    exportgraphics(fig, fullfile(fig_dir, 'fig_tp2_comparaison_gradateurs.png'), 'Resolution', 300);
    close(fig);
end
