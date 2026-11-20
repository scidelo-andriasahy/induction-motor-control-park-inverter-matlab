function stats = tp2_gradateur_train_ondes()
% TP2_GRADATEUR_TRAIN_ONDES Manipulation II - Partie II :
% Moteur asynchrone alimente par gradateur triphase a thyristors
% en commande par train d'ondes pleines (periodes entieres n_on / N_tot).

    script_dir = fileparts(mfilename('fullpath'));
    fig_dir = fullfile(script_dir, '..', 'figures');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

    C = set_pub_style();
    p = parametres_mas();
    h = 5e-5;
    tf = 5.0;

    % Configurations [n_on, N_tot] (T_c = N_tot * 20 ms)
    configs = [8, 8;
               7, 8;
               6, 8;
               5, 8;
               4, 8];
    cols = {C.black, C.navy, C.emerald, C.amber, C.crimson};
    res_to = cell(1, size(configs, 1));

    for k = 1:size(configs, 1)
        n_on  = configs(k, 1);
        N_tot = configs(k, 2);
        r = simuler_mas_gradateur('train_ondes', struct('n_on', n_on, 'N_tot', N_tot), p, h, tf);
        res_to{k} = r;
        tail = r.t >= 3.5;
        stats.n_on(k)       = n_on;
        stats.N_tot(k)      = N_tot;
        stats.beta(k)       = n_on / N_tot;
        stats.N_moy(k)      = mean(r.N_trmin(tail));
        stats.dN_pp(k)      = max(r.N_trmin(tail)) - min(r.N_trmin(tail));
        stats.isa_max(k)    = max(abs(r.isa));
        stats.Is_eff_max(k) = max(r.Is_eff);
        stats.Is_eff_moy(k) = mean(r.Is_eff(tail));
        stats.Ce_max(k)     = max(r.Ce);
        stats.Ce_min_ss(k)  = min(r.Ce(tail));
        stats.Ce_max_ss(k)  = max(r.Ce(tail));
        % FFT sur un multiple exact de T_c (ex. 8*T_c = 1.28 s = 64 periodes reseau)
        [~, ~, thd_i] = calc_fft_thd(r.t, r.isa, 3.68, p.fs, 64);
        stats.thd_i(k) = thd_i;
    end

    stats.r_68 = res_to{3}; % 6/8 (beta = 0.75) pour comparaison globale
    stats.r_48 = res_to{5}; % 4/8 (beta = 0.50)

    % --- Figure 1 : Grandeurs globales pour differents rapports cycliques beta ---
    fig1 = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 23.5 15.5]);
    tl1 = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl1, '\textbf{Manipulation II -- Gradateur en commande par train d''ondes pleines ($T_c = 8\,T = 160\,\mathrm{ms}$)}', ...
          'FontSize', 11.5, 'Interpreter', 'latex');

    lbls = {'$8/8$ ($\beta = 1{,}00$, direct)', '$7/8$ ($\beta = 0{,}875$)', ...
            '$6/8$ ($\beta = 0{,}75$)', '$5/8$ ($\beta = 0{,}625$)', '$4/8$ ($\beta = 0{,}50$)'};
    step_p = 10;

    % (a) Couple electromagnetique
    nexttile(tl1); hold on;
    for k = [1, 3, 5]
        plot(res_to{k}.t(1:step_p:end), res_to{k}.Ce(1:step_p:end), 'Color', cols{k}, 'LineWidth', 1.0);
    end
    xlabel('Temps $t$ (s)'); ylabel('Couple $C_e$ ($\mathrm{N\cdot m}$)');
    title('(a) Couple $C_e(t)$ (chocs \`a chaque salve)', 'FontSize', 9.8);
    legend(lbls([1, 3, 5]), 'Location', 'northeast', 'FontSize', 8);
    xlim([0 tf]); ylim([-22 34]);

    % (b) Vitesse mecanique
    nexttile(tl1); hold on;
    for k = 1:size(configs, 1)
        plot(res_to{k}.t(1:step_p:end), res_to{k}.N_trmin(1:step_p:end), 'Color', cols{k}, 'LineWidth', 1.3);
    end
    xlabel('Temps $t$ (s)'); ylabel('Vitesse $N$ (tr/min)');
    title('(b) Vitesse $N(t)$ (ondulations basse fr\''equence)', 'FontSize', 9.8);
    legend(lbls, 'Location', 'southeast', 'FontSize', 7.5);
    xlim([0 tf]); ylim([0 1620]);

    % (c) Courant efficace statorique
    nexttile(tl1); hold on;
    for k = [1, 3, 5]
        plot(res_to{k}.t(1:step_p:end), res_to{k}.Is_eff(1:step_p:end), 'Color', cols{k}, 'LineWidth', 1.0);
    end
    xlabel('Temps $t$ (s)'); ylabel('$I_{s,\mathrm{eff}}$ (A)');
    title('(c) Courant statorique $I_{s,\mathrm{eff}}(t)$ (r\''eappels)', 'FontSize', 9.8);
    legend(lbls([1, 3, 5]), 'Location', 'northeast', 'FontSize', 8);
    xlim([0 tf]); ylim([0 19]);

    % (d) Courant efficace rotorique
    nexttile(tl1); hold on;
    for k = [1, 3, 5]
        plot(res_to{k}.t(1:step_p:end), res_to{k}.Ir_eff(1:step_p:end), 'Color', cols{k}, 'LineWidth', 1.0);
    end
    xlabel('Temps $t$ (s)'); ylabel('$I_{r,\mathrm{eff}}$ (A)');
    title('(d) Courant rotorique $I_{r,\mathrm{eff}}(t)$', 'FontSize', 9.8);
    legend(lbls([1, 3, 5]), 'Location', 'northeast', 'FontSize', 8);
    xlim([0 tf]); ylim([0 46]);

    exportgraphics(fig1, fullfile(fig_dir, 'fig_tp2_train_ondes_global.png'), 'Resolution', 300);
    close(fig1);

    % --- Figure 2 : Zoom sur 2 cycles de modulation Tc = 320 ms pour beta = 6/8 ---
    r68 = res_to{3};
    iz = r68.t >= 3.84 & r68.t <= 4.16; % Exactement 2 cycles Tc = 2 * 0.16 s

    fig2 = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 23.5 15.0]);
    tl2 = tiledlayout(fig2, 2, 2, 'TileSpacing', 'normal', 'Padding', 'compact');
    title(tl2, '\textbf{Zoom sur deux cycles de modulation $T_c = 160\,\mathrm{ms}$ (Train d''ondes $6/8$, $\beta = 0{,}75$)}', ...
          'FontSize', 11.5, 'Interpreter', 'latex');

    % (a) Tension statorique vsa(t) et reseau ea(t)
    nexttile(tl2);
    plot(r68.t(iz), r68.ea_net(iz), ':', 'Color', C.slate, 'LineWidth', 0.9); hold on;
    plot(r68.t(iz), r68.vsa(iz),    '-', 'Color', C.navy,  'LineWidth', 1.25);
    xlabel('Temps $t$ (s)'); ylabel('Tension (V)');
    title('(a) Tension $v_{sa}(t)$ ($T_{\mathrm{on}}=120\,\mathrm{ms}$, $T_{\mathrm{off}}=40\,\mathrm{ms}$)', 'FontSize', 9.5);
    legend({'R\''eseau $e_a(t)$', 'Stator $v_{sa}(t)$ (f.\''e.m. en $T_{\mathrm{off}}$)'}, ...
           'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 8);
    xlim([3.84 4.16]); ylim([-350 350]);

    % (b) Courant statorique isa(t)
    nexttile(tl2);
    plot(r68.t(iz), r68.isa(iz), 'Color', C.crimson, 'LineWidth', 1.2); hold on;
    yline(0, ':', 'Color', C.black, 'LineWidth', 0.8);
    xlabel('Temps $t$ (s)'); ylabel('Courant $i_{sa}$ (A)');
    title('(b) Courant $i_{sa}(t)$ (extinction au z\''ero puis pic)', 'FontSize', 9.5);
    xlim([3.84 4.16]); ylim([-30 30]);

    % (c) Couple electromagnetique Ce(t)
    nexttile(tl2);
    plot(r68.t(iz), r68.Ce(iz), 'Color', C.emerald, 'LineWidth', 1.25); hold on;
    yline(0, ':', 'Color', C.black, 'LineWidth', 0.8);
    xlabel('Temps $t$ (s)'); ylabel('Couple $C_e$ ($\mathrm{N\cdot m}$)');
    title('(c) Pulsations de couple $C_e(t)$ (roue libre puis choc)', 'FontSize', 9.5);
    xlim([3.84 4.16]); ylim([-18 28]);

    % (d) Ondulation de vitesse N(t) et courant rotorique ira(t)
    nexttile(tl2);
    yyaxis left;
    plot(r68.t(iz), r68.N_trmin(iz), 'Color', C.navy, 'LineWidth', 1.5);
    ylabel('Vitesse $N$ (tr/min)', 'Color', C.navy);
    ax = gca; ax.YAxis(1).Color = C.navy;

    yyaxis right;
    plot(r68.t(iz), r68.ira(iz), 'Color', C.amber, 'LineWidth', 1.1);
    ylabel('Courant rotorique $i_{ra}$ (A)', 'Color', C.amber);
    ax.YAxis(2).Color = C.amber;

    xlabel('Temps $t$ (s)');
    title('(d) Vitesse $N(t)$ et courant rotorique $i_{ra}(t)$', 'FontSize', 9.5);
    xlim([3.84 4.16]);

    exportgraphics(fig2, fullfile(fig_dir, 'fig_tp2_train_ondes_zoom.png'), 'Resolution', 300);
    close(fig2);
end
