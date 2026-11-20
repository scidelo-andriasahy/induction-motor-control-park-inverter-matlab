function stats = tp1_en_charge()
% TP1_EN_CHARGE Manipulation I - Travail facultatif :
% Simulation en charge avec application d'un echelon de couple resistant Cr a t = 2 s,
% pour differentes valeurs (zone stable + etude du decrochage).

    script_dir = fileparts(mfilename('fullpath'));
    fig_dir = fullfile(script_dir, '..', 'figures');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

    C = set_pub_style();
    p = parametres_mas();
    h = 5e-4;
    tf = 4.5;
    t_ch = 2.0;

    Cr_vals = [0, 3, 6, 9, 12];
    cols = {C.slate, C.navy, C.emerald, C.amber, C.crimson};
    res_list = cell(1, length(Cr_vals));

    stats.Cr_vals = Cr_vals;
    for k = 1:length(Cr_vals)
        Cr_k = Cr_vals(k);
        r = simuler_mas_rk4('stator', p, h, tf, @(t) Cr_k * (t >= t_ch));
        res_list{k} = r;
        stats.N_fin(k)      = r.N_trmin(end);
        stats.s_fin_pct(k)  = 100 * (p.Ns - r.N_trmin(end)) / p.Ns;
        stats.Ce_fin(k)     = r.Ce(end);
        stats.Is_eff_fin(k) = r.Is_eff(end);
        stats.Ir_eff_fin(k) = r.Ir_eff(end);
        stats.Pmec_fin(k)   = r.Ce(end) * r.Omega(end);
    end

    % --- Figure 1 : Reponse a l'echelon de charge Cr = 0, 3, 6, 9, 12 N.m ---
    fig1 = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 23 15.5]);
    tl1 = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl1, '\textbf{Manipulation I (Facultatif) -- Comportement en charge (\''echelon de $C_r$ \`a $t = 2\,\mathrm{s}$)}', ...
          'FontSize', 12, 'Interpreter', 'latex');

    lbls = arrayfun(@(x) sprintf('$C_r = %d\\,\\mathrm{N\\cdot m}$', x), Cr_vals, 'UniformOutput', false);

    % (a) Couple electromagnetique
    nexttile(tl1); hold on;
    for k = 1:length(Cr_vals)
        plot(res_list{k}.t, res_list{k}.Ce, 'Color', cols{k}, 'LineWidth', 1.3);
    end
    xline(t_ch, '--', 'Color', C.black, 'LineWidth', 1.0);
    text(t_ch-0.78, 28, '$t = 2\,\mathrm{s}$', 'FontSize', 8.5);
    xlabel('Temps $t$ (s)'); ylabel('Couple $C_e$ ($\mathrm{N\cdot m}$)');
    title('(a) Couple \''electromagn\''etique $C_e(t)$', 'FontSize', 10);
    legend(lbls, 'Location', 'northeast', 'FontSize', 7.8);
    xlim([0 tf]); ylim([-3 33]);

    % (b) Vitesse de rotation
    nexttile(tl1); hold on;
    for k = 1:length(Cr_vals)
        plot(res_list{k}.t, res_list{k}.N_trmin, 'Color', cols{k}, 'LineWidth', 1.4);
    end
    xline(t_ch, '--', 'Color', C.black, 'LineWidth', 1.0);
    xlabel('Temps $t$ (s)'); ylabel('Vitesse $N$ (tr/min)');
    title('(b) Vitesse $N(t)$ et augmentation du glissement', 'FontSize', 10);
    legend(lbls, 'Location', 'southeast', 'FontSize', 7.8);
    xlim([0 tf]); ylim([0 1600]);

    % (c) Courant efficace statorique
    nexttile(tl1); hold on;
    for k = 1:length(Cr_vals)
        plot(res_list{k}.t, res_list{k}.Is_eff, 'Color', cols{k}, 'LineWidth', 1.4);
    end
    xline(t_ch, '--', 'Color', C.black, 'LineWidth', 1.0);
    xlabel('Temps $t$ (s)'); ylabel('$I_{s,\mathrm{eff}}$ (A)');
    title('(c) Courant efficace statorique $I_{s,\mathrm{eff}}(t)$', 'FontSize', 10);
    legend(lbls, 'Location', 'northeast', 'FontSize', 7.8);
    xlim([0 tf]); ylim([0 15]);

    % (d) Courant efficace rotorique
    nexttile(tl1); hold on;
    for k = 1:length(Cr_vals)
        plot(res_list{k}.t, res_list{k}.Ir_eff, 'Color', cols{k}, 'LineWidth', 1.4);
    end
    xline(t_ch, '--', 'Color', C.black, 'LineWidth', 1.0);
    xlabel('Temps $t$ (s)'); ylabel('$I_{r,\mathrm{eff}}$ (A)');
    title('(d) Courant efficace rotorique $I_{r,\mathrm{eff}}(t)$', 'FontSize', 10);
    legend(lbls, 'Location', 'northeast', 'FontSize', 7.8);
    xlim([0 tf]); ylim([0 44]);

    % Zoom vitesse apres t=2s place APRES tous les nexttile
    ax_v = axes(fig1, 'Position', [0.635 0.61 0.17 0.13]);
    box on; hold(ax_v, 'on');
    idx_v = res_list{1}.t >= 1.9 & res_list{1}.t <= 4.5;
    for k = 1:length(Cr_vals)
        plot(ax_v, res_list{k}.t(idx_v), res_list{k}.N_trmin(idx_v), 'Color', cols{k}, 'LineWidth', 1.2);
    end
    xlim(ax_v, [1.9 4.5]); ylim(ax_v, [1270 1510]);
    title(ax_v, 'Zoom en charge', 'FontSize', 7.5);
    set(ax_v, 'FontSize', 6.8);

    exportgraphics(fig1, fullfile(fig_dir, 'fig_tp1_en_charge.png'), 'Resolution', 300);
    close(fig1);

    % --- Figure 2 : Etude de la limite de stabilite et du decrochage ---
    r_12 = simuler_mas_rk4('stator', p, h, 5.5, @(t) 12 * (t >= t_ch));
    r_15 = simuler_mas_rk4('stator', p, h, 5.5, @(t) 15 * (t >= t_ch));
    r_16 = simuler_mas_rk4('stator', p, h, 5.5, @(t) 16 * (t >= t_ch));
    stats.N_fin_15 = r_15.N_trmin(end);
    stats.N_fin_16 = r_16.N_trmin(end);

    % Caracteristique statique theorique
    N_grid = linspace(0, p.Ns - 0.1, 500);
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

    fig2 = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 22 9.0]);
    tl2 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl2, '\textbf{Limite de stabilit\''e en charge et ph\''enom\`ene de d\''ecrochage ($C_r > C_{\max}$)}', ...
          'FontSize', 11.5, 'Interpreter', 'latex');

    nexttile; hold on;
    plot(r_12.t, r_12.N_trmin, 'Color', C.emerald, 'LineWidth', 1.5);
    plot(r_15.t, r_15.N_trmin, 'Color', C.amber,   'LineWidth', 1.5);
    plot(r_16.t, r_16.N_trmin, 'Color', C.crimson, 'LineWidth', 1.6);
    xline(t_ch, '--', 'Color', C.black, 'LineWidth', 1.0);
    xlabel('Temps $t$ (s)'); ylabel('Vitesse $N$ (tr/min)');
    title('(a) Effondrement de la vitesse lors d''une surcharge');
    legend({'$C_r = 12\,\mathrm{N\cdot m}$ (stable)', ...
            '$C_r = 15\,\mathrm{N\cdot m}$ (proche $C_{\max}$)', ...
            '$C_r = 16\,\mathrm{N\cdot m}$ (d\''ecrochage !)'}, ...
           'Location', 'southwest', 'FontSize', 8.5);
    xlim([0 5.5]); ylim([0 1600]);

    nexttile; hold on;
    plot(N_grid, Ce_stat, '--', 'Color', C.black, 'LineWidth', 1.4);
    idx_12 = r_12.t >= t_ch;
    idx_16 = r_16.t >= t_ch;
    plot(r_12.N_trmin(idx_12), r_12.Ce(idx_12), 'Color', C.emerald, 'LineWidth', 1.6);
    plot(r_16.N_trmin(idx_16), r_16.Ce(idx_16), 'Color', C.crimson, 'LineWidth', 1.5);
    xlabel('Vitesse $N$ (tr/min)'); ylabel('Couple $C_e$ ($\mathrm{N\cdot m}$)');
    title('(b) Trajectoires dans le plan $(N, C_e)$ apr\`es $t = 2\,\mathrm{s}$');
    legend({'Caract\''eristique statique', 'Trajectoire $C_r = 12\,\mathrm{N\cdot m}$', ...
            'Trajectoire $C_r = 16\,\mathrm{N\cdot m}$ (d\''ecrochage)'}, ...
           'Location', 'northwest', 'FontSize', 8.5);
    xlim([0 1550]); ylim([0 20]);

    exportgraphics(fig2, fullfile(fig_dir, 'fig_tp1_decrochage.png'), 'Resolution', 300);
    close(fig2);
end
