function res = simuler_mas_gradateur(mode_cmd, params_cmd, p_mas, h, tf, Cr_func)
% SIMULER_MAS_GRADATEUR Simule le moteur asynchrone alimente par un gradateur
% triphase a thyristors tete-beche (montage etoile a neutre isole) dans le
% referentiel de Park lie au stator (wa = 0), resolu par RK4.
%
% Modes de commande (mode_cmd) :
%   'angle_phase' : angle d'amorcage constant alpha_deg (params_cmd.alpha_deg)
%   'rampe_alpha' : demarrage progressif, rampe lineaire d'angle d'amorcage
%                   de alpha_init_deg vers alpha_fin_deg sur t_rampe (s)
%   'train_ondes' : commande par train d'ondes pleines
%                   params_cmd.n_on (nb periodes ON), params_cmd.N_tot (nb periodes totales)
%
% Modele physique exact de commutation naturelle des thyristors :
%   - Mode 3 (3 phases conductrices) : vs = vs_reseau
%   - Mode 2 (phase k bloquee, i_sk=0) : tension composee imposee sur les 2 autres phases,
%             composante selon u_k egale a la f.e.m. interne e_s pour maintenir di_sk/dt = 0
%   - Mode 0 (toutes phases bloquees, i_s=0) : stator en circuit ouvert, decroissance libre
%             du flux rotorique et roue libre mecanique.

    if nargin < 6 || isempty(Cr_func)
        Cr_func = @(t) 0;
    end

    Nsteps = floor(tf / h) + 1;

    % Axes unitaires des 3 phases dans le plan (d,q) lie au stator
    theta_ph = [0, 2*pi/3, 4*pi/3];
    u_ax = [cos(theta_ph); sin(theta_ph)];   % 2x3 : vecteurs unitaires de phase
    t_ax = [-sin(theta_ph); cos(theta_ph)];  % 2x3 : vecteurs orthogonaux (tensions composees)

    % Preallocation
    t        = zeros(1, Nsteps);
    isd      = zeros(1, Nsteps);
    isq      = zeros(1, Nsteps);
    ird      = zeros(1, Nsteps);
    irq      = zeros(1, Nsteps);
    vsd      = zeros(1, Nsteps);
    vsq      = zeros(1, Nsteps);
    vsa      = zeros(1, Nsteps);
    ea_net   = zeros(1, Nsteps);
    Ce       = zeros(1, Nsteps);
    Omega    = zeros(1, Nsteps);
    theta    = zeros(1, Nsteps);
    alpha_tr = zeros(1, Nsteps);
    nb_cond  = zeros(1, Nsteps);

    % Etat initial : X = [isd; isq; ird; irq; w; theta]
    X = zeros(6, 1);
    cond = [false, false, false];  % Etats de conduction des 3 phases
    ti = 0;

    c23 = sqrt(2/3);
    c32 = sqrt(3/2);
    Vpeak = p_mas.V * sqrt(2);
    T_res = 1 / p_mas.fs;

    for n = 1:Nsteps
        % Tensions reseau sinusoïdales (passage par zero a t=0 pour phase a)
        e_abc = Vpeak * [sin(p_mas.ws * ti);
                         sin(p_mas.ws * ti - 2*pi/3);
                         sin(p_mas.ws * ti - 4*pi/3)];
        vs_net = p_mas.P0 * e_abc;

        % Determination des impulsions de gachette selon le mode de commande
        [gate_pos, gate_neg, alpha_rad] = get_gates(mode_cmd, params_cmd, ti, p_mas.ws, T_res);

        % Calcul de la f.e.m. statorique equivalente e_s_emf
        e_emf = calc_emf(X, p_mas);

        % --- Transitions d'amorcage (blocage -> conduction) ---
        % Mode 0 -> Mode 2 (amorcage simultane de deux phases)
        if sum(cond) == 0
            pairs = [1 2; 2 3; 3 1];
            best_v = 0;
            best_pair = [];
            for idx = 1:3
                j = pairs(idx, 1);
                m = pairs(idx, 2);
                v_th_jm = c23 * ((u_ax(:,j) - u_ax(:,m))' * (vs_net - e_emf));
                can_pos = (v_th_jm > 1e-3)  && gate_pos(j) && gate_neg(m);
                can_neg = (v_th_jm < -1e-3) && gate_neg(j) && gate_pos(m);
                if (can_pos || can_neg) && abs(v_th_jm) > best_v
                    best_v = abs(v_th_jm);
                    best_pair = [j, m];
                end
            end
            if ~isempty(best_pair)
                cond(best_pair) = true;
            end
        end

        % Mode 2 -> Mode 3 (amorcage de la 3eme phase)
        if sum(cond) == 2
            k_off = find(~cond, 1);
            v_th_k = c32 * (u_ax(:,k_off)' * (vs_net - e_emf));
            if (v_th_k > 1e-3 && gate_pos(k_off)) || (v_th_k < -1e-3 && gate_neg(k_off))
                cond(k_off) = true;
            end
        end

        % Tension statorique effectivement appliquee dans le mode courant
        vs_app = calc_vs_mode(X, vs_net, cond, p_mas, u_ax, t_ax);

        % Stockage
        t(n)        = ti;
        isd(n)      = X(1);
        isq(n)      = X(2);
        ird(n)      = X(3);
        irq(n)      = X(4);
        vsd(n)      = vs_app(1);
        vsq(n)      = vs_app(2);
        vsa(n)      = c23 * vs_app(1);
        ea_net(n)   = e_abc(1);
        Ce(n)       = p_mas.p * p_mas.M * (X(3)*X(2) - X(4)*X(1));
        Omega(n)    = X(5) / p_mas.p;
        theta(n)    = X(6);
        alpha_tr(n) = rad2deg(alpha_rad);
        nb_cond(n)  = sum(cond);

        if n < Nsteps
            % Courants de phase avant le pas RK4
            i_abc_prev = c23 * (u_ax' * X(1:2));

            % Integration RK4 sur le pas h dans la topologie courante
            X_trial = rk4_step_grad(X, ti, h, cond, p_mas, Cr_func, u_ax, t_ax, Vpeak);
            ti_next = ti + h;

            % --- Detection exacte d'extinction naturelle (passage par zero du courant) ---
            i_abc_new = c23 * (u_ax' * X_trial(1:2));
            [gate_pos_n, gate_neg_n, ~] = get_gates(mode_cmd, params_cmd, ti_next, p_mas.ws, T_res);

            if sum(cond) == 3
                % Recherche de la premiere phase traversant zero sans gachette opposee
                gamma_min = 2.0;
                k_ext = 0;
                for k = 1:3
                    if i_abc_prev(k) * i_abc_new(k) <= 0 && abs(i_abc_prev(k)) > 1e-8
                        no_rev = (i_abc_prev(k) > 0 && ~gate_neg_n(k)) || ...
                                 (i_abc_prev(k) < 0 && ~gate_pos_n(k));
                        if no_rev
                            gam = abs(i_abc_prev(k)) / (abs(i_abc_prev(k)) + abs(i_abc_new(k)) + eps);
                            if gam < gamma_min
                                gamma_min = gam;
                                k_ext = k;
                            end
                        end
                    end
                end
                if k_ext > 0
                    % Sous-pas exact jusqu'a l'annulation du courant, puis fin du pas en mode 2
                    h1 = gamma_min * h;
                    h2 = h - h1;
                    X_mid = rk4_step_grad(X, ti, h1, cond, p_mas, Cr_func, u_ax, t_ax, Vpeak);
                    cond(k_ext) = false;
                    % Annulation residuelle a flux rotorique constant
                    d_is = -(u_ax(:,k_ext)' * X_mid(1:2)) * u_ax(:,k_ext);
                    X_mid(1:2) = X_mid(1:2) + d_is;
                    X_mid(3:4) = X_mid(3:4) - (p_mas.M / p_mas.Lr) * d_is;
                    if h2 > 1e-12
                        X = rk4_step_grad(X_mid, ti + h1, h2, cond, p_mas, Cr_func, u_ax, t_ax, Vpeak);
                    else
                        X = X_mid;
                    end
                else
                    X = X_trial;
                end

            elseif sum(cond) == 2
                idx_on = find(cond);
                j = idx_on(1);
                m = idx_on(2);
                if i_abc_prev(j) * i_abc_new(j) <= 0 && abs(i_abc_prev(j)) > 1e-8
                    can_rev = (i_abc_prev(j) > 0 && gate_neg_n(j) && gate_pos_n(m)) || ...
                              (i_abc_prev(j) < 0 && gate_pos_n(j) && gate_neg_n(m));
                    if ~can_rev
                        gam = abs(i_abc_prev(j)) / (abs(i_abc_prev(j)) + abs(i_abc_new(j)) + eps);
                        h1 = gam * h;
                        h2 = h - h1;
                        X_mid = rk4_step_grad(X, ti, h1, cond, p_mas, Cr_func, u_ax, t_ax, Vpeak);
                        cond = [false, false, false];
                        d_is = -X_mid(1:2);
                        X_mid(1:2) = [0; 0];
                        X_mid(3:4) = X_mid(3:4) - (p_mas.M / p_mas.Lr) * d_is;
                        if h2 > 1e-12
                            X = rk4_step_grad(X_mid, ti + h1, h2, cond, p_mas, Cr_func, u_ax, t_ax, Vpeak);
                        else
                            X = X_mid;
                        end
                    else
                        X = X_trial;
                    end
                else
                    X = X_trial;
                    k_off = find(~cond, 1);
                    d_is = -(u_ax(:,k_off)' * X(1:2)) * u_ax(:,k_off);
                    X(1:2) = X(1:2) + d_is;
                    X(3:4) = X(3:4) - (p_mas.M / p_mas.Lr) * d_is;
                end
            else
                X = X_trial;
                X(1:2) = [0; 0];
            end
            ti = ti_next;
        end
    end

    % Reconstruction des courants de phase
    isa = c23 * isd;
    isb = c23 * (-0.5 * isd + (sqrt(3)/2) * isq);
    isc = c23 * (-0.5 * isd - (sqrt(3)/2) * isq);

    th_r = -theta;
    ira = c23 * (ird .* cos(th_r)          - irq .* sin(th_r));
    irb = c23 * (ird .* cos(th_r - 2*pi/3) - irq .* sin(th_r - 2*pi/3));
    irc = c23 * (ird .* cos(th_r - 4*pi/3) - irq .* sin(th_r - 4*pi/3));

    res.mode     = mode_cmd;
    res.t        = t;
    res.isd      = isd;
    res.isq      = isq;
    res.ird      = ird;
    res.irq      = irq;
    res.isa      = isa;
    res.isb      = isb;
    res.isc      = isc;
    res.ira      = ira;
    res.irb      = irb;
    res.irc      = irc;
    res.Is_eff   = sqrt((isd.^2 + isq.^2) / 3);
    res.Ir_eff   = sqrt((ird.^2 + irq.^2) / 3);
    res.vsd      = vsd;
    res.vsq      = vsq;
    res.vsa      = vsa;
    res.ea_net   = ea_net;
    res.Ce       = Ce;
    res.Omega    = Omega;
    res.N_trmin  = Omega * 60 / (2*pi);
    res.theta    = theta;
    res.alpha_tr = alpha_tr;
    res.nb_cond  = nb_cond;
end

function [gate_pos, gate_neg, alpha_rad] = get_gates(mode_cmd, params_cmd, t, ws, T_res)
    switch lower(mode_cmd)
        case 'angle_phase'
            alpha_rad = deg2rad(params_cmd.alpha_deg);
            [gate_pos, gate_neg] = gates_from_alpha(t, ws, alpha_rad);

        case 'rampe_alpha'
            if t <= params_cmd.t_rampe
                alpha_deg = params_cmd.alpha_init_deg + ...
                    (params_cmd.alpha_fin_deg - params_cmd.alpha_init_deg) * (t / params_cmd.t_rampe);
            else
                alpha_deg = params_cmd.alpha_fin_deg;
            end
            alpha_rad = deg2rad(alpha_deg);
            [gate_pos, gate_neg] = gates_from_alpha(t, ws, alpha_rad);

        case 'train_ondes'
            T_c  = params_cmd.N_tot * T_res;
            T_on = params_cmd.n_on  * T_res;
            t_mod = mod(t, T_c);
            if t_mod < T_on
                gate_pos = [true, true, true];
                gate_neg = [true, true, true];
                alpha_rad = 0;
            else
                gate_pos = [false, false, false];
                gate_neg = [false, false, false];
                alpha_rad = pi;
            end

        otherwise
            error('Mode de commande gradateur inconnu : %s', mode_cmd);
    end
end

function [gate_pos, gate_neg] = gates_from_alpha(t, ws, alpha_rad)
    % Commande classique par impulsions larges (largeur 90 deg) pour gradateur
    % triphase etoile a neutre isole : assure le recouvrement a alpha et alpha+60 deg
    % sans jamais superposer T_k+ et T_k-.
    pulse_width = pi / 2; % 90 deg
    if alpha_rad <= 1e-4
        gate_pos = [true, true, true];
        gate_neg = [true, true, true];
        return;
    end
    gate_pos = false(1, 3);
    gate_neg = false(1, 3);
    phi_k = [0, 2*pi/3, 4*pi/3];
    for k = 1:3
        psi = mod(ws * t - phi_k(k), 2*pi);
        d_pos = mod(psi - alpha_rad, 2*pi);
        d_neg = mod(psi - (pi + alpha_rad), 2*pi);
        gate_pos(k) = (d_pos < pulse_width);
        gate_neg(k) = (d_neg < pulse_width);
    end
end

function e_emf = calc_emf(X, p_mas)
    I = X(1:4);
    w = X(5);
    rr = [p_mas.Rr * I(3) + w * (p_mas.M * I(2) + p_mas.Lr * I(4));
          p_mas.Rr * I(4) - w * (p_mas.M * I(1) + p_mas.Lr * I(3))];
    e_emf = p_mas.Rs * I(1:2) - (p_mas.M / p_mas.Lr) * rr;
end

function vs = calc_vs_mode(X, vs_net, cond, p_mas, u_ax, t_ax)
    nc = sum(cond);
    if nc == 3
        vs = vs_net;
    elseif nc == 2
        k_off = find(~cond, 1);
        e_emf = calc_emf(X, p_mas);
        u = u_ax(:, k_off);
        tv = t_ax(:, k_off);
        vs = (tv' * vs_net) * tv + (u' * e_emf) * u;
    else
        vs = calc_emf(X, p_mas);
    end
end

function X_next = rk4_step_grad(X, ti, dt, cond, p_mas, Cr_func, u_ax, t_ax, Vpeak)
    K1 = dt * deriv_grad(X,        ti,          cond, p_mas, Cr_func, u_ax, t_ax, Vpeak);
    K2 = dt * deriv_grad(X + K1/2, ti + dt/2,   cond, p_mas, Cr_func, u_ax, t_ax, Vpeak);
    K3 = dt * deriv_grad(X + K2/2, ti + dt/2,   cond, p_mas, Cr_func, u_ax, t_ax, Vpeak);
    K4 = dt * deriv_grad(X + K3,   ti + dt,     cond, p_mas, Cr_func, u_ax, t_ax, Vpeak);
    X_next = X + (1/6) * (K1 + 2*K2 + 2*K3 + K4);
end

function dX = deriv_grad(X, t, cond, p_mas, Cr_func, u_ax, t_ax, Vpeak)
    I = X(1:4);
    w = X(5);

    e_abc = Vpeak * [sin(p_mas.ws * t);
                     sin(p_mas.ws * t - 2*pi/3);
                     sin(p_mas.ws * t - 4*pi/3)];
    vs_net = p_mas.P0 * e_abc;
    vs_app = calc_vs_mode(X, vs_net, cond, p_mas, u_ax, t_ax);

    R = [p_mas.Rs,       0,            0,          0;
         0,              p_mas.Rs,     0,          0;
         0,              w*p_mas.M,    p_mas.Rr,   w*p_mas.Lr;
        -w*p_mas.M,      0,           -w*p_mas.Lr, p_mas.Rr];

    U = [vs_app(1); vs_app(2); 0; 0];
    dI = p_mas.Lu * (U - R * I);

    Ce_val = p_mas.p * p_mas.M * (I(3)*I(2) - I(4)*I(1));
    Cr_val = Cr_func(t);
    Omega_val = w / p_mas.p;

    if abs(Omega_val) < 1e-4 && abs(Ce_val - Cr_val) <= p_mas.Cp
        dw = 0;
    else
        if abs(Omega_val) >= 1e-4
            Cp_eff = p_mas.Cp * sign(Omega_val);
        else
            Cp_eff = p_mas.Cp * sign(Ce_val - Cr_val);
        end
        dw = (p_mas.p / p_mas.J) * (Ce_val - p_mas.f * Omega_val - Cp_eff - Cr_val);
    end

    dX = [dI; dw; w];
end
