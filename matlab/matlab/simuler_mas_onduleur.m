function res = simuler_mas_onduleur(mode_ond, params_ond, p_mas, h, tf, Cr_func)
% SIMULER_MAS_ONDULEUR Simule le moteur asynchrone alimente par un onduleur de tension
% triphase a deux niveaux dans le referentiel de Park lie au stator (wa = 0), par RK4.
%
% Modes d'onduleur (mode_ond) :
%   'pleine_onde' : commande pleine onde 180 deg (6 pas), frequence f_ref
%                   E_dc ajuste pour meme fondamental (ou specifie dans params_ond.Edc)
%   'mli'         : commande MLI sinus-triangle (SPWM) avec loi scalaire V/f = cte
%                   params_ond.f_func : @(t) frequence de reference f_ref(t) (Hz)
%                   params_ond.fp     : frequence porteuse triangulaire (Hz), ex. 1050 Hz
%                   params_ond.Edc    : tension bus continu (defaut 2*sqrt(2)*220 = 622.25 V)

    if nargin < 6 || isempty(Cr_func)
        Cr_func = @(t) 0;
    end

    Nsteps = floor(tf / h) + 1;

    t      = zeros(1, Nsteps);
    isd    = zeros(1, Nsteps);
    isq    = zeros(1, Nsteps);
    ird    = zeros(1, Nsteps);
    irq    = zeros(1, Nsteps);
    vsd    = zeros(1, Nsteps);
    vsq    = zeros(1, Nsteps);
    vsa    = zeros(1, Nsteps);
    uab    = zeros(1, Nsteps);
    Ce     = zeros(1, Nsteps);
    Omega  = zeros(1, Nsteps);
    theta  = zeros(1, Nsteps);
    f_tr   = zeros(1, Nsteps);

    X = zeros(6, 1);  % [isd; isq; ird; irq; w; theta]
    ti = 0;
    phi_ref = 0;      % Phase integree de la reference sinusoidale (pour f(t) variable)

    % Matrice de conversion etats d'interrupteurs (Sa,Sb,Sc) -> tensions simples (va,vb,vc)
    M_ond = (1/3) * [ 2, -1, -1;
                     -1,  2, -1;
                     -1, -1,  2];

    for n = 1:Nsteps
        f_cur = params_ond.f_func(ti);
        [vs_abc, vs_dq] = calc_onduleur(mode_ond, params_ond, ti, phi_ref, f_cur, p_mas, M_ond);

        t(n)     = ti;
        isd(n)   = X(1);
        isq(n)   = X(2);
        ird(n)   = X(3);
        irq(n)   = X(4);
        vsd(n)   = vs_dq(1);
        vsq(n)   = vs_dq(2);
        vsa(n)   = vs_abc(1);
        uab(n)   = vs_abc(1) - vs_abc(2);
        Ce(n)    = p_mas.p * p_mas.M * (X(3)*X(2) - X(4)*X(1));
        Omega(n) = X(5) / p_mas.p;
        theta(n) = X(6);
        f_tr(n)  = f_cur;

        if n < Nsteps
            % Pour un pas h tres petit (ex. 25 us), la tension d'onduleur est constante
            % sur le pas d'echantillonnage (ou evaluee aux points RK4).
            f_mid  = params_ond.f_func(ti + h/2);
            f_end  = params_ond.f_func(ti + h);
            phi_mid = phi_ref + 2*pi * (0.5*(f_cur + f_mid)) * (h/2);
            phi_end = phi_ref + 2*pi * (0.5*(f_cur + f_end)) * h;

            [~, dq1] = calc_onduleur(mode_ond, params_ond, ti,       phi_ref, f_cur, p_mas, M_ond);
            [~, dq2] = calc_onduleur(mode_ond, params_ond, ti + h/2, phi_mid, f_mid, p_mas, M_ond);
            [~, dq4] = calc_onduleur(mode_ond, params_ond, ti + h,   phi_end, f_end, p_mas, M_ond);

            K1 = h * deriv_ond(X,        ti,       dq1, p_mas, Cr_func);
            K2 = h * deriv_ond(X + K1/2, ti + h/2, dq2, p_mas, Cr_func);
            K3 = h * deriv_ond(X + K2/2, ti + h/2, dq2, p_mas, Cr_func);
            K4 = h * deriv_ond(X + K3,   ti + h,   dq4, p_mas, Cr_func);

            X = X + (1/6) * (K1 + 2*K2 + 2*K3 + K4);
            phi_ref = phi_end;
            ti = ti + h;
        end
    end

    c23 = sqrt(2/3);
    isa = c23 * isd;
    isb = c23 * (-0.5 * isd + (sqrt(3)/2) * isq);
    isc = c23 * (-0.5 * isd - (sqrt(3)/2) * isq);

    th_r = -theta;
    ira = c23 * (ird .* cos(th_r)          - irq .* sin(th_r));
    irb = c23 * (ird .* cos(th_r - 2*pi/3) - irq .* sin(th_r - 2*pi/3));
    irc = c23 * (ird .* cos(th_r - 4*pi/3) - irq .* sin(th_r - 4*pi/3));

    res.mode    = mode_ond;
    res.t       = t;
    res.isd     = isd;
    res.isq     = isq;
    res.ird     = ird;
    res.irq     = irq;
    res.isa     = isa;
    res.isb     = isb;
    res.isc     = isc;
    res.ira     = ira;
    res.irb     = irb;
    res.irc     = irc;
    res.Is_eff  = sqrt((isd.^2 + isq.^2) / 3);
    res.Ir_eff  = sqrt((ird.^2 + irq.^2) / 3);
    res.vsd     = vsd;
    res.vsq     = vsq;
    res.vsa     = vsa;
    res.uab     = uab;
    res.Ce      = Ce;
    res.Omega   = Omega;
    res.N_trmin = Omega * 60 / (2*pi);
    res.theta   = theta;
    res.f_ref   = f_tr;
end

function [vs_abc, vs_dq] = calc_onduleur(mode_ond, params_ond, t, phi_ref, f_cur, p_mas, M_ond)
    switch lower(mode_ond)
        case 'pleine_onde'
            % Commande 180 deg : Sa = 1 si sin(phi) >= 0, etc.
            ref = [sin(phi_ref);
                   sin(phi_ref - 2*pi/3);
                   sin(phi_ref - 4*pi/3)];
            S = double(ref >= 0);
            % Tension bus continu adaptee a la loi V/f si souhaite, sinon Edc fixe
            if isfield(params_ond, 'use_vf') && params_ond.use_vf
                Edc_eff = params_ond.Edc * (f_cur / p_mas.fs);
            else
                Edc_eff = params_ond.Edc;
            end
            vs_abc = Edc_eff * (M_ond * S);

        case 'mli'
            % Porteuse triangulaire symetrique dans [-1, +1] de frequence fp
            fp = params_ond.fp;
            tau = mod(t * fp, 1);
            if tau < 0.25
                tri = 4 * tau;
            elseif tau < 0.75
                tri = 2 - 4 * tau;
            else
                tri = 4 * tau - 4;
            end

            % Indice de modulation r suivant la loi V/f = cte (avec boost Rs faible si besoin)
            r = min(1.0, max(0.05, f_cur / p_mas.fs));
            ref = r * [sin(phi_ref);
                       sin(phi_ref - 2*pi/3);
                       sin(phi_ref - 4*pi/3)];
            S = double(ref >= tri);
            vs_abc = params_ond.Edc * (M_ond * S);

        otherwise
            error('Mode onduleur inconnu : %s', mode_ond);
    end
    vs_dq = p_mas.P0 * vs_abc;
end

function dX = deriv_ond(X, t, vs_dq, p_mas, Cr_func)
    I = X(1:4);
    w = X(5);

    R = [p_mas.Rs,       0,            0,          0;
         0,              p_mas.Rs,     0,          0;
         0,              w*p_mas.M,    p_mas.Rr,   w*p_mas.Lr;
        -w*p_mas.M,      0,           -w*p_mas.Lr, p_mas.Rr];

    U = [vs_dq(1); vs_dq(2); 0; 0];
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
