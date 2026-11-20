function res = simuler_mas_rk4(repere, p_mas, h, tf, Cr_func)
% SIMULER_MAS_RK4 Simule le moteur asynchrone par la methode de Runge-Kutta d'ordre 4
% dans l'un des trois referentiels de Park :
%   'stator'       : axes lies au stator        (wa = 0,  theta_s = 0)
%   'rotor'        : axes lies au rotor         (wa = w,  theta_s = theta)
%   'synchronisme' : axes lies au synchronisme  (wa = ws, theta_s = ws*t)
%
% Entrees :
%   repere  : 'stator', 'rotor' ou 'synchronisme'
%   p_mas   : structure des parametres retournee par parametres_mas()
%   h       : pas de calcul (s)
%   tf      : temps final de simulation (s)
%   Cr_func : fonction handle @(t) retournant le couple de charge Cr(t) (N.m)
%
% Sortie :
%   res     : structure contenant t, isd, isq, ird, irq, isa, isb, isc,
%             ira, irb, irc, Is_eff, Ir_eff, Ce, Omega, N_trmin, theta, Cr, vsd, vsq, vsa

    if nargin < 5 || isempty(Cr_func)
        Cr_func = @(t) 0;
    end

    Nsteps = floor(tf / h) + 1;

    % Preallocation des vecteurs de stockage
    t      = zeros(1, Nsteps);
    isd    = zeros(1, Nsteps);
    isq    = zeros(1, Nsteps);
    ird    = zeros(1, Nsteps);
    irq    = zeros(1, Nsteps);
    Ce     = zeros(1, Nsteps);
    Omega  = zeros(1, Nsteps);
    theta  = zeros(1, Nsteps);
    Cr_vec = zeros(1, Nsteps);
    vsd    = zeros(1, Nsteps);
    vsq    = zeros(1, Nsteps);

    % Conditions initiales : X = [isd; isq; ird; irq; w; theta]
    X = zeros(6, 1);
    ti = 0;

    for n = 1:Nsteps
        t(n)     = ti;
        isd(n)   = X(1);
        isq(n)   = X(2);
        ird(n)   = X(3);
        irq(n)   = X(4);
        w_cur    = X(5);
        th_cur   = X(6);
        Omega(n) = w_cur / p_mas.p;
        theta(n) = th_cur;
        Ce(n)    = p_mas.p * p_mas.M * (X(3)*X(2) - X(4)*X(1));
        Cr_vec(n)= Cr_func(ti);
        U_cur    = get_U(repere, p_mas, ti, th_cur);
        vsd(n)   = U_cur(1);
        vsq(n)   = U_cur(2);

        if n < Nsteps
            % Integration RK4 couplee (courants, vitesse electrique, angle)
            K1 = h * deriv_mas(X,          ti,         repere, p_mas, Cr_func);
            K2 = h * deriv_mas(X + K1/2,   ti + h/2,   repere, p_mas, Cr_func);
            K3 = h * deriv_mas(X + K2/2,   ti + h/2,   repere, p_mas, Cr_func);
            K4 = h * deriv_mas(X + K3,     ti + h,     repere, p_mas, Cr_func);
            X  = X + (1/6) * (K1 + 2*K2 + 2*K3 + K4);
            ti = ti + h;
        end
    end

    % Reconstruction des grandeurs de phase (transformation de Park inverse)
    switch lower(repere)
        case 'stator'
            th_s = zeros(size(t));
            th_r = -theta;
        case 'rotor'
            th_s = theta;
            th_r = zeros(size(t));
        case 'synchronisme'
            th_s = p_mas.ws * t;
            th_r = p_mas.ws * t - theta;
        otherwise
            error('Referentiel inconnu : %s', repere);
    end

    c23 = sqrt(2/3);
    isa = c23 * (isd .* cos(th_s)          - isq .* sin(th_s));
    isb = c23 * (isd .* cos(th_s - 2*pi/3) - isq .* sin(th_s - 2*pi/3));
    isc = c23 * (isd .* cos(th_s - 4*pi/3) - isq .* sin(th_s - 4*pi/3));

    ira = c23 * (ird .* cos(th_r)          - irq .* sin(th_r));
    irb = c23 * (ird .* cos(th_r - 2*pi/3) - irq .* sin(th_r - 2*pi/3));
    irc = c23 * (ird .* cos(th_r - 4*pi/3) - irq .* sin(th_r - 4*pi/3));

    vsa = c23 * (vsd .* cos(th_s)          - vsq .* sin(th_s));

    % Modules efficaces instantanes (conservation de puissance : ||i_dq|| = sqrt(3)*I_eff)
    Is_eff = sqrt((isd.^2 + isq.^2) / 3);
    Ir_eff = sqrt((ird.^2 + irq.^2) / 3);

    % Remplissage de la structure de sortie
    res.repere  = repere;
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
    res.Is_eff  = Is_eff;
    res.Ir_eff  = Ir_eff;
    res.Ce      = Ce;
    res.Omega   = Omega;
    res.N_trmin = Omega * 60 / (2*pi);
    res.theta   = theta;
    res.Cr      = Cr_vec;
    res.vsd     = vsd;
    res.vsq     = vsq;
    res.vsa     = vsa;
end

function dX = deriv_mas(X, t, repere, p_mas, Cr_func)
    I  = X(1:4);
    w  = X(5);
    th = X(6);

    switch lower(repere)
        case 'stator'
            wa = 0;
        case 'rotor'
            wa = w;
        case 'synchronisme'
            wa = p_mas.ws;
    end
    wr = wa - w;

    % Matrice [R] generale du modele de Park
    R = [p_mas.Rs,       -wa*p_mas.Ls,   0,            -wa*p_mas.M;
         wa*p_mas.Ls,     p_mas.Rs,      wa*p_mas.M,    0;
         0,              -wr*p_mas.M,    p_mas.Rr,     -wr*p_mas.Lr;
         wr*p_mas.M,      0,             wr*p_mas.Lr,   p_mas.Rr];

    U = get_U(repere, p_mas, t, th);

    % Systeme d'etat electrique : dI/dt = A*I + B
    A  = -p_mas.Lu * R;
    B  =  p_mas.Lu * U;
    dI = A * I + B;

    % Equation mecanique avec frottement sec physique Cp
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

    dth = w;
    dX = [dI; dw; dth];
end

function U = get_U(repere, p_mas, t, th)
    Vamp = sqrt(3) * p_mas.V;
    switch lower(repere)
        case 'stator'
            U = [Vamp * cos(p_mas.ws * t);
                 Vamp * sin(p_mas.ws * t);
                 0; 0];
        case 'rotor'
            U = [Vamp * cos(p_mas.ws * t - th);
                 Vamp * sin(p_mas.ws * t - th);
                 0; 0];
        case 'synchronisme'
            U = [Vamp;
                 0;
                 0; 0];
    end
end
