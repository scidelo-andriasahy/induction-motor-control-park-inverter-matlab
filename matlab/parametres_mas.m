function p_mas = parametres_mas()
% PARAMETRES_MAS Retourne les parametres electriques et mecaniques du moteur asynchrone
% TP de Commande des Machines Electriques - ENSA Tanger (3eme annee GSEA)
% Caracteristiques nominales : 50 Hz ; 220/380 V ; 3 kVA ; 1500 tr/min ; p = 2

    % Alimentation et plaque signaletique
    p_mas.fs = 50;                     % Frequence du reseau (Hz)
    p_mas.ws = 2 * pi * p_mas.fs;      % Pulsation statorique (rad/s) = 100*pi
    p_mas.V  = 220;                    % Tension simple efficace (V)
    p_mas.U  = 380;                    % Tension composee efficace (V)
    p_mas.Sn = 3000;                   % Puissance apparente nominale (VA)
    p_mas.p  = 2;                      % Nombre de paires de poles
    p_mas.Ns = 60 * p_mas.fs / p_mas.p;% Vitesse de synchronisme (tr/min) = 1500
    p_mas.Ws_mec = p_mas.ws / p_mas.p; % Vitesse angulaire synchrone (rad/s) = 50*pi

    % Parametres electriques (modele de Park a conservation de puissance)
    p_mas.Rs = 1.090;                  % Resistance statorique (Ohm)
    p_mas.Ls = 0.2165;                 % Inductance cyclique statorique (H)
    p_mas.Rr = 0.543;                  % Resistance rotorique (Ohm)
    p_mas.Lr = 0.0246;                 % Inductance cyclique rotorique (H)
    p_mas.M  = 0.0615;                 % Inductance mutuelle cyclique (H)

    % Parametres mecaniques
    p_mas.J  = 0.08;                   % Moment d'inertie (kg.m^2)
    p_mas.f  = 0.0064;                 % Coefficient de frottement visqueux (N.m.s/rad)
    p_mas.Cp = 0.65;                   % Couple de pertes mecaniques / frottement sec (N.m)

    % Matrices d'inductance constantes
    p_mas.L  = [p_mas.Ls, 0,        p_mas.M,  0;
                0,        p_mas.Ls, 0,        p_mas.M;
                p_mas.M,  0,        p_mas.Lr, 0;
                0,        p_mas.M,  0,        p_mas.Lr];
    p_mas.Lu = inv(p_mas.L);

    % Parametres derives utiles
    p_mas.Delta = p_mas.Ls * p_mas.Lr - p_mas.M^2;
    p_mas.sigma = 1 - (p_mas.M^2) / (p_mas.Ls * p_mas.Lr);
    p_mas.Ts = p_mas.Ls / p_mas.Rs;
    p_mas.Tr = p_mas.Lr / p_mas.Rr;

    % Matrice de transformation directe Concordia / Park (theta = 0, puissance conservee)
    p_mas.P0 = sqrt(2/3) * [1, -1/2,       -1/2;
                            0,  sqrt(3)/2, -sqrt(3)/2];
end
