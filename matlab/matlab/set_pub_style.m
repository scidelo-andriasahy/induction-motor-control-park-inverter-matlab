function C = set_pub_style()
% SET_PUB_STYLE Configure les reglages graphiques MATLAB pour des figures
% de qualite publication (police, tailles, couleurs harmonieuses) et retourne
% une palette de couleurs adaptee aux rapports scientifiques.

    set(0, 'DefaultFigureColor', 'w');
    set(0, 'DefaultAxesFontName', 'Helvetica');
    set(0, 'DefaultAxesFontSize', 10);
    set(0, 'DefaultAxesLineWidth', 0.9);
    set(0, 'DefaultAxesBox', 'on');
    set(0, 'DefaultAxesXGrid', 'on');
    set(0, 'DefaultAxesYGrid', 'on');
    set(0, 'DefaultAxesGridAlpha', 0.25);
    set(0, 'DefaultLineLineWidth', 1.35);
    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');

    % Palette ingenieur (lisible en couleur et en niveaux de gris)
    C.navy    = [0.043, 0.235, 0.365];  % #0B3C5D Bleu marine profond
    C.crimson = [0.753, 0.224, 0.169];  % #C0392B Rouge brique
    C.emerald = [0.118, 0.518, 0.286];  % #1E8449 Vert emeraude
    C.amber   = [0.839, 0.537, 0.063];  % #D68910 Ambre / Orange
    C.purple  = [0.490, 0.235, 0.596];  % #7D3C98 Violet
    C.teal    = [0.086, 0.627, 0.522];  % #16A085 Turquoise
    C.slate   = [0.365, 0.427, 0.494];  % #5D6D7E Gris ardoise
    C.black   = [0.120, 0.120, 0.120];  % Anthracite
end
