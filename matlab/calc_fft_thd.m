function [f_axis, amp, thd_pct, f1_amp] = calc_fft_thd(t, y, t_start, f0, n_cycles)
% CALC_FFT_THD Calcule le spectre harmonique (FFT unilatere) et le THD (%)
% sur un nombre entier exact de periodes du fondamental f0 (aucune fuite spectrale).

    dt = t(2) - t(1);
    T_win = n_cycles / f0;
    N_pts = round(T_win / dt);
    idx0 = find(t >= t_start, 1);
    if idx0 + N_pts - 1 > length(y)
        idx0 = length(y) - N_pts + 1;
    end
    seg = y(idx0 : idx0 + N_pts - 1);
    seg = seg - mean(seg); % Suppression composante continue eventuelle

    Y = fft(seg);
    P2 = abs(Y / N_pts);
    P1 = P2(1 : floor(N_pts/2) + 1);
    P1(2:end-1) = 2 * P1(2:end-1);

    fs_ech = 1 / dt;
    f_axis = fs_ech * (0 : floor(N_pts/2)) / N_pts;

    % Fondamental a f0
    [~, idx_f1] = min(abs(f_axis - f0));
    f1_amp = P1(idx_f1);

    % Energie harmonique jusqu'a 40*f0 (ou Nyquist)
    f_max = min(50 * f0, f_axis(end));
    mask_h = (f_axis > 1.5 * f0) & (f_axis <= f_max);
    % Calcul du THD a partir de la valeur efficace totale hors fondamental (y compris interharmoniques)
    harm_rms_sq = sum(P1(mask_h).^2);
    % Si la fenetre contient des sous-harmoniques (ex. train d'ondes), on inclut 0 < f < f0
    mask_sub = (f_axis > 0) & (f_axis < 0.8 * f0);
    harm_rms_sq = harm_rms_sq + sum(P1(mask_sub).^2);

    if f1_amp > 1e-9
        thd_pct = 100 * sqrt(harm_rms_sq) / f1_amp;
    else
        thd_pct = NaN;
    end
    amp = P1;
end
