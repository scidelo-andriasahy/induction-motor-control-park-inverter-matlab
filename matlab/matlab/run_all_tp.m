% RUN_ALL_TP Script maitre executant l'ensemble des simulations MATLAB des
% Manipulations I et II (y compris tous les travaux facultatifs), generant
% les 13 figures haute resolution (300 DPI) dans ../figures/ et sauvegardant
% le resume numerique exact dans results_summary.txt.

clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);

tic;
fprintf('=================================================================\n');
fprintf('  TP DE COMMANDE DES MACHINES ELECTRIQUES (ENSA Tanger - GSEA3)\n');
fprintf('  Execution complete des simulations MATLAB (Manip. I & II)\n');
fprintf('=================================================================\n\n');

p = parametres_mas();

fprintf('[1/8] Manipulation I-I   : Demarrage direct (repere stator)...\n');
st1 = tp1_repere_stator();

fprintf('[2/8] Manipulation I-II  : Demarrage direct (repere rotor)...\n');
st2 = tp1_repere_rotor();

fprintf('[3/8] Manipulation I-III : Demarrage direct (repere synchronisme)...\n');
st3 = tp1_repere_synchronisme();

fprintf('[4/8] Manipulation I-IV  : Comparaison des 3 referentiels...\n');
st_comp = tp1_comparaison_reperes(st1.res, st2.res, st3.res);

fprintf('[5/8] Manipulation I-3   : Fonctionnement en charge & decrochage...\n');
st_ch = tp1_en_charge();

fprintf('[6/8] Manipulation II-I  : Gradateur - Commande par angle de phase...\n');
st_ap = tp2_gradateur_angle_phase();

fprintf('[7/8] Manipulation II-II : Gradateur - Commande par train d ondes...\n');
st_to = tp2_gradateur_train_ondes();

fprintf('[8/8] Manipulation II-III & II-3 : Comparaison & Onduleur (PO/MLI/Vf)...\n');
tp2_comparaison_gradateurs(st_ap, st_to);
st_ond = tp2_onduleur();

elapsed = toc;

% Ecriture du resume numerique complet dans results_summary.txt
fid = fopen(fullfile(script_dir, 'results_summary.txt'), 'w');

fprintf(fid, '=== PARAMETRES DERIVES DU MOTEUR ===\n');
fprintf(fid, 'sigma = %.6f | Delta = %.6e H^2 | Ts = %.4f s | Tr = %.4f s\n', p.sigma, p.Delta, p.Ts, p.Tr);
fprintf(fid, 'Ns = %.1f tr/min | Ws_mec = %.4f rad/s | ws = %.4f rad/s\n', p.Ns, p.Ws_mec, p.ws);
fprintf(fid, 'Ce_dem_stat = %.2f Nm | Ce_decrochage_stat = %.2f Nm (a N = %.1f tr/min)\n\n', ...
    st1.Ce_dem_stat, st1.Ce_decrochage_stat, st1.N_decrochage_stat);

fprintf(fid, '=== MANIPULATION I : DEMARRAGE DIRECT A VIDE (3 REPERES) ===\n');
fprintf(fid, 't_dem_95 = %.3f s | Ce_max = %.2f Nm | Ce_min = %.2f Nm | Ce_fin = %.3f Nm\n', ...
    st1.t_dem_95, st1.Ce_max, st1.Ce_min, st1.Ce_fin);
fprintf(fid, 'N_fin = %.2f tr/min | Omega_fin = %.2f rad/s | s_fin = %.2f %% | f_gliss = %.3f Hz\n', ...
    st1.N_fin, st1.Omega_fin, st1.s_fin_pct, st2.f_gliss);
fprintf(fid, 'Is_max_pic_3ph = %.2f A | isa_max (cos) = %.2f A | Is_eff_max = %.2f A | Is_eff_fin = %.2f A (rapport %.2f)\n', ...
    st1.Is_max_pic, st1.isa_max, st1.Is_eff_max, st1.Is_eff_fin, st1.Is_eff_max/st1.Is_eff_fin);
fprintf(fid, 'Ir_eff_max = %.2f A | Ir_eff_fin = %.2f A | ira_max = %.2f A\n', ...
    st1.Ir_eff_max, st1.Ir_eff_fin, st1.ira_max);
fprintf(fid, 'Synchronisme DC final: isd=%.3f A, isq=%.3f A, ird=%.3f A, irq=%.3f A\n', ...
    st3.isd_fin, st3.isq_fin, st3.ird_fin, st3.irq_fin);
fprintf(fid, 'Ecart max Stator-Rotor : dCe=%.2e Nm, dN=%.2e tr/min, disa=%.2e A\n', ...
    st_comp.err_Ce_st_ro, st_comp.err_N_st_ro, st_comp.err_isa_st_ro);
fprintf(fid, 'Ecart max Stator-Sync  : dCe=%.2e Nm, dN=%.2e tr/min, disa=%.2e A\n\n', ...
    st_comp.err_Ce_st_sy, st_comp.err_N_st_sy, st_comp.err_isa_st_sy);

fprintf(fid, '=== MANIPULATION I (FACULTATIF) : FONCTIONNEMENT EN CHARGE (t=2s) ===\n');
for k = 1:length(st_ch.Cr_vals)
    fprintf(fid, 'Cr=%2d Nm | N_fin=%7.2f tr/min | s=%5.2f %% | Ce_fin=%5.2f Nm | Is_eff=%5.2f A | Ir_eff=%5.2f A | Pmec=%7.1f W\n', ...
        st_ch.Cr_vals(k), st_ch.N_fin(k), st_ch.s_fin_pct(k), st_ch.Ce_fin(k), ...
        st_ch.Is_eff_fin(k), st_ch.Ir_eff_fin(k), st_ch.Pmec_fin(k));
end
fprintf(fid, 'Surcharge: Cr=15 Nm -> N_fin=%.2f tr/min | Cr=16 Nm -> N_fin=%.2f tr/min (decrochage)\n\n', ...
    st_ch.N_fin_15, st_ch.N_fin_16);

fprintf(fid, '=== MANIPULATION II-I : GRADATEUR ANGLE DE PHASE ===\n');
for k = 1:length(st_ap.alphas)
    fprintf(fid, 'alpha=%3d deg | N_fin=%7.2f tr/min | t95=%.3f s | Vsa_rms=%5.1f V | isa_max=%5.2f A | Is_eff_max=%5.2f A | Ce_max=%5.2f Nm | THD_i=%5.1f %%\n', ...
        st_ap.alphas(k), st_ap.N_fin(k), st_ap.t95(k), st_ap.Vsa_rms_fin(k), ...
        st_ap.isa_max(k), st_ap.Is_eff_max(k), st_ap.Ce_max(k), st_ap.thd_i(k));
end
fprintf(fid, 'Rampe 110->70 (2s) | N_fin=%7.2f tr/min | t95=%.3f s | isa_max=%5.2f A | Is_eff_max=%5.2f A | Ce_max=%5.2f Nm | THD_i=%5.1f %%\n\n', ...
    st_ap.rampe.N_fin, st_ap.rampe.t95, st_ap.rampe.isa_max, st_ap.rampe.Is_eff_max, st_ap.rampe.Ce_max, st_ap.rampe.thd_i);

fprintf(fid, '=== MANIPULATION II-II : GRADATEUR TRAIN D ONDES ===\n');
for k = 1:length(st_to.beta)
    fprintf(fid, 'Train %d/%d (beta=%.3f) | N_moy=%7.2f tr/min | dN_pp=%5.1f tr/min | isa_max=%5.2f A | Is_eff_max=%5.2f A | Ce_max=%5.2f Nm | Ce_ss=[%6.2f, %6.2f] | THD_i=%5.1f %%\n', ...
        st_to.n_on(k), st_to.N_tot(k), st_to.beta(k), st_to.N_moy(k), st_to.dN_pp(k), ...
        st_to.isa_max(k), st_to.Is_eff_max(k), st_to.Ce_max(k), st_to.Ce_min_ss(k), st_to.Ce_max_ss(k), st_to.thd_i(k));
end
fprintf(fid, '\n');

fprintf(fid, '=== MANIPULATION II-3 (FACULTATIF) : ONDULEUR ===\n');
fprintf(fid, 'Pleine onde 180 (50Hz) | N_fin=%7.2f tr/min | Ce_max=%5.2f Nm | dCe_pp=%5.3f Nm | std_Ce=%5.3f Nm | isa_max=%5.2f A | THD_i=%5.1f %%\n', ...
    st_ond.po.N_fin, st_ond.po.Ce_max, st_ond.po.dCe_pp, st_ond.po.std_Ce, st_ond.po.isa_max, st_ond.po.thd_i);
fprintf(fid, 'MLI SPWM 2kHz  (50Hz) | N_fin=%7.2f tr/min | Ce_max=%5.2f Nm | dCe_pp=%5.3f Nm | std_Ce=%5.3f Nm | isa_max=%5.2f A | THD_i=%5.1f %%\n', ...
    st_ond.mli.N_fin, st_ond.mli.Ce_max, st_ond.mli.dCe_pp, st_ond.mli.std_Ce, st_ond.mli.isa_max, st_ond.mli.thd_i);
fprintf(fid, 'MLI V/f=cte | N(25Hz)=%6.2f tr/min | N(37.5Hz)=%7.2f tr/min | N(50Hz)=%7.2f tr/min | Is_eff_max=%5.2f A | Ce_max=%5.2f Nm\n', ...
    st_ond.vf.N_25, st_ond.vf.N_375, st_ond.vf.N_50, st_ond.vf.Is_max, st_ond.vf.Ce_max);

fclose(fid);

fprintf('\nTermine avec succes en %.1f s !\n', elapsed);
fprintf('Toutes les figures sont enregistrees dans ../figures/\n');
fprintf('Resume numerique enregistre dans matlab/results_summary.txt\n');
