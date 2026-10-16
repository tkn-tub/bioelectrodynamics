% run_D1_IFFL  Drive the IFFL differentiator with m_in and plot m_d1.
clc; clear; close all;

dt = 0.5;                                  % s  (keep <= ~1 s for the r_I step)
T  = 120*60;  t = (0:dt:T).';  N = numel(t);

% --- synthetic input: rest, a downward change at 5-15 min, then hold ---
m_in_rest = 0.00125;                       % mM
m_in = m_in_rest*ones(N,1);
ramp = min(max((t/60-5)/10,0),1);          % 0 before 5 min, 1 after 15 min
m_in = m_in_rest*(1 - 0.30*ramp);

% --- to use your real data instead, comment the block above and do: ---
% load Dataset_cell2cell.mat; Dataset = rmmissing(Dataset);
% m_in = <your m_in column at the analysed distance>;  dt = dt0;  t = (0:numel(m_in)-1).'*dt;

[m_d1, aux] = D1_IFFL(m_in, dt);

fprintf('baseline m_d1 = %.3g mM,  peak = %.3g mM (%.0fx baseline),  r_I enzyme = %.3g uM\n', ...
        aux.baseline, max(m_d1), max(m_d1)/aux.baseline, aux.E_rI_uM);

figure;
yyaxis left;  plot(t/60, m_d1, 'LineWidth',2);  ylabel('m_{d1} [mM]');
yyaxis right; plot(t/60, -gradient(m_in,t), ':','LineWidth',1.5); ylabel('-dm_{in}/dt [mM/s]');
xlabel('Time [min]'); grid on; xlim([0 min(120,t(end)/60)]);
legend('m_{d1} (IFFL output)','input falling slope','Location','best');
title('Block 1 differentiator — IFFL primitive');
