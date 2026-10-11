%[text] # Simulation parameters
% distance between cells
dist=(2:0.01:10)*1e-6;
dd0=dist(2)-dist(1);
% Tumor cell radious
a_T=650e-8;
m_T=200;%total of released molecules
% Immune radious
a_I=4e-6;

%Diffusion Coefficient of Tumor extracellular vesicle (EV)
D_T=66e-13;
%Diffusion Coefficient of Immune EV
D_I=33e-13;
% Degradation rate of tumor EV
kd_T=.2;
% Degradation rate of immune EV
kd_I=0.1;

% Reaction constants on the tumor cell
kb_T=110e-6;
% kb_T=2.8e-5;
kr_T=0.2;
ki_T=0.02;

% Release coefficient of tumor EV
kre_T=-0.5;
% Release coefficient of immune EV
kre_I=kre_T/2;
% Reaction constants on tumor cells
kb_I=60e-6;
% kb_I=2.8e-5;

% Recycling rate of tumor derived extracellular vesicles (TDEV) on immune cells
kr_I=0.1;
% Recycling rate of inmmune derived extracellular vesicles (IDEV) on tumor cells
ki_I=0.005;

% Frequency range  and step
domega=1e-4;
omega0=1e-10:domega:10;
% domega=omega0(2)-omega0(1);
% Time range and step
t0=linspace(0,1/domega,length(omega0));
dt0=t0(2)-t0(1);

% Natural release rate of tumor in time and frequency
a1=-6642;
b1=21.92;
d1=8.771;
a2=6667;
b2=21.94;
d2=8.788;

%Natural release mechanism of TDEV in time and frequency
R_T=m_T*exp(-3e-5*t0);
R_Tf=(fft(R_T,2*length(t0)+1));

%Natural release mechanism of IDEV in time and frequency
R_I=50*exp(-1.5e-5*t0);
R_If=(fft(R_I,2*length(t0)+1));

%[text] #### Plotting parameters
fontsize=15;

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
