clc
clear
close all

%% Figures
h1 = figure(1);
ha1 = axes;

figure(2);
t = tiledlayout(2,2); 

% Top-left
ha2 = axes(t);
ha2.Layout.Tile = 1;

% Top-right
ha4 = axes(t);
ha4.Layout.Tile = 2;

% Bottom-left
ha6 = axes(t);
ha6.Layout.Tile = 3;

% Bottom-right
ha8 = axes(t);
ha8.Layout.Tile = 4;


figure(3);
t = tiledlayout(2,2);   % 2x2 grid

% Top-left
ha3 = polaraxes(t);
ha3.Layout.Tile = 1;

% Top-right
ha5 = polaraxes(t);
ha5.Layout.Tile = 2;

% Bottom-left
ha7 = polaraxes(t);
ha7.Layout.Tile = 3;

% Bottom-right
ha9 = polaraxes(t);
ha9.Layout.Tile = 4;

%% Fixed variables
dt = 0.001; % sampling time interval
fs = 1/dt; % sampling freq
t = 0:dt:1; % time axis
ls = length(t); % length of signal
freq_components = 0:1:25;
n_fcomp = length(freq_components);
% phi_components = zeros(size(freq_components)); 
phi_components = pi/4*linspace(0, 1, size(freq_components, 2)); % linear
 % phase components
 % phi_components = 2*pi*rand(size(freq_components)); % Random phase components


%% Construct real signal

% phi_components = 2*pi*rand(size(freq_components)); % Random phase components

% signal = @(f, phi) exp(1i*2*pi*f.'*t + phi.');
signal = @(f, phi) cos(2*pi*f.'*t + phi.');

signal_components = signal(freq_components, phi_components);
real_signal = real(sum(signal_components, 1)).';

pd_signal = real_signal - min(real_signal); % positive valued signal (DC offset)

plot(ha1, t, real(signal_components), 'LineWidth', 0.5, 'LineStyle','--');
hold(ha1, "on")
plot(ha1, t, real_signal, 'LineWidth', 3);
plot(ha1, t, pd_signal, 'LineWidth', 3, 'LineStyle', ':', 'DisplayName', ...
    'PD signal');
xlabel(ha1, "Time");
hold(ha1, "off");
ylim(ha1, [-10 30])


%% Hilbert transform

% First, look at FFT of the real signal

fft_phasor_analysis(pd_signal,ls,fs,ha2,ha3, n_fcomp);
title(ha2, 'FFT of real Signal');
title(ha3, 'Polar representation of FFT of real Signal');


% Second, Hilbert transform the PD signal and look at FFT

analytic_signal = hilbert(pd_signal);

fft_phasor_analysis(analytic_signal,ls,fs,ha4,ha5, n_fcomp);
title(ha4, 'FFT of Analytic Signal');
title(ha5, 'Polar representation of FFT of Analytic Signal');


% Understand with distribitive property FFT(x+y) = FFT(x) + FFT(y): Here,
% the FFT(x) represents the FFT of the real signal and FFT(y) represents
% the FFT of the Hilbert transform of the real signal at quadrature.

% https://en.wikipedia.org/wiki/Hilbert_transform#Relationship_with_the_Fourier_transform
% The Hilbert transform of u can be thought of as the convolution of u(t) 
% with the function h(t) = ⁠1/(πt), known as the Cauchy kernel.
% H(u)(t) has the effect of shifting the phase of the negative frequency 
% components of u(t) by +90° (π⁄2 radians) and the phase of the positive 
% frequency components by −90°, and i·H(u)(t) has the effect of restoring 
% the positive frequency components while shifting the negative frequency 
% ones an additional +90°, resulting in their negation (i.e., a 
% multiplication by −1).

analytic_signal_real_comp = real(hilbert(pd_signal));

fft_phasor_analysis(analytic_signal_real_comp,ls,fs,ha6,ha7, n_fcomp);
title(ha6, 'FFT of real part of Analytic Signal');
title(ha7, ["Polar representation of FFT" "of real part of Analytic Signal"]);


analytic_signal_1i_imag_comp = 1i*imag(hilbert(pd_signal));

fft_phasor_analysis(analytic_signal_1i_imag_comp,ls,fs,ha8,ha9, n_fcomp);
title(ha8, 'FFT of imaginary part of Analytic Signal');
title(ha9, ["Polar representation of FFT" "of imaginary part of Analytic Signal"]);

% Note here how the negative frequency components are out of phase and thus
% cancel each oter out, while the positive frequency components have the
% entire signal energy.

%% Helper functions
function fft_phasor_analysis(signal,ls,fs,h_axes_cart,h_axes_polar, n_fcomp)
fft_pd_signal = fftshift(fft(signal));
fa = (-(ls-1)/2:(ls-1)/2)*fs/ls;
yyaxis(h_axes_cart, "left")
stem(h_axes_cart, fa, abs(fft_pd_signal), 'LineWidth', 2, 'Marker', '+');
ylim(h_axes_cart, [0, 10e5]);
hold(h_axes_cart, "on");
yyaxis(h_axes_cart, "right")
stem(h_axes_cart, fa, rad2deg(angle(fft_pd_signal)), 'LineWidth', 2, ...
    'Marker', 'o');
hold(h_axes_cart, "off")
xlim(h_axes_cart, [-n_fcomp+0.5, n_fcomp-0.5]);
ylim(h_axes_cart, [-200, 200]);

magXneg   = abs(fft_pd_signal(1:(ls-1)/2));         % magnitude
phaseXneg = angle(fft_pd_signal(1:(ls-1)/2));       % phase in radians
magXpos   = abs(fft_pd_signal((ls+1)/2:end));         % magnitude
phaseXpos = angle(fft_pd_signal((ls+1)/2:end));       % phase in radians

polarplot(h_axes_polar, phaseXpos, magXpos, '.');   % polar representation
hold(h_axes_polar, "on")
polarplot(h_axes_polar, phaseXneg, magXneg, 'o');   % polar representation
hold(h_axes_polar, "off")
rlim(h_axes_polar, [0 2*max(abs(fft_pd_signal(fa~=0)))])

end
