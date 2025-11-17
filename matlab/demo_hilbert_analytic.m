clc
clear
close all

%% Figures
h1 = figure(1); % real signal
ha1 = axes;

h2 = figure(2); % FFTs
t = tiledlayout(2,2); 

% Top-left
ha2 = axes(t);
ha2.Layout.Tile = 1;

% Top-right
ha103 = axes(t);
ha103.Layout.Tile = 2;

% Bottom-left
ha104 = axes(t);
ha104.Layout.Tile = 3;

% Bottom-right
ha105 = axes(t);
ha105.Layout.Tile = 4;


h3 = figure(3); % FFT Polar plots
t = tiledlayout(2,2);   % 2x2 grid

% Top-left
ha106 = polaraxes(t);
ha106.Layout.Tile = 1;

% Top-right
ha107 = polaraxes(t);
ha107.Layout.Tile = 2;

% Bottom-left
ha108 = polaraxes(t);
ha108.Layout.Tile = 3;

% Bottom-right
ha109 = polaraxes(t);
ha109.Layout.Tile = 4;


h4 = figure(4); % FFT quiver plot
t = tiledlayout(2,2);   % 2x2 grid

% Top Left
ha110 = axes(t);
ha110.Layout.Tile = 1;

% Top Right
ha111 = axes(t);
ha111.Layout.Tile = 2;

% Bottom Left
ha112 = axes(t);
ha112.Layout.Tile = 3;

% Bottom Right
ha113 = axes(t);
ha113.Layout.Tile = 4;

h5 = figure(5); % envelope, phase
t = tiledlayout(3,1);   % 3x1 grid

% Top 
ha114 = axes(t);
ha114.Layout.Tile = 1;

% Middle
ha115 = axes(t);
ha115.Layout.Tile = 2;

% Bottom
ha116 = axes(t);
ha116.Layout.Tile = 3;



%% Fixed variables
dt = 0.001; % sampling time interval
fs = 1/dt; % sampling freq
t = 0:dt:1; % time axis
ls = length(t); % length of signal
freq_components = 0:1:2;
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

fft_phasor_analysis(pd_signal,ls,fs,n_fcomp,ha2,ha106,ha110);
title(ha2, 'FFT of real Signal');
title(ha106, 'Polar representation of FFT of Real Signal');
title(ha110, ["Cartesian representation of "
 "FFT of real part of Real Signal"]);


% Second, Hilbert transform the PD signal and look at FFT

analytic_signal = hilbert(pd_signal);

fft_phasor_analysis(analytic_signal,ls,fs,n_fcomp,ha103,ha107,ha111);
title(ha103, 'FFT of Analytic Signal');
title(ha107, 'Polar representation of FFT of Analytic Signal');
title(ha111, ["Cartesian representation of" "FFT of Analytic Signal"]);


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

fft_phasor_analysis(analytic_signal_real_comp,ls,fs,n_fcomp,ha104,ha108,ha112);
title(ha104, 'FFT of real part of Analytic Signal');
title(ha108, ["Polar representation of FFT" "of real part of Analytic Signal"]);
title(ha112, ["Cartesian representation FFT of" "real part of Analytic Signal"]);

analytic_signal_1i_imag_comp = 1i*imag(hilbert(pd_signal));

fft_phasor_analysis(analytic_signal_1i_imag_comp,ls,fs,n_fcomp,ha105,ha109,ha113);
title(ha105, 'FFT of imaginary part of Analytic Signal');
title(ha109, ["Polar representation of FFT" "of imaginary part of Analytic Signal"]);
title(ha113, ["Cartesian representation FFT of" "imaginary part of Analytic Signal"]);

% Note here how the negative frequency components are out of phase and thus
% cancel each oter out, while the positive frequency components have the
% entire signal energy.


%% Show envelope and phase

envelope = abs(analytic_signal);
phase = unwrap(angle(analytic_signal));

plot(ha114, t, pd_signal); title(ha114, 'Real Signal');
xlabel(ha114, 'Time (s)'); ylabel(ha114, 'x(t)');

plot(ha115, t, envelope); title(ha115, 'Envelope = |analytic|');
xlabel(ha115, 'Time (s)'); ylabel(ha115, '|x(t)|');

plot(ha116, t, phase); title(ha116, 'Instantaneous Phase');
xlabel(ha116, 'Time (s)'); ylabel(ha116, 'angle(x(t))');




%% Helper functions
function fft_phasor_analysis(signal,ls,fs, n_fcomp, h_axes_fft,h_axes_polar, h_axes_cart)
fft_pd_signal = fftshift(fft(signal));
fa = (-(ls-1)/2:(ls-1)/2)*fs/ls;
yyaxis(h_axes_fft, "left")
stem(h_axes_fft, fa, abs(fft_pd_signal), 'LineWidth', 2, ...
    'Marker', '^');
ylabel(h_axes_fft, "Magnitude");
ylim(h_axes_fft, [0, 5*max(abs(fft_pd_signal(fa~=0)))]);
hold(h_axes_fft, "on");
yyaxis(h_axes_fft, "right")
stem(h_axes_fft, fa, rad2deg(angle(fft_pd_signal)), 'LineWidth', 2, ...
    'Marker', 'o');
ylabel(h_axes_fft, "Phase");
hold(h_axes_fft, "off")
xlim(h_axes_fft, [-n_fcomp+0.5, n_fcomp-0.5]);
ylim(h_axes_fft, [-200, 200]);
xlabel(h_axes_fft, "Frequency (Hz)");

magXneg   = abs(fft_pd_signal(1:(ls-1)/2));         % magnitude
phaseXneg = angle(fft_pd_signal(1:(ls-1)/2));       % phase in radians
magXpos   = abs(fft_pd_signal((ls+1)/2:end));         % magnitude
phaseXpos = angle(fft_pd_signal((ls+1)/2:end));       % phase in radians

polarplot(h_axes_polar, phaseXpos, magXpos, '.', 'MarkerSize', 12);   % polar representation
hold(h_axes_polar, "on")
polarplot(h_axes_polar, phaseXneg, magXneg, 'o');   % polar representation
hold(h_axes_polar, "off")
rlim(h_axes_polar, [0 2*max(abs(fft_pd_signal(fa~=0)))])

realXneg = (real(fft_pd_signal(1:(ls-1)/2)));     % real
imagXneg = (imag(fft_pd_signal(1:(ls-1)/2)));     % imag
realXpos = (real(fft_pd_signal((ls+1)/2:end)));     % real
imagXpos = (imag(fft_pd_signal((ls+1)/2:end)));     % imag

quiver(h_axes_cart, zeros((ls+1)/2,1), zeros((ls+1)/2,1), ...
    realXpos, imagXpos, "LineWidth", 1, "AutoScale", "on");
   % cartesian representation
hold(h_axes_cart, "on")
quiver(h_axes_cart, zeros((ls-1)/2,1), zeros((ls-1)/2,1), ...
    realXneg, imagXneg, "LineWidth", 1, "AutoScale", "on");
   % cartesian representation
hold(h_axes_cart, "off")
xlabel(h_axes_cart, "Real axis");
ylabel(h_axes_cart, "Imaginary axis");

% a = max(max(max(realXpos, imagXpos)), max(max(realXneg, imagXneg)));
a = 2*max(abs(fft_pd_signal(fa~=0)));
xlim(h_axes_cart, [-a, a]);
ylim(h_axes_cart, [-a, a]);
axis(h_axes_cart, "square");
grid(h_axes_cart, "on");

end