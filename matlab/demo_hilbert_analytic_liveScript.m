clc
clear
close all

%% Figures
h1 = figure(1);
ha1 = axes; %[output:0a4bb1e8]

h2 = figure(2);
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


h3 = figure(3);
t = tiledlayout(2,2);   % 2x2 grid %[output:0a4bb1e8]

% Top-left
ha3 = polaraxes(t); %[output:0a4bb1e8]
ha3.Layout.Tile = 1; %[output:0a4bb1e8]

% Top-right
ha5 = polaraxes(t); %[output:0a4bb1e8]
ha5.Layout.Tile = 2; %[output:0a4bb1e8]

% Bottom-left
ha7 = polaraxes(t); %[output:0a4bb1e8]
ha7.Layout.Tile = 3; %[output:0a4bb1e8]

% Bottom-right
ha9 = polaraxes(t); %[output:0a4bb1e8]
ha9.Layout.Tile = 4; %[output:0a4bb1e8]


figure(4);
t = tiledlayout(3,1);   % 2x2 grid %[output:0a4bb1e8]

% Top 
ha10 = axes(t); %[output:0a4bb1e8]
ha10.Layout.Tile = 1; %[output:0a4bb1e8]

% Middle
ha11 = axes(t); %[output:0a4bb1e8]
ha11.Layout.Tile = 2; %[output:0a4bb1e8]

% Bottom
ha12 = axes(t); %[output:0a4bb1e8]
ha12.Layout.Tile = 3; %[output:0a4bb1e8]


%% Fixed variables
dt = 0.001; % sampling time interval
fs = 1/dt; % sampling freq
t = 0:dt:1; % time axis
ls = length(t); % length of signal
freq_components = 0:1:6;
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


%% Show envelope and phase

envelope = abs(analytic_signal);
phase = unwrap(angle(analytic_signal));

plot(ha10, t, pd_signal); title('Real Signal'); ylabel('x(t)');

plot(ha11, t, envelope); title('Envelope = |analytic|');

plot(ha12, t, phase); title('Instantaneous Phase');
xlabel('Time (s)');




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

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":23.1}
%---
%[output:0a4bb1e8]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAbYAAAEICAYAAAAzydF1AAAAAXNSR0IArs4c6QAAIABJREFUeF7tnXmsTVf7xx9DUf6gN1xctJW3YkrfmoVrSuMVVaUUiV6VmmK+iOFWQtuUmEpKkSCoPwSpVs3kvdUbU0pMkYiQEhJzK695vIZfnpXfPtnnOuc6995z9jl77c+TNNhnD+v5PE\/P96y113pWqbS0tJeCQQACEIAABCwhUAphsySSuAEBCEAAAoYAwkYiQAACEICAVQQQNqvCiTMQgAAEIICwkQMQgAAEIGAVAYTNqnDiDAQgAAEIIGzkAAQgAAEIWEUAYbMqnDgDAQhAAAIIGzkAAQhAAAJWEUDYrAonzkAAAhCAAMJGDkAAAhCAgFUEEDarwokzEIAABCCAsJEDEIAABCBgFQGEzapw4gwEIAABCCBs5AAEIAABCFhFAGGzKpw4AwEIQAACCBs5AAEIQAACVhFA2KwKJ85AAAIQgADCRg5AAAIQgIBVBBA2q8KJMxCAAAQggLCRAxCAAAQgYBUBhM2qcOIMBCAAAQggbOQABCAAAQhYRQBhsyqcOAMBCEAAAggbOQABCEAAAlYR8L2w1axZU8qWLSuXLl2yKjA4AwEIQAACxSPga2F76623ZOfOnfLy5Utp27Zt8QhwFQQgAAEIWEXA18I2f\/58yczMRNisSkmcgQAEIFAyAr4Ttl69eknTpk2lefPmkpaWZrynx1ayJOBqCEAAAjYR8J2w7d27V8qVKxcWA4TNppTEFwhAAAIlI+A7YRsxYoSkp6cbr+vVqyfvvfcePbaS5QBXQwACELCKgO+EzU1\/8ODBMmzYMITNqpTEGQhAAAIlIxAIYWvWrFkYpePHj5eMGldDAAIQgEDKErBe2IYOHSpDhgwJC8Do0aMFcUvZnKRhEIAABEpEIDDC1qZNmxKB4mIIQAACEPAHAc+FrXTp0vLOO+\/IgwcP5O+\/\/y4RpVjesTk9NoStRKi5GAIQgIBvCHgmbCpos2bNkg4dOkipUqUMoKdPn8rp06dl3Lhx5u+x2OLFi8V5Z6b3ce714sULc\/n27dtl9uzZoVshbLFQ5RwIQAAC9hDwTNhWrlwpjRs3DpFTIVKxU7t48aJkZWWJI06F4d26datUq1Yt6in79u2TnJwchM2eHMUTCEAAAkUi4ImwtW\/fXubNm2cadvLkSRk\/frw8fvxYBgwYIDqRQ00Fy93TiuaFs0B7+vTp0qpVK\/nkk0\/MdP+RI0eaS27cuCHXr19H2IqUBpwMAQhAwB4Cngjb6tWrpWHDhvLw4UPp0qWLPH\/+PETQ+ezOnTvStWvXQsm+8cYboj0yp9II79jsSUQ8gQAEIBAvAp4Im1bg10r8BYcJ1Yn+\/ftLdna28UfrQLp7WwWd1KFMHdLUiSedO3cWR9h0CFOLIUcy3rHFK1W4DwQgAAF\/EPBE2FTQtLelQ4065Oi22rVry8aNG80hHaI8fPhwVHL9+vWTCRMmiPburly5InXq1DF1I\/XfR48elTlz5kh+fn7Y9QibPxKRVkIAAhCIF4GEC1vlypVl9+7dUYXLGV7UE\/Q93G+\/\/RbVt6lTp0qPHj2ifq49uc8\/\/zxsGQELtOOVKtwHAhCAgD8IJFzYWrduLQsXLjQ0VHQuXLjwCpk\/\/\/zTHFu+fLmsWbMmKrlly5bJBx98YD7Pzc2VHTt2mB7aoEGDpEWLFua49uT69OkTuocjbM4kFf2AqiP+SE5aCQEIQKA4BBIubCpEKkhqX375pZw9ezasnRUqVJC8vDxzLNJQpfvk7t27mzVshw4dkt9\/\/z1sobcOQ3bs2NGcrsKmAqfGUGRx0oJrIAABCPiXQMKFTSeN6OQRtUmTJsnBgwfDaGkVkg0bNphjr3vHpudEW+j9119\/hdbJ6c7av\/76K8Lm37yk5RCAAASKTSDhwqYtc4Ya169fLz\/++GNYY3WK\/zfffGOOdevWTW7duhXVmSpVqsgPP\/wgDRo0CJ3jXujtHHT3\/OixFTs3uBACEICALwl4Imzae8rIyDCLpz\/99NMwUE5FEhU0FbZo5u7Z6TkFF3qPGjUqVF6rb9++cvnyZXpsvkxJGg0BCECgZAQ8ETZ9tzZ8+HDTUvcEEV2s\/e233xpB+uWXX2TBggXmHJ0pqTMky5QpIz\/\/\/LMcOHDAHNc\/9Zgu8P7444\/NNH81nfavvUH9rOCaNnpsJUsQrvaGwNKlS82D3JOcvHkyT4FA7AR09G3VqlVmPXEqmyfCpmvNVKCqV69uWNy7d8+U1HJqPmrvSntZjrnfy+n7t0WLFpmP9uzZIxUrVjR\/16LJV69eNaKowubUndSqJL179w4t9EbYUjn9aJtDAGEjF\/xAAGErECUVpHXr1oXEzfn4n3\/+MTUj7969G7rCvfbNLWzOQm\/tsWnvzG0qlDrDUs09CQVh88P\/LrQRYSMH\/EAAYYsSpfT0dPnwww9NT2vXrl1y+\/btmOLpFruvvvrKCFv9+vXNYuwTJ07IpUuXTMkuNfdCb0fYWLsWE2ZOShIBZysm8jRJAeCxMRHQPGUoMiZUsZ1U3IXeKmxNmzaN7SGcBQEIQAAChRJQYUv1H2CevGOLR57Ec6F3PNrDPSAAAQhAIDUJ+EbY4r3QOzXDQasgAAEIQKCkBHwjbOpovBZ6lxQa10MAAhCAQOoS8JWwxWOhd+qGgpZBAAIQgEA8CPhK2Iq60DsegLgHBCAAAQj4i4CvhK2oC739FQpaCwEIQAAC8SDgK2FTh4uy0FuF8N133zUVSu7fvx8PXtwDAkUioBVxtM6pboKray4xCPiZgF\/y2XfC5iRFYQu9a9SoYRZp16tXL5RD+sWiC8KdepR+Ti7anvoEom2vdPr0aRk3bpwpCReL6ZZPWpwgmk2ZMuWVraBiuS\/nQKA4BHr16iWacwU3dC7OvRJ5jW+FLRqUSpUqyZYtW0T\/VNOiyFrlRP9T27Rpk3z\/\/feJZMq9IWCKxDZu3DhEwr290sWLFyUrK8vk5uvMmQkc7bzvvvvO\/GDDIJBoAlqcfvPmzZKWloawJRp2wft\/\/fXX8tFHH5nDP\/30k\/mC0fJb2lNr2bKlOR7LhqZet5vn2UOgffv2ZsRAreD2Sk71\/q1bt5od4wuzunXrmvqqT548Mb28SHb+\/HmG2e1JnZTzRMVs0KBBZg\/MJk2ayJtvvmnaSI\/N41Dt3btX9N3aqVOnZNiwYaGnq7jl5uaawBw5ckSys7M9bhmPCwqB1atXS8OGDeXhw4eiWzNp0W7HnM90yyXdZLcw062Zpk2bJteuXTM7VmAQ8JpArVq1zJZiBQ1h8zAS+uvCKYQcaYhmyZIl0rx5c3n06JEpxIxBIBEE9L2YVsrRXMzJyQl7RP\/+\/UM\/qvR9xfXr16M2YeLEidKnTx85evSojB07NhFN5Z4QKJSAfqdOnTo1dE6bNm2kSpUq9Ni8zJv3339fVqxYYR7ZoUMHyc\/PD3u89uAGDx4sumdb27ZtvWwazwoQAWd7JR1q1CFHt9WuXVs2btxoDr1uSHzx4sXSokUL0XdyajphqmzZsqJbPWmxAt1cF4OAlwScd8f02Dyk7vwajiZczowebVLHjh1jnpnmoQs8yucE3NsrRRIu96iCe3ulSG7rEJAOBUWzc+fOyRdffOFzYjTfTwQQtiREa\/r06dKtWzfzsr1Tp06vtKBz584yY8YMc7xnz56sK0pCjGx\/ZHG3V4rExdkxXkcedBLJ\/v37zUa9OgElIyPDXBLLJBTbmeOfdwQQNu9Yh540adIk+eyzz0xPTHtkBU0\/03PU2rVrF\/ZSPwnN5ZEWEojn9kojRowwQrZ27VrR2Y9u27Ztm1StWjVqrluIFpdSgADCloQgDBgwwPyajTYUOXLkSBk4cKBZP5SZmZmEFvJI2wl4tb3SkCFDRDfRVdORCC1AgEEg0QQQtkQTjnB\/HX501gZpweSzZ8+GnTV37lwzqUS\/BPTLAINAIgjEY3slXbKi5eOePXsWcZ2ae\/Qh0kSpRPjFPSGAsCUhB\/TLIC8vT7SckVZj0Cn\/bvvjjz\/MOrYDBw7I5MmTk9BCHhkEAvHYXskZXdDRBx1WLzjDd9GiRdKqVSuGIoOQUCnkI8KWpGA44PULQVfMO702XYvRo0cP0yodrjx+\/HiSWshjbSdQ1O2VdJmKU0xAf4zdvHlTnKojyurgwYOhd8P6b12DOXPmTFMmLtJaOdv54l\/yCCBsSWKvXwhaSqt8+fKmBbrmR3tyTiFZrSM5Z86cJLWOxwaBQFG3V3LeDSsb\/TF25swZg2nVqlXSqFEj83cdPte1Q5rHOqFETWf\/du\/enZJaQUiqFPERYUtiILSu2fLly42gOaY9uGPHjlHBIYlxCdKji7K9khZEHjNmjMHjFjbNX61xqou0C5qK3KhRo1iyEqSkSgFfHWG7fPmy9O3bNwVaFLkJ1lX3d7upAqcVRrTWnq4JinWrkJSNFg3zHYHCtleK1Rm9h75Pe\/vtt820\/xMnTiBoscLjvEASsFrYAhlRnIYABCAQcAIIW8ATAPchAAEI2EYAYbMtovgDAQhAIOAEELaAJwDuQwACELCNAMJmW0TxBwIQgEDACSBsAU8A3IcABCBgGwGEzbaI4g8EIACBgBNA2AKeALgPAQhAwDYCCJttEcUfCEAAAgEngLAFPAFwHwIQgIBtBBA22yKKPxCAAAQCTgBhC3gC4D4EIAAB2wggbLZFFH8gAAEIBJwAwhbwBMB9CEAAArYRQNhsiyj+QAACEAg4AYQt4AmA+xCAAARsI4Cw2RZR\/IEABCAQcAIIW8ATAPchAAEI2EYAYbMtovgDAQhAIOAEELaAJwDuQwACELCNAMJmW0TxBwIQgEDACSBsAU8A3IcABCBgGwGEzbaI4g8EIACBgBNA2AKeALgPAQhAwDYCCJttEcUfCEAAAgEngLAFPAFwHwIQgIBtBBA22yKKPxCAAAQCTgBhC3gC4D4EIAAB2wggbLZFFH8gAAEIBJwAwhbwBMB9CEAAArYRQNhsiyj+QAACEAg4Ad8LW82aNaVs2bJy6dKlgIcS9yEAAQhAQAn4Wtjeeust2blzp7x8+VLatm1LRCEAAQhAAAL+Frb58+dLZmYmwkYiQwACEIBAiIDvemy9evWSpk2bSvPmzSUtLc04Qo+NjIYABCAAAYeA74Rt7969Uq5cubAIImwkNAQgAAEI+FbYRowYIenp6ab99erVk\/fee48eG\/kMAQhAAAL+HYp0x27w4MEybNgwhI2EhgAEIACB4Ajb0KFDZciQIWEhHz16tBw\/fpw0gAAEIAABCwn47h1bUXtsjrCpmDmGqFmYybgEAQhA4P8JeC5spUuXlnfeeUcePHggf\/\/9d4kCEctQpCNsbdq0KdGzuBgCEIAABPxBwDNhU0GbNWuWdOjQQUqVKmXoPH36VE6fPi3jxo0zf4\/FFi9eLM2aNTOn6n2ce7148cIc2759u8yePTt0K4QtFqqcAwEIQMAeAp4J28qVK6Vx48YhcipEKnZqFy9elKysLHHEqTC8W7dulWrVqkU9Zd++fZKTk4Ow2ZOjeAIBCECgSAQ8Ebb27dvLvHnzTMNOnjwp48ePl8ePH8uAAQPEefelguXuaUXzwlnHNn36dGnVqpV88sknZlbkyJEjzSU3btyQ69evI2xFSgNOhgAEIGAPAU+EbfXq1dKwYUN5+PChdOnSRZ4\/fx4i6Hx2584d6dq1a6Fk33jjDdEembMgm3ds9iQinkAAAhCIFwFPhE0LFWvB4oLDhOpE\/\/79JTs72\/ij5bLcva2CTupQpg5p6sSTzp07C8IWrzTgPhCAAATsIeCJsKmgaW9Lhxp1yNFttWvXlo0bN5pDOkR5+PDhqHT79esnEyZMEO3dXblyRf71r39J+fLlzfkqnnPmzJH8\/Pyw65k8Yk+y4gkEIACBWAgkXNgqV64su3fvjipczvCinqDv4X777beo7Z46dar06NEj6ufak\/v888\/DlhE4wuZeu7Zq1SoWaMeSHZwDAQhAwIcEEi5srVu3loULFxo0KjoXLlx4BdOff\/5pji1fvlzWrFkTFeOyZcvkgw8+MJ\/n5ubKjh07TA9t0KBB0qJFC3Nce3J9+vQJ3QNh82FW0mQIQAACJSCQcGFTIVJBUvvyyy\/l7NmzYc2tUKGC5OXlmWORhirdJ3fv3t2sYTt06JD8\/vvvYQu9dRiyY8eO5nQVNhU4NYYiS5AdXAoBCEDAhwQSLmzOLtfKZtKkSXLw4MEwTFqFZMOGDebY696x6TnRFnr\/9ddfoXVyugHpr7\/+irD5MCFpMgQgAIGSEki4sGkDnaHG9evXy48\/\/hjWZp3i\/80335hj3bp1k1u3bkX1qUqVKvLDDz9IgwYNQue4F3o7B909P3psJU0RrocABCDgLwKeCJv2njIyMszi6U8\/\/TSMkFORRAVNhS2auXt2ek7Bhd6jRo0Kldfq27evXL58mR6bv3KR1kIAAhCICwFPhE3frQ0fPtw02D1BRBdrf\/vtt0aQfvnlF1mwYIE5R2dK6gzJMmXKyM8\/\/ywHDhwwx\/VPPaYLvD\/++GMz7V+tTp06or1B\/Ux7cJmZmSE49NjikifcJMEENE\/V9IceBoFUJbB06VI5ceJEyuepJ8JWrlw5I1DVq1c38bp3754pqeXUfNTelfayHHO\/l9P3b4sWLTIf7dmzRypWrGj+rkWTr169akRRhc2pO6lVSXr37h1a6I2wper\/IrTLTUC\/MNTc2ytBCAKpRkBfK+lyqVT\/AeaJsGlwVJDWrVsXEjcnYP\/884+pGXn37t1QDN1r39zC5iz01h6b9s7cpkKpMyzV3JNQELZU+1+D9kQigLCRF34ggLBFiVJ6erp8+OGHpqe1a9cuuX37dkzxdIvdV199ZYStfv36ZjG2do0vXbpkSnapuRd6R9poNKYHchIEPCTg7PKuv4YxCKQqAf0BRo8tjtEp7kJvR9ji2BRuBQEIQCCwBBC2OIY+ngu949gsbgUBCEAAAilGwLN3bCX1O94LvUvaHq6HAAQgAIHUJOAbYVN88VronZqhoFUQgAAEIBAPAr4Stngs9I4HNO4BAQhAAAKpS8BXwlbUhd6pi52WQQACEIBAogj4StiKutA7UdC4LwQgAAEIpC4BXwmbYizKQu\/UxU7LIAABCEAgUQR8J2wOiFgWemsP79133zWlt+7fv58ohtwXAlEJaKk3LeCtu7trMQEMAn4m4Jd89q2wFZYcNWrUMNVH6tWrFzpNv1i00olTaNnPyUXbU59AtH0DT58+LePGjTO1TmOxnTt3ilbdiWZTpkx5ZY\/DWO7LORAoDoFevXqJ5pxu5KwbOqeqWSdslSpVki1btoj+qabV\/rV8l\/6ntmnTJvn+++9TNR60yxICznZMjjvufQMvXrwoWVlZJjdfZ84Sl2jnfffdd+YHGwaBRBPQXVc2b94saWlpCFuiYRe8\/9dffy0fffSROfzTTz+ZKtRaV1J7ai1btjTHY9mp2+t28zx7CLRv396MGKgV3DfQqd6\/detW0Q1xC7O6deuawuFPnjwxvbxIdv78eYbZ7UmdlPNExWzQoEFmc+cmTZrIm2++adpIj83jUO3du1f03dqpU6dk2LBhoaeruOXm5prAHDlyRLKzsz1uGY8LCoHVq1dLw4YN5eHDh6J7DupuFI45n+legrp7fGGmew5OmzZNrl27ZrZiwiDgNYFatWqZvTILGsLmYST014VT4T\/SEM2SJUukefPm8ujRI7PDAAaBRBDQ92JaAk5zMScnJ+wR\/fv3D\/2o0vcV169fj9qEiRMnmvcYR48elbFjxyaiqdwTAoUS0O\/UqVOnhs5p06aNVKlShR6bl3nz\/vvvy4oVK8wjO3ToIPn5+WGP1x7c4MGDRTcjbdu2rZdN41kBIuDsG6hDjTrk6LbatWvLxo0bzaHXDYkvXrxYWrRoIfpOTk1nApctW1Z0D0OtwqO7xmMQ8JKA8+6YHpuH1J1fw9GEy5nRo03q2LFjzDPTPHSBR\/mcgHvfwEjC5R5VcO8bGMltHQLSoaBodu7cOfniiy98Tozm+4kAwpaEaE2fPl26detmXrZ36tTplRZ07txZZsyYYY737NmTdUVJiJHtjyzuvoGRuOzZs8cUJNCRB51Esn\/\/frMDvU5AycjIMJfEMgnFdub45x0BhM071qEnTZo0ST777DPTE9MeWUHTz\/QctXbt2oW91E9Cc3mkhQTiuW\/giBEjjJCtXbtWdPaj27Zt2yZVq1aNmusWosWlFCCAsCUhCAMGDDC\/ZqMNRY4cOVIGDhxo1g9lZmYmoYU80nYCXu0bOGTIENHd4dV0JEILEGAQSDQBhC3RhCPcX4cfnbVBuhPA2bNnw86aO3eumVSiXwL6ZYBBIBEE4rFvoC5Z0WHIZ8+eRVyn5h59iDRRKhF+cU8IIGxJyAH9MsjLyxMtZ6TVGHTKv9v++OMPs47twIEDMnny5CS0kEcGgUA89g10Rhd09EGH1QvO8F20aJG0atWKocggJFQK+YiwJSkYDnj9QtAV806vTddi9OjRw7RKhyuPHz+epBbyWNsJFHXfQF2m4hQT0B9jN2\/eFKfqiLI6ePBg6N2w\/lvXYM6cOdOUiYu0Vs52vviXPAIIW5LY6xeCltIqX768aYGu+dGenFNIVutIzpkzJ0mt47FBIFDUfQOdd8PKRn+MnTlzxmBatWqVNGrUyPxdh8917ZDmsU4oUdPZv927d6ekVhCSKkV8RNiSGAita7Z8+XIjaI5pD+7YsWNUcEhiXIL06KLsG6gFkceMGWPwuIVN81drnOoi7YKmIjdq1CiWrAQpqVLAV0fYLl++LH379k2BFkVugnXV\/d1uqsBphRGttadrgmLdKiRlo0XDfEcgln0DX+eU3kPfp7399ttm2v+JEycQtNdB4\/NAE7Ba2AIdWZyHAAQgEFACCFtAA4\/bEIAABGwlgLDZGln8ggAEIBBQAghbQAOP2xCAAARsJYCw2RpZ\/IIABCAQUAIIW0ADj9sQgAAEbCWAsNkaWfyCAAQgEFACCFtAA4\/bEIAABGwlgLDZGln8ggAEIBBQAghbQAOP2xCAAARsJYCw2RpZ\/IIABCAQUAIIW0ADj9sQgAAEbCWAsNkaWfyCAAQgEFACCFtAA4\/bEIAABGwlgLDZGln8ggAEIBBQAghbQAOP2xCAAARsJYCw2RpZ\/IIABCAQUAIIW0ADj9sQgAAEbCWAsNkaWfyCAAQgEFACCFtAA4\/bEIAABGwlgLDZGln8ggAEIBBQAghbQAOP2xCAAARsJYCw2RpZ\/IIABCAQUAIIW0ADj9sQgAAEbCWAsNkaWfyCAAQgEFACCFtAA4\/bEIAABGwlgLDZGln8ggAEIBBQAghbQAOP2xCAAARsJYCw2RpZ\/IIABCAQUAK+F7aaNWtK2bJl5dKlSwENIW5DAAIQgICbgK+F7a233pKdO3fKy5cvpW3btkQWAhCAAAQgIL4Wtvnz50tmZibCRiJDAAIQgECIgO+ErVevXtK0aVNp3ry5pKWlGUfosZHREIAABCDgEPCdsO3du1fKlSsXFkGEjYSGAAQgAAHfCtuIESMkPT3dtL9evXry3nvv0WMjnyEAAQhAwL9Dke7YDR48WIYNG4awkdAQgAAEIBAcYRs6dKgMGTIkLOSjR4+W48ePkwYQgAAEIGAhAd+9Yytqj80RNhUzxxA1CzMZlyAAAQj8P4HACFubNm0IOgQgAAEIBICA58JWunRpeeedd+TBgwfy999\/lwhxLO\/YnB4bwlYi1FwMAQhAwDcEPBM2FbRZs2ZJhw4dpFSpUgbQ06dP5fTp0zJu3Djz91hs8eLF0qxZM3Oq3se514sXL8yx7du3y+zZs0O3Qthioco5EIAABOwh4JmwrVy5Uho3bhwip0KkYqd28eJFycrKEkecCsO7detWqVatWtRT9u3bJzk5OQibPTmKJxCAAASKRMATYWvfvr3MmzfPNOzkyZMyfvx4efz4sQwYMECcSR0qWO6eVjQvnAXa06dPl1atWsknn3xipvuPHDnSXHLjxg25fv06wlakNOBkCEAAAvYQ8ETYVq9eLQ0bNpSHDx9Kly5d5Pnz5yGCzmd37tyRrl27Fkr2jTfeEO2ROZVGeMdmTyLiCQQgAIF4EfBE2LQCv1biLzhMqE70799fsrOzjT9aB9Ld2yropA5l6pCmTjzp3LmzOMKmQ5haDDmS8Y4tXqnCfSAAAQj4g4AnwqaCpr0tHWrUIUe31a5dWzZu3GgO6RDl4cOHo5Lr16+fTJgwQbR3d+XKFalTp46pG6n\/Pnr0qMyZM0fy8\/PDrkfY\/JGItBICEIBAvAgkXNgqV64su3fvjipczvCinqDv4X777beovk2dOlV69OgR9XPtyX3++edhywioPBKvVOE+EIAABPxBIOHC1rp1a1m4cKGhoaJz4cKFV8j8+eef5tjy5ctlzZo1UcktW7ZMPvjgA\/N5bm6u7Nixw\/TQBg0aJC1atDDHtSfXp0+f0D2oPOKPRKSVEIAABOJFIOHCpkKkgqT25ZdfytmzZ8PaXqFCBcnLyzPHIg1Vuk\/u3r27WcN26NAh+e9\/\/xt2Hx2G7NixozmmwqYCp8ZQZLxShftAAAIQ8AeBhAubThrRySNqkyZNkoMHD4aR0SokGzZsMMde947NfWHBCiY1atQIDWPqztq\/\/vorwuaPHKSVEIAABOJKIOHCpq11hhrXr18vP\/74Y5gDOsX\/m2++Mce6desmt27diupglSpVTKURXYBdsILJmTNn5N\/\/\/vcrPT96bHHNF24GAQhAIOUJeCJs2nvKyMgwi6c\/\/fTTMChORRIVNBW2aObu2bnPcSqY6No2p7xW37595fLly\/TYUj79aCAEIACB+BPwRNj03drw4cNN690TRHSx9rfffmsE6ZdffpEFCxaYc3SmpM6QLFOmjPz8889y4MABc9xZNqBipjUmx44dayqYjBkzxkxM0fs8efJEOnXqFCJFjy3+ScMd409A81S3U2JLpfiz5Y6dBlYbAAAH7ElEQVTxI+CXPPVE2HStmQpU9erVDeF79+4ZQXJqPmrvSntZjrnfy+n7t0WLFpmP3HUitWjy1atXjZjpejan7uT9+\/flP\/\/5D8IWv1zmTh4QWLp0qXmKe99ADx7LIyBQJAL6WmnVqlWmUEYqmyfCpgAqVqwo69atC4mbA+Wff\/4xNSPv3r0b4uRe++YWNqeCiZbk0t6c21ToVEDV3BVM6LGlcvrRNocAwkYu+IEAwhYlSunp6fLhhx+antauXbvk9u3bMcfTGYrUIcv\/\/e9\/Ur9+fbMY+8SJE2bbm0gVTBC2mPFyYhIJIGxJhM+jYyaAsMWMKrYTi1vBJFLlkdieyFkQgAAEIFCQAEORccyJ4lYwcTYljWNTuBUEIACBwBLwwwQnz96xlTQL4lnBpKRt4XoIQAACEEhdAr4RtkRVMEnd0NAyCEAAAhAoDgHfCJs6F68KJsUBxTUQgAAEIOAPAr4StnhUMPFHWGglBCAAAQgUl4CvhK2oFUyKC4XrIAABCEDAvwR8JWxFrWDi37DQcghAAAIQKC4BXwmbOlmUCiYqhO+++64pvaWltjAIeE2g4PZKXj+f50EgngT8ks++EzYnSIVVMNG92bSIcr169UIxffDggal04hRajmewuRcEChLQL4BZs2a9sr2SFu8eN26cqZQTi2kZOS1OEM2mTJnyyh6HsdyXcyBQHAJarlBzTjdy1g2dU9V8K2zRgFaqVEm2bNki+qea7gSg5bucLW02bdok33\/\/farGg3ZZQsDZjslxx9leSf998eJFycrKMrn5OnNmAkc777vvvjM\/2DAIJJqA7rqyefNmSUtLQ9gSDbvg\/b\/++mv56KOPzOGffvrJVKHWgsnaU2vZsqU5XpSdur1uP8\/zP4H27dubEQO1kydPmnzT3Sy02LdTvV93qpg9e3ahztatW9cUDtetmLSXF8nOnz\/PMLv\/UyZlPVAxGzRokDRo0ECaNGkib775pmkrPTaPQ7Z3715T5f\/UqVMybNiw0NNV3HJzc01gjhw5ItnZ2R63jMcFhcDq1aulYcOG8vDhQ9E9B3U3Csecz+7cuSO6e3xh9vHHH8u0adPk2rVr0rt376Dgw88UIlCrVi2zV2ZBQ9g8DJL+utAdANQiDdEsWbJEmjdvLo8ePTI7DGAQSAQBZ3slzcWcnJywR\/Tv3z\/0o8q9vVKkdkycONG8xzh69KjZVBeDgNcE9Dt16tSpoce2adNGqlSpQo\/Ny0C8\/\/77smLFCvPIDh06SH5+ftjjtQc3ePBgefnypbRt29bLpvGsABFwtlfSoUYdcnRb7dq1I26vFAnP4sWLpUWLFuadnJpOmCpbtqzoHoZarGD9+vUBooqrqUDAeXdMj83DaDi\/hqMJlzOjR5vUsWPHmGemeegCj\/I5geJurxTJbR0C0qGgaHbu3Dn54osvfE6M5vuJAMKWhGhNnz5dunXrZl62d+rU6ZUWdO7cWWbMmGGO9+zZ02xSikEgngSKu71SpDbs2bPHrNvUkQedRLJ\/\/36zA71OQMnIyDCXxDIJJZ7+ca9gE0DYkhD\/SZMmyWeffWZ6YtojK2j6mZ6j1q5du7CX+kloLo+0kEA8t1caMWKEEbK1a9eKzn5027Zt26Rq1apRc91CtLiUAgQQtiQEwZlOHW0ocuTIkTJw4ECzfigzMzMJLeSRthPwanulIUOGiO4Or6YjEVqAAINAogkgbIkmHOH+OvzorA3Sgslnz54NO2vu3LlmUol+CeiXAQaBRBCIx\/ZKumRFhyGfPXsWcZ2ae\/Qh0kSpRPjFPSGAsCUhB\/TLIC8vT7SckVZj0Cn\/bvvjjz\/MOrYDBw7I5MmTk9BCHhkEAvHYXskZXdDRBx1WLzjDd9GiRdKqVSuGIoOQUCnkI8KWpGA44PULQVfMO702XYvRo0cP0yp9+X78+PEktZDH2k6gqNsr6TIVp5iA\/hi7efOmOFVHlNXBgwdD74b137oGc+bMmaZMXKS1crbzxb\/kEUDYksRevxC0lFb58uVNC3TNj\/bknEKyWkdyzpw5SWodjw0CgaJur+QutaU\/xs6cOWMwrVq1Sho1amT+rsPnunZI81gnlKjp7N\/u3btTUisISZUiPiJsSQyE1jVbvny5ETTHtAd37NgxKjgkMS5BenRRtlfSgshjxowxeNzCpvmrNU51kXZBU5EbNWoUS1aClFQp4KsjbJcvX5a+ffumQIsiN8G66v5uN1XgtMKI1trTNUGxbhWSstGiYb4jUNj2SrE6o\/fQ92lvv\/22mfZ\/4sQJBC1WeJwXSAJWC1sgI4rTEIAABAJOAGELeALgPgQgAAHbCCBstkUUfyAAAQgEnADCFvAEwH0IQAACthFA2GyLKP5AAAIQCDgBhC3gCYD7EIAABGwjgLDZFlH8gQAEIBBwAghbwBMA9yEAAQjYRgBhsy2i+AMBCEAg4AQQtoAnAO5DAAIQsI0AwmZbRPEHAhCAQMAJIGwBTwDchwAEIGAbAYTNtojiDwQgAIGAE0DYAp4AuA8BCEDANgIIm20RxR8IQAACASeAsAU8AXAfAhCAgG0EEDbbIoo\/EIAABAJOAGELeALgPgQgAAHbCCBstkUUfyAAAQgEnADCFvAEwH0IQAACthFA2GyLKP5AAAIQCDgBhC3gCYD7EIAABGwjgLDZFlH8gQAEIBBwAghbwBMA9yEAAQjYRgBhsy2i+AMBCEAg4AT+D0l0IbDrI0n1AAAAAElFTkSuQmCC","height":132,"width":219}}
%---
