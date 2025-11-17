# hilbert-analytic-signal-intuition [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=arkaddas/hilbert-analytic-signal-intuition&file=README.md)

<img width="902/6" height="687/6" alt="image" src="https://github.com/user-attachments/assets/48f079d2-e936-4fcf-8910-2128d360a3c4" />
<img width="936/6" height="674/6" alt="image" src="https://github.com/user-attachments/assets/f0971291-af8b-436b-90df-8af900772013" />
<img width="958/6" height="684/6" alt="image" src="https://github.com/user-attachments/assets/b3257016-533b-4310-8608-c9d53b02bb35" />
<img width="885/6" height="677/6" alt="image" src="https://github.com/user-attachments/assets/f2633a3f-199a-40d0-9714-56e1a8fc6eb5" />



Minimal MATLAB &amp; Python demos for intuition behind analytic signals and the Hilbert transform.

Minimal MATLAB demos to build intuition for the Hilbert transform and analytic signals: FFT, cancellation of negative frequencies, envelopes, instantaneous phase, and positive-frequency representation.


Usage:
matlab -r "demo_hilbert_analytic"

Short Elxplanation
The Hilbert transform allows us to create an analytic signal whose Fourier transform contains only positive frequencies.

From this analytic signal:
	•	Magnitude → amplitude envelope
	•	Angle → instantaneous phase
	•	Derivative of angle → instantaneous frequency

This repo contains minimal, highly commented demos showing:
✔ Real → analytic signal conversion
✔ Envelopes
✔ Phase
✔ Chirps / AM examples
✔ Spectrum before vs after
