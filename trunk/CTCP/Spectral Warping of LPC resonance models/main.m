% Load a speech waveform
[d,sr] = wavread('we_were_away.wav');

% Fit the original LPC model (high-order)
[a,g,e] = lpcfit(d,20);

% Warp the poles up a little
% (warppoles modifies every frame - rows of a - at the same time)
alpha = -0.2;
[bhat, ahat]  = warppoles(a, alpha);

% Resynthesize with the new LPC
% (fortunately, bhat is the same for all frames)
dw = filter(bhat(1,:), 1, lpcsynth(ahat, g, e));
