clear all
close all

L = 2;
D = 256;
t = linspace(0, L, D+1);
t = t(1:end-1);

q = 2*sin(8/L*2*pi*t) - 1*sin(16/L*2*pi*t);
Q = fft(q);

E = [-100 800];
grid_spacing = 0.01;
R = 1024;

[floq_det, al21]      = mex_fnft_kdvp_floquet(q, [0 L], E, 'L', R);
[main_spec, aux_spec] = mex_fnft_kdvp(q, [0 L], E, 'grid_spacing', grid_spacing, 'keep_degenerate');

ampmodfreqs = mex_fnft_kdvp_ampmodfreq(main_spec);
A = ampmodfreqs(1:3:end)    % amplitudes
m = ampmodfreqs(2:3:end)    % moduli
F = ampmodfreqs(3:3:end)  

dt = t(2)-t(1);
f = fftfreq(D, dt);
figure
stem(F, A, 'o')
hold on 
plot(f(1:D/2), abs(Q(1:D/2))/D, 'x')
hold off
grid on
xlabel('f')
ylabel('A')
legend('NFT', 'FFT')
xlim([0 20])

ix = abs(floq_det) > 1;
floq_scaled = real(floq_det);
floq_scaled(ix) = sign(real(floq_det(ix))).*(1+log(abs(floq_det(ix))));

ix = abs(al21) > 1;
al21_scaled = real(al21);
al21_scaled(ix) = sign(real(al21(ix))).*(1+log(abs(al21(ix))));

Es = linspace(E(1), E(2), R);

figure
subplot(2,1,1)
plot(t, q);

subplot(2,1,2)
plot(Es, floq_scaled, Es, al21_scaled)
xlim(E);
hold on
plot(main_spec(1:2:end), main_spec(2:2:end), 'ro')
plot(aux_spec, 0*aux_spec, 'rs')
plot([E(1) E(2)], [1 1], '-k')
plot([E(1) E(2)], [-1 -1], '-k')
hold off
grid on
xlabel('E');
legend('Floquet det. (scaled)', '\alpha_{21} (scaled)', 'Main spectrum', 'Aux. spectrum');
