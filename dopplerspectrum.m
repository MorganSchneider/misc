
%close all

iqh_mod = reshape([squeeze(iqh(10,19,:)); squeeze(iqh(14,31,:))], [1 1 200 1]);


figure(1);

%for xi = 1:21

xi = 1;
yi = 25;

%iqh_mod = iqh(10,19,:,1) + iqh(14,31,:,1);

fs = 2*dat.params.va/dat.params.lambda;
[ss,ff] = periodogram(squeeze(iqh(xi,yi,:)),[],linspace(-fs,fs,size(iqh,3)),dat.params.prf);
vv = ff/2*dat.params.lambda;
vvv = dat.params.va / pi * angle(iqh(:,:,2:end) .* conj(iqh(:,:,1:end-1)));

subplot(2,1,1)
histogram(squeeze(vvv(xi,yi,:)),30)
xlim([-dat.params.va,dat.params.va])
xlabel('vr')
title(['Fake Doppler spectrum? x=' num2str(xi) ' y=' num2str(yi) ' \theta_e=0.5'])
subplot(2,1,2)
semilogy(vv,ss)
xlim([-dat.params.va,dat.params.va])
xlabel('vr')
title(['Actual periodogram x=' num2str(xi) ' y=' num2str(yi) ' \theta_e=0.5'])
shg