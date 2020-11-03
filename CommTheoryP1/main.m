clear;clc;close all;

%loading audio and basic variables
[m_o, fs_o] = audioread('Dame_Da_NE.mp3'); 
m_o = transpose(m_o(:,1)/max(abs(m_o(:,1))));
T_o = length(m_o)/fs_o;
t_o = linspace(0,T_o,length(m_o));
scale = 40; %just a scaling factor I chose 800/44*2 rounded up
m = interp(m_o,scale);
fs = fs_o * scale;
N_o = fs_o*T_o;

w_o = linspace(-pi,pi,N_o);
f_o = w_o*fs_o/(2*pi);

T = length(m)/fs;
N = fs*T;
t = linspace(0,T,N); 
w = linspace(-pi,pi,N);
f = w*fs/(2*pi);

fc = 800000;
c = cos(2*pi*fc*t);
s = sin(2*pi*fc*t);

%plot original
figure
plot(t_o, m_o)
title("original norm'd signal in time");
xlabel("time")
ylabel("magnitude")
M_O = fftshift(fft(m_o)/fs_o);

figure
semilogy(f_o, abs(M_O))
title("original norm'd signal in freq");
xlabel("freq")
ylabel("magnitude (logged)")

%SSB
SSB = SSB_Mod(m,s,c,f);
f_SSB = fftshift(fft(SSB));
demod_SSB = SSB_demod(SSB, c, 10000, fs);
DEMOD_SSB = fftshift(fft(demod_SSB)/fs);
figure
subplot(2,2,1)
plot(t,SSB)
title('SSB AM in time');
xlabel("time")
ylabel("magnitude")
subplot(2,2,3)
semilogy(f,abs(f_SSB))
xlabel("frequency")
ylabel("magnitude(logged)")
title('SSB AM in freq');    
subplot(2,2,2)
plot(t_o,demod_SSB)
xlabel("time")
ylabel("magnitude")
title("SSB demod in time")
subplot(2,2,4)
semilogy(f_o, abs(DEMOD_SSB))
xlabel("frequency")
ylabel("magnitude(logged)")
title("SSb demod in freq")

%AM
conv_am = conven_AM(m, fc, t,1);
f_conv = fftshift(fft(conv_am));
demod_conv = conv_demod(conv_am,fs,10000);
demod_AM = fftshift(fft(demod_conv)/fs);

figure
subplot(2,2,1)
plot(t,conv_am)
xlabel("time")
ylabel("magnitude")
title('conv in time')
subplot(2,2,3)
semilogy(f,abs(f_conv))
title('conv in freq')
xlabel("frequency")
ylabel("magnitude(logged)")
subplot(2,2,2)
plot(t_o, demod_conv)
title("demod AM in time")
xlabel("time")
ylabel("magnitude")
subplot(2,2,4)
semilogy(f_o, abs(demod_AM))
xlabel("frequency")
ylabel("magnitude(logged)")
title("am_demod in freq")

%FM
fm_mod = FM_mod(m,t,fc,fs,14000);
FM_mod = fftshift(fft(fm_mod)/fs);
figure
plot(f, abs(FM_mod))
xlabel("frequency")
f_fm_mod = fftshift(fft(fm_mod));
demod_fm = FM_demod(fm_mod, fs, 10000, f);
demod_FM = fftshift(fft(demod_fm)/fs);
figure
subplot(2,2,1)
plot(t, fm_mod);
xlabel("time")
ylabel("magnitude")
title("fm_mod in time")
subplot(2,2,3)
semilogy(f, abs(f_fm_mod));
title("fm mod in freq")
xlabel("frequency")
ylabel("magnitude(logged)")
subplot(2,2,2)
plot(t_o, demod_fm)
title("demod fm in time")
xlabel("time")
ylabel("magnitude")
subplot(2,2,4)
semilogy(f_o,abs(demod_FM))
xlabel("frequency")
ylabel("magnitude(logged)")
title("fm demod in freq")

%PM
pm_mod = PM_mod(m,T,fc, 2);
PM_mod = fftshift(fft(pm_mod)/fs);
figure
plot(f, abs(PM_mod))
xlabel("frequency")

f_pm_mod = fftshift(fft(pm_mod));
demod_pm = PM_demod(pm_mod, s, 10000,fs);
demod_PM = fftshift(fft(demod_pm)/fs);
figure
subplot(2,2,1)
plot(t, pm_mod)
xlabel("time")
ylabel("magnitude")
title("pm_mod in time")
subplot(2,2,3)
semilogy(f,abs(f_pm_mod))
xlabel("frequency")
ylabel("magnitude(logged)")
title("pm mod in freq")
subplot(2,2,2)
plot(t_o, demod_pm)
xlabel("time")
ylabel("magnitude")
title("pm demod in time")
subplot(2,2,4)
semilogy(f_o,abs(demod_PM))
xlabel("frequency")
ylabel("magnitude(logged)")
title("pm demod in freq")

%power (matching)
Orgpower = rms(m_o)^2;
fprintf('Power of original Signal: %.4f\n',Orgpower);
Convpower = rms(demod_conv)^2;
fprintf('Power of Conv Signal: %.4f\n',Convpower);
SSBpower = rms(demod_SSB)^2;
fprintf('Power of SSB Signal: %.4f\n',SSBpower);
FMpower = rms(demod_fm)^2;
fprintf('Power of FM Signal: %.4f\n',FMpower);
PMpower = rms(demod_pm)^2;
fprintf('Power of PM Signal: %.4f\n',PMpower);

%P(clean demod_org)/p(noise -clean) P(mod, + noise - clean)
