clear;clc;close all;

%loading audio and basic variables
[m_o, fs_o] = audioread('Dame_Da_NE.mp3'); 
m_o = transpose(m_o(:,1)/max(abs(m_o(:,1))));
T_o = length(m_o)/fs_o;
t_o = linspace(0,T_o,length(m_o));
N_o = fs_o*T_o;
w_o = linspace(-pi,pi,N_o);
f_o = w_o*fs_o/(2*pi);


scale = 40; %just a scaling factor I chose 800/44*2 rounded up
m = interp(m_o,scale);
fs = fs_o * scale;

T = length(m)/fs;
N = fs*T;
t = linspace(0,T,N); 
w = linspace(-pi,pi,N);
f = w*fs/(2*pi);


Ac = 1;
fc = 800000;
c = cos(2*pi*fc*t);
s = sin(2*pi*fc*t);

%SSB
var = [0.01, 0.05, 0.1];
M = fftshift(fft(m_o)/fs_o);
ssb_mod = SSB_Mod(m,s,c,f);
f_SSB = fftshift(fft(ssb_mod)/fs);
demod_ssb = SSB_demod(ssb_mod, c, 10000, fs);
noiseless_ssb = (rms(demod_ssb)).^2;

k = 0;
for i = 1:3
    v = var(i);
    noise = sqrt(var(i)) * randn(1, length(m));
    noise_ssb = ssb_mod + noise;
    MOD_Noise_SSB = fftshift(fft(noise_ssb)/fs);
    
    demod_ssb_noise = SSB_demod(noise_ssb,c,10000, fs);
    demod_ssb_noise = lowpass(demod_ssb_noise, 5000, fs_o);
    demod_SSB_noise = fftshift(fft(demod_ssb_noise)/fs_o);
    
    if (k < i)
        figure
    end
    
    subplot (2,2,1)
    semilogy(f_o, abs(M));
    title("original signal in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    subplot (2,2,2)
    semilogy(f, abs(MOD_Noise_SSB));
    title("SSB mod'd with noise variance of "+ var(i));
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    subplot (2,2,4)
    semilogy(f_o, abs(demod_SSB_noise));
    title("SSB demod'd with noise variance of "+ var(i));
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    subplot (2,2,3)
    semilogy(f, abs(f_SSB));
    title("SSB_mod'd in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    noisy = rms(demod_ssb_noise - demod_ssb).^2;
    snrconv = 10* log10(noiseless_ssb/noisy);
    theoretical_ssb = 10* log10(noiseless_ssb * 1.9^2 / (var(i) *2 * 7000/fs ));
    disp("SNR for SSB with variance " + var(i) + ": " + snrconv);
    disp("Theo_SNR for SSB with variance " + var(i) + ": " + theoretical_ssb);
end

k = 0;
%AM
am_mod = conven_AM(m, fc, t,2);
f_conv = fftshift(fft(am_mod)/fs);
demod_am = conv_demod(am_mod,fs,10000);
noiseless_am = (rms(demod_am)).^2;

for i = 1:3
    v = var(i);
    noise = sqrt(var(i)) * randn(1, length(m));
    noise_am = am_mod + noise;
    MOD_Noise_am =  fftshift(fft(noise_am)/fs);
    
    demod_am_noise = conv_demod(noise_am, fs, 10000);
    demod_am_noise = lowpass(demod_am_noise, 5000, fs_o);
    demod_AM_noise = fftshift(fft(demod_am_noise)/fs_o);
    if (k < i)
        figure
    end
    
    subplot (2,2,1)
    semilogy(f_o, abs(M));
    title("original signal in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    subplot (2,2,2)
    semilogy(f, abs(MOD_Noise_am));
    title("AM mod'd with noise variance of "+ var(i));
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    subplot (2,2,4)
    semilogy(f_o, abs(demod_AM_noise));
    title("AM demod'd with noise variance of "+ var(i));
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    subplot (2,2,3)
    semilogy(f, abs(f_conv));
    title("AM_mod'd in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    noisy = rms(demod_am_noise - demod_am).^2;
    snrconv = 10* log10(noiseless_am/noisy);
    theoretical_am = 10*log10(2^2*.6^2*noiseless_am/(2*2*var(i)*7000/fs));
    
    disp("SNR for AM with variance " + var(i) + ": " + snrconv);
    disp("Theo_SNR for AM with variance " + var(i) + ": " + theoretical_am);
    k = i;
end

k = 0;
%fm mod part with noise
fm_mod = FM_mod(m,t,fc,fs,14000);
demod_fm = FM_demod(fm_mod, fs, 10000, f);
MOD_fm = fftshift(fft(fm_mod/fs));
noiseless_fm = (rms(demod_fm)).^2;

for i = 1:3
    v = var(i);
    noise = sqrt(var(i)) * randn(1, length(m));
    noise_fm = fm_mod + noise;
    MOD_Noise_fm = fftshift(fft(noise_fm)/fs);
    
    demod_fm_noise = FM_demod(noise_fm, fs, 10000, f);
    demod_fm_noise = lowpass(demod_fm_noise, 5000, fs_o);
    demod_FM_noise = fftshift(fft(demod_fm_noise)/fs_o);
    
    if (k < i)
        figure
    end
    
    subplot (2,2,1)
    semilogy(f_o, abs(M));
    title("original signal in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    subplot (2,2,2)
    semilogy(f, abs(MOD_Noise_fm));
    xlabel("frequency")
    ylabel("magnitude(logged)")
    title("FM mod'd with noise variance of "+ var(i));
    
    subplot (2,2,4)
    semilogy(f_o, abs(demod_FM_noise));
    title("FM demod'd with noise variance of "+ var(i));
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    subplot (2,2,3)
    semilogy(f, abs(MOD_fm));
    title("FM_mod'd in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    noisy = rms(demod_fm_noise - demod_fm).^2;
    snrconv = 10* log10(noiseless_fm/noisy);
    theoretical_fm = 10 * log10( 3 * 14000^2 * ((11/5)^2) * (noiseless_fm) / (2* (7000^2)*( 2 * (7000/fs * (var(i))))) );
    
    disp("SNR for FM with variance " + var(i) + ": " + snrconv);
    disp("Theo_SNR for FM with variance " + var(i) + ": " + theoretical_fm);
    k = i;
end

%PM
pm_mod = PM_mod(m,T,fc, 2);
f_pm_mod = fftshift(fft(pm_mod));

demod_pm = PM_demod(pm_mod, s, 10000,fs);
MOD_pm = fftshift(fft(pm_mod/fs));
noiseless_pm = (rms(demod_pm)).^2;

k = 0;
for i = 1:3
    v = var(i);
    noise = sqrt(var(i)) * randn(1, length(m));
    noise_pm = pm_mod + noise;
    MOD_Noise_pm = fftshift(fft(noise_pm)/fs);
    
    demod_pm_noise = PM_demod(noise_pm,s, 10000, fs);
    demod_pm_noise = lowpass(demod_pm_noise, 5000, fs_o);
    demod_PM_noise = fftshift(fft(demod_pm_noise)/fs_o);
    
    if (k < i)
        figure
    
    end
    subplot (2,2,1)
    semilogy(f_o, abs(M));
    title("original signal in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    subplot (2,2,2)
    semilogy(f, abs(MOD_Noise_pm));
    title("PM mod'd with noise variance of "+ var(i));
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    subplot (2,2,4)
    semilogy(f_o, abs(demod_PM_noise));
    title("PM demod'd with noise variance of "+ var(i));
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    subplot (2,2,3)
    semilogy(f, abs(MOD_pm));
    title("PM_mod'd in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    
    noisy = rms(demod_pm_noise - demod_pm).^2;
    snrconv = 10* log10(noiseless_pm/noisy);
    theoretical_pm = 10* log10(2^2 * 1.05^2 * noiseless_pm /(2*var(i)*7000/fs));
   
    disp("SNR for PM with variance " + var(i) + ": " + snrconv);
    disp("Theo_SNR for PM with variance " + var(i) + ": " + theoretical_pm);
    k = i;
end

k= 0;
%we will be choosing the best variance (0.01)
index = [0.5, 1, 2];
variance = 0.01;
noise = sqrt(variance) * randn(1, length(m));

figure
for i = 1:3
    a = 2*index(i);
    str = strcat('a = ' , int2str(a));
    pm_mod = PM_mod(m,T,fc, index(i)*2);
    MOD_pm = fftshift(fft(pm_mod)/fs);
    demod_pm = PM_demod(pm_mod, s, 10000,fs); 
    
    pm_noise = pm_mod + noise;
    MOD_Noise_pm = fftshift(fft(pm_noise)/fs);
    
    demod_pm_noise = PM_demod(pm_noise,s, 10000, fs);
    demod_pm_noise = lowpass(demod_pm_noise, 5000, fs_o);
    demod_PM_noise = fftshift(fft(demod_pm_noise)/fs_o);
    
    subplot (2,2,1)
    semilogy(f_o, abs(M),'DisplayName',str);
    title("original signal in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    legend
    hold on
    
    subplot (2,2,2)
    semilogy(f, abs(MOD_Noise_pm),'DisplayName',str);    
    title("PM mod'd with noise");
    xlabel("frequency")
    ylabel("magnitude(logged)")
    legend
    hold on
    
    subplot (2,2,4)
    semilogy(f_o, abs(demod_PM_noise),'DisplayName',str);
    title("PM demod'd with noise");
    xlabel("frequency")
    ylabel("magnitude(logged)")
    legend
    hold on
    
    subplot (2,2,3)
    semilogy(f, abs(MOD_pm),'DisplayName',str);
    title("PM_mod'd in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    legend
    hold on
    
    noiseless_pm = (rms(demod_pm)).^2;
    noisy_pm = rms(demod_pm_noise - demod_pm).^2;
    snrconv_pm = 10* log10(noiseless_pm/noisy_pm);
    theoretical_pm = 10* log10((index(i)*2)^2 * 1.05^2 * noiseless_pm /(2*variance*7000/fs));
    
    disp("SNR for PM with index " + index(i)*2 + ": " + snrconv_pm);
    disp("Theo_SNR for PM with index " + index(i)*2 + ": " + theoretical_pm);
end
hold off

figure
for i = 1:3
    fm_mod = FM_mod(m,t,fc,fs,index(i)*14000);
    MOD_fm = fftshift(fft(fm_mod)/fs);
    demod_fm = FM_demod(fm_mod, fs, 10000, f);
    a = 2*index(i);
    str = strcat('a = ' , int2str(a));
    fm_noise = fm_mod + noise;
    
    MOD_Noise_fm = fftshift(fft(fm_noise)/fs);
    
    demod_fm_noise = FM_demod(fm_noise,fs, 10000, f);
    demod_fm_noise = lowpass(demod_fm_noise, 5000, fs_o);
    demod_FMnoise = fftshift(fft(demod_fm_noise)/fs_o);
    
    subplot (2,2,1)
    semilogy(f_o, abs(M),'DisplayName',str);
    title("original signal in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    legend
    hold on
    
    subplot (2,2,2)
    semilogy(f, abs(MOD_Noise_fm),'DisplayName',str);
    title("FM mod'd with Noise");
    xlabel("frequency")
    ylabel("magnitude(logged)")
    legend
    hold on
    
    subplot (2,2,4)
    semilogy(f_o, abs(demod_FM_noise),'DisplayName',str);
    title("FM demod'd with Noise");
    xlabel("frequency")
    ylabel("magnitude(logged)")
    legend
    hold on
    
    subplot (2,2,3)
    semilogy(f, abs(MOD_fm),'DisplayName',str);
    title("FM_mod'd in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    legend
    hold on
    
    noiseless_fm = (rms(demod_fm)).^2;
    noisy_fm = rms(demod_fm_noise - demod_fm).^2;
    snrconv_fm = 10* log10(noiseless_fm/noisy_fm);
    theoretical_fm = 10 * log10( 3 * (index(i)*14000)^2 * ((11/5)^2) * (noiseless_fm) / (2* (7000^2)*( 2 * (7000/fs * (variance)))) );
   
    disp("SNR for FM with index of " + index(i)*2 + ": " + snrconv_fm);
    disp("Theo_SNR for fM with index " + index(i)*2 + ": " + theoretical_fm);
end
hold off

figure
for i = 1:3

    a = 2*index(i);
    str = strcat('a = ' , int2str(a));
    am_mod = conven_AM(m, fc, t, 2* index(i));
    MOD_am = fftshift(fft(am_mod)/fs);
    demod_am = conv_demod(am_mod,fs,10000);
    
    am_noise = am_mod + noise;
    MOD_Noise_am = fftshift(fft(am_noise)/fs);
    
    
    demod_am_noise = conv_demod(am_noise,fs, 10000);
    demod_am_noise = lowpass(demod_am_noise, 5000, fs_o);
    demod_AM_noise = fftshift(fft(demod_am_noise)/fs_o);

    subplot (2,2,1)
    semilogy(f_o, abs(M),'DisplayName',str);
    title("original signal in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    legend
    hold on
    subplot (2,2,2)
    semilogy(f, abs(MOD_Noise_am),'DisplayName',str);
    title("AM mod'd with Noise");
    xlabel("frequency")
    ylabel("magnitude(logged)")
    legend
    hold on
    subplot (2,2,4)
    semilogy(f_o, abs(demod_AM_noise),'DisplayName',str);
    title("AM demod'd with Noise");
    xlabel("frequency")
    ylabel("magnitude(logged)")
    legend
    hold on
    
    subplot (2,2,3)
    semilogy(f, abs(MOD_am),'DisplayName',str);
    title("AM_mod'd in freq")
    xlabel("frequency")
    ylabel("magnitude(logged)")
    legend
    hold on
    
    noiseless_am = (rms(demod_am)).^2;
    noisy_am = rms(demod_am_noise - demod_am).^2;
    
    snrconv_am = 10* log10(noiseless_am/noisy_am);
    theoretical_am = 10* log10((index(i)*2)^2*.6^2*noiseless_am/(2*2*variance*7000/fs));
    
    disp("SNR for AM with index " + index(i)*2 + ": " + snrconv_am);
    disp("Theo_SNR for AM with index " + index(i)*2 + ": " + theoretical_am);   

end
hold off
% I do not want to call the sound function because it will just be a mess so you may call the sound of the demod'd function by 
%enter sound(demod_"type"_noise, fs_o) and the sound should appear.      