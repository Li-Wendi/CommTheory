clear;clc;close all;

[mt,fs] = audioread('Dame_Da_NE.mp3'); 
mt = transpose(mt(:,1));
mt = mt./max(abs(mt));
N = length(mt);
T = fs/N;

up_N = N*40;
up_fs = fs*40;
up_t = linspace(0,T,up_N);

fc = 800000;
f_cutoff = fs/2;
w = linspace(-pi,pi,N);
up_w = linspace(-pi,pi, up_N);

f = up_w*up_fs/(2*pi);

m_up = interp1(t, mt, up_t);  
s = sin(2*pi*fc*up_t);
c = cos(2*pi*fc*up_t);

conv_am = conven_AM(m_up, fc, up_t);
conv_AM = fft(conv_am);
conv_AM_s = fftshift(conv_AM/up_fs);
   
%y = cumsum(m_up/up_fs)
%result = 1*cos(2*pi*fc*up_t + 2*pi*1*cumsum(m_up/up_fs))/fs;

figure;
plot (f,abs(conv_AM_s));
