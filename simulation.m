

% radio communication system simulation script
%
clc; clear variables; % clear all variables
addpath('subfunctions'); % add directory "subfunctions" to path
PLOTS_ON = true; 
% global simulation parameters
ebN0dB = 0:3:30; % SNR (per bit) in dB
nBitsPerLoop = 10000000;%10e3; % simulate nBitsPerLoop bits per simlation loop
%nMinErr = 100; % simulate at least nMinErr bit errors...
%nMaxBits = 100 * nBitsPerLoop; % ... or stop after nMaxBits
const = [-1-1j, 1-1j, -1+1j, 1+1j]; % constellation of the
% modulation format here: QPSK wiht Gray mapping)
nSamp = 100;
snrdB = 15;

% here starts the main code.

%% Bits plotting generateBits
bits = generateBits(nBitsPerLoop);
if PLOTS_ON==true
    figure(1)
    stem( bits, 'linewidth',2)
   axis([0 70 0 1])
end


%% Mapper function

ydataSymb = mapper(bits, const);
figure(2);
plot(real(ydataSymb),imag(ydataSymb), 'ored','linewidth',3);
xlim([-2 2]);
ylim([-2 2]);
line(xlim, [0 0], 'color','k','linewidth',1)
line([0 0],ylim,'color','k','linewidth',1)
%  axis ([-5 5 -5 5]);
hold on
fplot([xlim,ylim])  %ezplot  fplot fimplicit
grid on
title('QPSK constellation');
xlabel('real part');
ylabel('imaginary part');

%% Rayleigh Channel

KanalRayleigh = radioFadingChannel(nSamp);
figure(3);
plot(20*log10(abs(KanalRayleigh))); %For converting amplitude (Amp) in dB
grid on;
xlabel('time in s');
ylabel('Amplituden in dB');
%axis ([-5 5 -5 5]);
title('Rayleigh Fading');


%% SetSNR

x = ydataSymb;
ynSnr = setSNR(x, snrdB);

figure(4);
% plot(ynSnr,'ogreen'); % plot constellation with noise
% xlim([-2 2]);
% ylim([-2 2]);
% line(xlim, [0 0], 'color','k','linewidth',1)
% line([0 0],ylim,'color','k','linewidth',1)
% %  axis ([-5 5 -5 5]);
% hold on
% fplot([xlim,ylim])
% xlabel('real'); ylabel('imag');
% title('QPSK constellation with noise');

% Plot the noised signal and the original signal
plot(real(ynSnr),imag(ynSnr),'ko','MarkerFaceColor',[0 0 0],'MarkerSize',1);
axis([-2 2 -2 2]);
title(strcat('SNR=', num2str(snrdB), ' dB'));
hold on;
plot(real(x),imag(x),'ro','MarkerFaceColor',[1 0 0],'MarkerSize',8);
axis([-2 2 -2 2]);
hold off;

xlim([-2 2]);
ylim([-2 2]);
line(xlim, [0 0], 'color','k','linewidth',1)
line([0 0],ylim,'color','k','linewidth',1)
hold on
fplot([xlim,ylim])
xlabel('real'); ylabel('imag');
%title('QPSK constellation with noise');



%% Decision
%Receiver
yEntsch = decision(x, const);
figure (5);
plot(real(yEntsch),imag(yEntsch), 'o','MarkerFaceColor','b');
xlim([-2 2]);
ylim([-2 2]);
line(xlim, [0 0], 'color','k','linewidth',1)
line([0 0],ylim,'color','k','linewidth',1)
%  axis ([-5 5 -5 5]);
hold on
fplot([xlim,ylim])
title('demodulated symbols after hard decisions');
xlabel('real part');
ylabel('imaginary part');

%  figure(6)
ydemapping = demapper(x, const);
%  stem( ydemapping, 'linewidth',2)
%  axis([0 70 0 1])
%  title('Recieved BitStream');
% %  stem(ydemapping(1:200),'filled');


%% Impletation from BER
% [nError, indx, ~] = countErrors(x, bits);
SNRProSymbol = 10 .* log10 (log2(length(const)) .* 10 .^ (ebN0dB ./ 10));
berr = zeros(1, length(SNRProSymbol));
berr_ch = zeros(1, length(SNRProSymbol));



%% Impltementation from BER in Rayleigh Channel


figure(7), hold on, clf
for i = 1:length(SNRProSymbol)

        y = generateBits(nBitsPerLoop);
        yc = mapper(y, const);
        n = radioFadingChannel(length(yc));

%         % awgn
        y_noise = setSNR(yc, SNRProSymbol(i));
        [idx] = decision(y_noise, const);
        y_demap = demapper(idx, const);
        [~, ~, b_err] = countErrors(y, y_demap);
        berr(i) = b_err;

%         % rayleigh
        y_ray = yc .* n;
        y_ray_noise = setSNR(y_ray, SNRProSymbol(i));
        y_unray = y_ray_noise ./ n;

        [idx] = decision(y_unray, const);
        y_unray_demap = demapper(idx, const);
        [~,~, b_err] = countErrors(y, y_unray_demap);
        berr_ch(i) = b_err;
end
figure(7), drawnow, hold off

semilogy(ebN0dB, berr,'linewidth',2);
hold on;

semilogy(ebN0dB, berr_ch, '-o','linewidth',2);
grid on;
set(gca,'YScale','log'); % workaround for Matlab bug
title('BER over Rayleigh Fading Channels with AWGN noise');
xlabel('NR*SNR/dB');
ylabel('Bit Error Rate');

%% Uncomment until here

%% Implemaenting Rice Channel

K = [1 10 20 40]; % K parameter = Plos/PNlos
RiceKanal = radioFadingRiceChannel(nSamp, K);
ebN0abs = 10.^(ebN0dB./10);
%set ored, ogreen etc. 
error_colors = {'or-'; 'og-'; 'ob-'; 'om-'};
theor_colors = {'-r'; '-g'; '-b'; '-m'};

figure(8), hold off , clf

ideal = zeros(1, length(ebN0abs));  %integral von theoretical, formel unten
for k = 1:length(K)
    for i = 1:length(SNRProSymbol)
        disp(['k=',num2str(k),', i=',num2str(i)])
        y = generateBits(nBitsPerLoop);
        yc = mapper(y, const);
        n = radioFadingRiceChannel(length(yc), K(k));
        
        % Rice
        y_rice = (yc(:) .* n(:)).'; % nonconjugate transpose to get row vector
        y_rice_noise = setSNR(y_rice, SNRProSymbol(i));
        y_unrice = (y_rice_noise(:) ./ n(:)).';
        
        [index] = decision(y_unrice, const);
        y_unrice_demap = demapper(index, const);
        [~, ~, b_err] = countErrors(y, y_unrice_demap);
        berr_ch(k,i) = b_err;
        
        %Theortical
%         figure(9), drawnow, hold off
        
        theoretical = @(x, snr) ((((1+K(k)) .* sin(x).^2) ./ ((1+K(k)) .* sin(x).^2 + snr)) .* exp ( - (K(k) .* snr) ./ ((1+K(k)) .* sin(x).^2 + snr)));
        ideal(k,i) = (integral(@(x)theoretical(x, ebN0abs(i)), 0, (pi / 2))) / pi;
        
         figure(9)
%     
         semilogy(ebN0abs, ideal(k,:), (theor_colors{mod(k-1,length(theor_colors))+1}),'linewidth',2);
         axis([0 30 -10000 10])
         grid on;
         set(gca,'YScale','log'); % workaround for Matlab bug
         ylim([0.0001 0.1]);
         title('Theoretical BER over Rician Fading Channels with AWGN noise');
        xlabel('NR*SNR/dB');
        ylabel('Bit Error Rate');
        
    end

    figure(8), hold on
    semilogy(ebN0dB, berr_ch(k,:), (error_colors{mod(k-1,length(error_colors))+1}),'linewidth',2);drawnow
    grid on;
    set(gca,'YScale','log'); % workaround for Matlab bug
    title('BER over Rician Fading Channels with AWGN noise');
    xlabel('NR*SNR/dB');
    ylabel('Bit Error Rate');
    figure(9),hold on
           
 

   
end
 figure(8), drawnow, hold off
 figure(9), drawnow, hold off





%% Histogram con Kanalkoeffezient

figure(10)
subplot(1,2,1);
histogram(abs(RiceKanal), 100)
xlabel('amplitude of channel coefficients')
ylabel('No. of same samples')
title('Rayleighverteilte Amplitude der Kanalkoeffizienten')
subplot(1,2,2);
histogram(angle(RiceKanal), 100)
xlabel('phase of channel coefficients')
ylabel('No. of same samples')
title('Gleichverteilte Phase der Kanalkoeffizienten')



%% Hilfsmitteln 
%https://de.mathworks.com/matlabcentral/fileexchange/25293-matlab-for-digital-communication
%http://inhousernd.blogspot.com/2014/03/rotatedconstellation-demapper-specifics.html
%http://www.rfwireless-world.com/source-code/MATLAB/BPSK-QPSK-16QAM-64QAM-modulation-matlab-code.html
%Book: SIMULATION OF DIGITAL COMMUNICATION SYSTEMS USING MATLAB :: Mathuranathan Viswanathan
