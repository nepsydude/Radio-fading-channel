function y = setSNR(x, snrdB)
%
% SNRdB = 10*log10(Ps/Pn)   

    Ps = mean(abs(x) .^ 2); %Signalleistung
    Pn = Ps * 10.^(-snrdB / 10); %Rauschleistung
    noise = randn(1, length(x)) + randn(1, length(x)) .* 1j;
    PnMeas = mean(abs(noise) .^ 2);
    
    noise = noise .* sqrt(Pn / PnMeas); %Pn/PnMeas: Spektralsdichte
    
%     noiseAverage = mean(abs(noise).^2); %Geräusch Mittelwert
%     disp(noiseAverage);

    y = x + noise;
end

