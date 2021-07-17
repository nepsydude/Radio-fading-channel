function y = radioFadingRiceChannel(nSamp, K)

coeffNLOS = randn(nSamp,1) + 1i*randn(nSamp,1); %Kanalkoeffizienten(NLOS)

% Normierung von LOS-Amplitude auf P_LOS = K * P_NLOS 
    PNLOS = mean(abs(coeffNLOS) .^ 2); % gegebene (lineare) mittlere Leistung NLOS-Komponente
	PLOS = K * PNLOS;
    coeffLOS = sqrt(PLOS) .*exp(1j.*2*pi.*rand(nSamp,1)); % Phase dreht sich, wegen Dopplerfrequenz; Amplitude konstant mit K-fachem der mittleren(!) NLOS-Leistung 

    % NLOS und LOS Komponente addieren
    coeff = coeffNLOS + coeffLOS;
    %checkPnAveragecoeff = mean(abs(coeff) .^ 2);
    
    % Normierung auf mittlere Leistung=1
    PnAverageCoeff = mean(abs(coeff) .^ 2); % mittlere Leistung berechnen
    normed_coeff = coeff ./ sqrt(PnAverageCoeff);  % Normierte, komplexe Kanalkoeffizienten
    %PnAverageNormedCoeff = mean(abs(normed_coeff) .^ 2);  % Ueberpruefen, ob Normierung funktioniert hat
    y = normed_coeff;
    
end