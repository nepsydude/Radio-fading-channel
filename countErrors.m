function [nErr, idx, ber] = countErrors(x, bits)
    idx = find(bits ~= x); % position / Indizes
    nErr = numel(idx); %  Anzahl der Element(bits)
    ber = nErr/numel(bits); % Bit Error Rate
end

