function y = radioFadingChannel(nSamp)
%RADIOFADINGCHANNEL Summary 
    y = (randn(1, nSamp) + randn(1, nSamp) .* 1j); %Generieren des nSamp=Kanalkoeffizient
    p0 = mean(abs(y) .^2);  %Power =1/N ?|y|2  P0 ->2
    alpha = sqrt(1./p0);  %Skalierungsfaktor alpha, Prozess für p1=1
    %p1 = mean (abs(alpha.*y).^2);  %Power entspricht 1
    
    y = y.*alpha;
    
    
%     power = mean(abs(y) .^2);
%     y = y ./ sqrt(power);
%     mean(abs(y) .^2);
%     disp(power)
end