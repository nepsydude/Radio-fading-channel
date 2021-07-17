function y = decision(x, const)



y = zeros(size(x));

%Berechnung der mittleren Leistung
PconstAverage = mean(abs(const).^2);
PxAverage = mean(abs(x).^2);
%disp(x)
%Normierung des Signals
Xnormed = x.*sqrt(PconstAverage/PxAverage);

%Test
XnormedAverage = mean(abs(Xnormed).^2); 
% disp(['XnormedAverage=',num2str(XnormedAverage)])

%Berechnung des kleinsten Abstands vom Konstellationspunkt
for kx = 1:length(x)
    %x = x(kx)   
    for kc=1:length(const)
        dist(kc) = abs(Xnormed(kx)-const(kc)); %Berechnung der unterschiedlichen Strecke
    end
    
    [~, index] = min(dist); % Entscheidung 
    y(kx)= const(index);

end





