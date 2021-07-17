
%Declaring function to generate random Bitsreams up to n number of bits
%nBits

function y = generateBits(nBits)
%y = round (rand(1,nBits));
y = randi([0 1], 1, nBits); %random (0 or 1) in 1XN matrix
% figure(1)
% stem( y, 'linewidth',2)
% axis([0 70 0 1])
end

% 
% y = randi([0 1], 1, nBits); %random (0 or 1) in 1XN matrix
% figure(1);
% end