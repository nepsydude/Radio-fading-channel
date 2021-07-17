function y = mapper(bits, const)
%The mapper functions as an encoder and arranges nBitsPerSymb = log2(nConst)
%binary symbols into one nConst-ary symbol

nConst = length(const);
nBits = length(bits);

nBitsPerSymb = log2(nConst); % number of bits per symbol

%Umwandlung in Matrix-Form eines Vektors 
bitMatrix = reshape(bits,nBitsPerSymb, (nBits/nBitsPerSymb)).'; %Transponiert
dec = bi2de(bitMatrix, 'left-msb'); % convert the bit stream into symbol stream
[bitMatrix dec];
index = dec + 1 ; % addieren von 1 um index 1-4 zu bekommen
y = const(index);
 
 end  