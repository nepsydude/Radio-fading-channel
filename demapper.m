function y = demapper(x, const)   
    % initialize bitstream of y
    k = log2(length(const));
    % y = zeros(1, length(x)*k);
     y = [0];
    
    % get index of constellation vector for each symbol in signal vector
    for i=1:length(x)
        index = find(const == x(i));
        index = uint16(index);      % integer value necessary for bitget()
        yBit = bitget(index - 1, k:-1:1, 'uint16');    % LSB left, get first k bits
        %yBit = yBit(k:-1:1);              % LSB
        y = cat(2, y, yBit); %concatenate yBit to the end of y along dimension 2 
    end
    y = y(2:end);   
end