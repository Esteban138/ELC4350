function m = QPSK_converter(str)
N=length(str);                           % length of string
%m=zeros(1,4*N);
m = [];
%str = '01234 I wish I were an Oscar Meyer wiener 56789';
for k=1:N 
    binary_string = dec2bin(str(k));
    value = sprintf('%08s', binary_string);
    for it = 1:2:8
        if(value(it:it+1) == '00')
            m(end+1) = -1 -1i; 
        elseif(value(it:it+1) == '01')
            m(end+1) = 1 + 1i; 
        elseif(value(it:it+1) == '10')
            m(end+1) = -1 + 1i;
        else(value(it:it+1) == '11')
            m(end+1) = 1 - 1i;
            
        end
        
    end

end
end