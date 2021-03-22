function reconstructed_message = QPSK_2_ascii(mprime)
    N = length(mprime);
    reconstructed_message = [];
    for iter = 1:N
       if(mprime(iter) == -1 -1i)
           reconstructed_message = [reconstructed_message,'00'];
       elseif(mprime(iter) == -1 + 1j)
           reconstructed_message = [reconstructed_message,'10'];
       elseif(mprime(iter) == 1 - 1j)
           reconstructed_message = [reconstructed_message,'11'];
       else
           reconstructed_message = [reconstructed_message,'01'];
       end
    end
    reconstructed_message = [reconstructed_message,'01'];
  
end