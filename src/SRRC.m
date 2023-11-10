function [modulatedBitStream, t] = SRRC(bitStream, T, sn, A2, alpha, K) %Square Root Raised Cosine
    Ts = T/sn; %Every bit has for example 1 sec duration and 32 samples.
    t = Ts:Ts:T*length(bitStream);    
    pulse = -K*T:Ts:K*T - Ts;

    % angle = pi * pulse * (1-alpha) / T;
    % 
    % 
    % 
    % x = zeros(1, length(pulse));
    % 
    % for i = 1 : length(pulse)
    %     if pulse(i) == 0
    %         x(i) = 1- alpha + 4*alpha/pi;
    %     elseif pulse(i) == T/(4*alpha) | pulse(i) == -T/(4*alpha)
    %         x(i) = (alpha/sqrt(2)) * ((1 + 2/pi) * sin(pi/(4*alpha)) + (1 - 2/pi) * cos(pi/(4*alpha)));
    %     else
    %         angle = pi * pulse(i) * (1-alpha) / T;
    %         x_num(i) = sin(angle) + (4*alpha*pulse(i)/T)*cos(angle);
    %         x_denom(i) = (pi*pulse(i)/T)*(1-(4*alpha *pulse(i)/T)^2);
    %         x(i) = x_num(i)/x_denom(i);
    %     end
    % end
    % 

    % x = x_num/x_denom;
    % for k = 1 : length(bitStream)
    %     if bitStream(k) == 1
    %         modulatedBitStream((k-1)*sn +1 : k*sn) =  A1*sin(pi*pulse/T);
    %     else
    %         modulatedBitStream((k-1)*sn +1 : k*sn) =  -A1*sin(pi*pulse/T);
    %     end
    % end


a = rcosdesign(alpha, 2*K, sn); % size of the output is 2 * K * sn

end