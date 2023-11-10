function [impulse_train] = bit_sequence_to_impulses(bit_sequence, sn)
    impulse_train = [];
    n = [];
    for i = 1:length(bit_sequence)
        if bit_sequence(i) == 1
            impulse_train = [impulse_train, +1, zeros(1, sn - 1)];  % positive impulse for 1
        elseif bit_sequence(i) == 0
            impulse_train = [impulse_train, -1, zeros(1, sn - 1)];  % negative impulse for 0
        else
            error('Input bit sequence should only contain 0s and 1s.');
        end
        n = [n, (i-1)*sn + (1:sn)];
    end
end