function [channel] = Channel3()
    h = zeros(1,8);
    h(1) = 1;
    h(2) = 0.4365;
    h(3) = 0.1905;
    h(4) = 0.0832;
    h(6) = 0.0158;
    h(8) = 0.003;
    channel = h;
end
