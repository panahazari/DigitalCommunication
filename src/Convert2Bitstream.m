function [numberOfGroups,BitStreams] = Convert2Bitstream(N,inputData, m, n, qbit)

numberOfGroups = 0;
BitStreams = [];
for i = 1:1
    % if N + i-1 <= m*n/64
    %     blockStream = reshape(inputData(:,:, i: i+ N-1), [], 1);    %Reshaping N blocks to a vector
    % 
    % else 
        blockStream = reshape(inputData(:,:, i: end), [], 1);
    % end
    blockBits = de2bi(blockStream);     %Convert the vector to bits
    
    BitStreams =reshape(blockBits.', 1, []);
    % bitStreamSize = N * 8 * 8 *qbit;
    % if length(blockBits) < bitStreamSize % To fix the last block's size problem
    %     blockBits = [blockBits; zeros(bitStreamSize - length(blockBits), 1)];
    % end
    % BitStreams(numberOfGroups + 1, :) = blockBits;  
    numberOfGroups  = numberOfGroups + 1;
end

