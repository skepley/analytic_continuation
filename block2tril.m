function blockM = block2tril(blocks)
% Converts blocks of triangular operators (block basis) to a lower triangular matrix (triL basis) for fast inversion. Invert this by blockM\eye(MN)

% Written by S.K. 07/2016
Nblocks = size(blocks);
blocksize = size(blocks{1,1});
blockM = zeros(Nblocks(1)*blocksize(1),Nblocks(2)*blocksize(2));

for j = 1:Nblocks(1)
    for k = 1:Nblocks(2)
        blockM(j:Nblocks(1):end,k:Nblocks(2):end) = blocks{j,k};
    end
end
end

