function blocks = tril2block(blockM,Nblocks)
% Converts a matrix in tril coordinate basis to a matrix in block (canonical) coordinates basis. (Used for fast operator inversion)

% Written by S.K. 07/2016
blocksize = [size(blockM,1)/Nblocks(1),size(blockM,2)/Nblocks(2)];
blocks = zeros(size(blockM)); 
for j = 0:Nblocks(1)-1;
    for k = 0:Nblocks(2)-1;
        idrow = j*blocksize(1)+1:(j+1)*blocksize(1);
        idcol = k*blocksize(2)+1:(k+1)*blocksize(2);
        blocks(idrow,idcol) = blockM(j+1:Nblocks(1):end,k+1:Nblocks(2):end);
    end
end
end

