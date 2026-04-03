function C = row_kron(A,B)
%ROW_KRON perform and return the tensor product to each row of two 
% matrices, B1 and B2, assuming both matrices have the same number of rows
    n = size(A,1);
    p_A = size(A,2);
    p_B = size(B,2);
    C = zeros([n,p_A*p_B]);
    for i = 1:n
        C(i,:) = kron(A(i,:),B(i,:));
    end
end

