function a = alpha(A)
    n = length(A);
    B = zeros(nnz(A)/2,n);
    b = ones(length(B),1);
    count = 1;
    options = optimoptions('intlinprog');options.Display = 'off';
    for i = 2:n
        for j = 1:i-1
            if(A(i,j)==1)
                B(count,[i j]) = [1 1];
                count = count + 1;
            end
        end
    end
    [~,a] = intlinprog(-ones(n,1),1:n,B,b,[],[],zeros(n,1),ones(n,1),options);
    a = -a;
end