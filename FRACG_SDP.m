% cd('/home/pro=/Desktop/Dropbox/Parthe RA/Numerical Simulations/SDPT3-4.0');
% startup;
%% N_+(FRAC(G))


A = [0 1 0 0 1; 1 0 1 0 0; 0 1 0 1 0; 0 0 1 0 1; 1 0 0 1 0]; n =5;
% n = 12;     A = graph(n);
disp([num2str(nnz(A)/2) ' edges in the graph'])

%% Relaxation

e = nnz(A)/2;
n2 = 0;
blk{1,1} = 's';      blk{1,2} = 1 + n;

As = cell(0);
b = [];
constraint_ind = 0;

% X_11 = 1
constraint_ind = constraint_ind + 1;
As{constraint_ind} = sparse(1,1,1,n+1,n+1);
b(constraint_ind) = 1;

% X_1i = X_ii
for i = 2:n+1
    constraint_ind = constraint_ind + 1;
    As{constraint_ind} = sparse([1  i],[i  i],[1/2 -1],n+1,n+1);
    b(constraint_ind) = 0;
    
end
eq_constraints = constraint_ind;

% inequality constraints
% linearized inequalities from 0 <= x_i <= 1
for i = 2:n+1
        
        % -X_ii <= 0
        constraint_ind = constraint_ind + 1;
        As{constraint_ind} = sparse(i,i,-1,n+1,n+1);
        b(constraint_ind) = 0;
        n2 = n2 + 1;
        
        % X_ii <= 1
        constraint_ind = constraint_ind + 1;
        As{constraint_ind} = sparse(i,i,1,n+1,n+1);
        b(constraint_ind) = 1;
        n2 = n2 + 1;
        
    for j = 2:i-1
        
        % -X_ji <= 0
        constraint_ind = constraint_ind + 1;
        As{constraint_ind} = sparse(j,i,-1,n+1,n+1);
        b(constraint_ind) = 0;
        n2 = n2 + 1;
        
        % X_ji <= X_ii
        constraint_ind = constraint_ind + 1;
        As{constraint_ind} = sparse([j j],[i j],[1/2 -1],n+1,n+1);
        b(constraint_ind) = 0;
        n2 = n2 + 1;
        
        % X_ji <= X_jj
        constraint_ind = constraint_ind + 1;
        As{constraint_ind} = sparse([j i],[i i],[1/2 -1],n+1,n+1);
        b(constraint_ind) = 0;
        n2 = n2 + 1;
        
        % X_1i + X_1j - X_ji <=  1
        constraint_ind = constraint_ind + 1;
        As{constraint_ind} = sparse([1 1 j],[j i i],[1/2 1/2 -1/2],n+1,n+1);
        b(constraint_ind) = 1;
        n2 = n2 + 1;
        
    end
end


% linearized inequalities from x_i + x_j <= 1 for (i,j) in E
for k = 2:n+1
    for i = 3:n+1
        for j = 2:i - 1
        
            if (A(i-1,j-1) == 1)

            % X_ik + X_jk <= X_1k
                constraint_ind = constraint_ind + 1;
                As{constraint_ind} = sparse([min(j,k) min(i,k) 1],[max(j,k) max(i,k) k],...
                    [ifelse(j==k,1,1/2) ifelse(i==k,1,1/2) -1/2],n+1,n+1);
                b(constraint_ind) = 0;
                n2 = n2 + 1;
                
            % X_1i + X_1j - X_ik - X_jk + X_1k <= 1
                constraint_ind = constraint_ind + 1;
                As{constraint_ind} = sparse([1 1 min(j,k) min(i,k) 1],[i j max(j,k) max(i,k) k],...
                    [1/2 1/2 -ifelse(j==k,1,1/2) -ifelse(i==k,1,1/2) 1/2],n+1,n+1); 
                b(constraint_ind) = 1;
                n2 = n2 + 1;
                
            end
        end
    end
end

blk{2,1} = 'l';      blk{2,2} = n2;

b = b';
At(1) = svec(blk(1,:),As,1);
At{2,1} = [sparse(n2,eq_constraints), speye(n2)];
C{1,1} = -sparse(2:n+1,2:n+1,ones(1,n),n+1,n+1);
C{2,1} = zeros(n2,1);
% tic
disp(['Alpha(G) is ' num2str(alpha(A))])
% toc
disp(['Theta(G) is ' num2str(theta(A))]);
% toc
[~,obj,X,y,Z] = evalc('sqlp(blk,At,C,b)');
disp(['TH_FRAC(G) is ' num2str(-obj(1))]);
% toc
clearvars -except X A