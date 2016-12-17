function reX = SSA_oneD(X,L)
%embedding with length of L
    nlen = length(X);
    if L >= nlen
        error('Hey! Too big a lag!');
    end
    k = nlen - L + 1;
    trX = zeros(L,k);
    for i = 1:k
        ck = X(i:i+L-1);
        trX(:,i) = ck;
    end  
    % ����ֵ�ֽ�
    % trX = USV'
    [U,S,V] = svds(trX,L);
   %δ���飬��ȡÿ��Xi
    recon = zeros(L,k,L);
    reX = zeros(L,nlen);
    for i  = 1:L
        xi = U(:,i)*V(:,i)'.*S(i,i);
        recon(:,:,i) = xi;
    end
    % б�ߴ�ƽ������ȡÿ��Xi��Ӧ��ԭʼһά����
    for i = 1:L
        for j = 1:nlen
            sum = 0;
            %****************************
            if j>=1 && j<=L
                for n = 1:j
                    sum = sum + recon(n,j-n+1,i);
                end
                reX(i,j) = sum./j;
            end
            %****************************
            if j>L && j<k
                for n = 1:L
                    sum = sum + recon(n,j-n+1,i);
                end
                reX(i,j) = sum./L;
            end
            %****************************
            if j>=k && j<=nlen
                for n = j-k+1:L
                    sum = sum + recon(n,j-n+1,i);
                end
                reX(i,j) = sum./(nlen-j+1);
            end
        end
    end
end