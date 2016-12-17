function ave = calave(x,lx,ly,Ny)
    [L,k] = size(x);
    nlen = Ny;
    reX = zeros(1,nlen);
    for j = 1:nlen
        sum = 0;
        %****************************
        if j>=1 && j<=L
            for n = 1:j
                sum = sum + x(n,j-n+1);
            end
            reX(j) = sum./j;
        end
        %****************************
        if j>L && j<k
            for n = 1:L
                sum = sum + x(n,j-n+1);
            end
            reX(j) = sum./L;
        end
        %****************************
        if j>=k && j<=nlen
            for n = j-k+1:L
                sum = sum + x(n,j-n+1);
            end
            reX(j) = sum./(nlen-j+1);
        end
    end
    ave = zeros(L,k);
    for i = 1:k
        ck = reX(i:i+L-1);
        ave(:,i) = ck;
    end 
 end