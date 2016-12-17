function re2DX = SSA_twoD(X,Lx,Ly)
%    对二维数据X进行奇异谱分析
%    X---原始二维数据
%    Lx，Ly----embedding窗口大小
    [Nx,Ny] = size(X);
    xx = Nx - Lx + 1;
    yy = Ny - Ly + 1;
    K = xx * yy;
    L = Lx * Ly;
    %******************获取X的轨迹矩阵***************************
    HH = zeros(Ly, yy, Nx);
    H = zeros(Ly, yy);
    for i = 1:Nx  %Nx个子矩阵
        for j = 1:Ly  %子矩阵的行数
            for k = 1:yy %子矩阵的列数
                H(j,k) = X(i,k+j-1);
            end
        end
        HH(:,:,i) = H;
    end
    %**********由子矩阵组成X的轨迹矩阵*************************
    for i = 1:Lx
        for j = 1:xx
            if j == 1
                A = HH(:,:,i+j-1);
            else
                Atemp = HH(:,:,i+j-1);
                A = [A,Atemp];   %子矩阵行连接
            end
        end
        if i==1
            B = A;
        else
            B = [B;A];           %子矩阵列连接，B即为X的轨迹矩阵
        end
    end
    %*****************奇异值分解******************************
    [U,S,V] = svds(B,L);
    %计算每个奇异值对应的Xi
    recon = zeros(L,K,L);
    for i  = 1:L
        xi = U(:,i)*V(:,i)'.*S(i,i);
        recon(:,:,i) = xi;
    end
    % 斜线处平均，获取每个Xi对应的原始的数据，分两步计算，首先需要计算每个子HbH矩阵的斜线平均值，在计算矩阵间的平均
    % 取出每个H，并对每个H求平均
    Hnum = Lx*xx;
    Htotal = zeros(Ly,yy,Hnum,Ny);
    Have = Htotal;
    re2DX = zeros(Nx,Ny,L);
    for i = 1:L%对每个特征值对应的xi处理
        for j = 1:Lx
            for k = 1:xx
                xi = recon(:,:,i);
                Htotal(:,:,(j-1)*xx+k) = xi((j-1)*Ly+1:j*Ly,(k-1)*yy+1:k*yy);
            end
        end
        for ii = 1:Hnum
            htemp = Htotal(:,:,ii);
            Have(:,:,ii) = calave(htemp,Lx,Ly);
        end
        %求取每个H的平均值之后，将H重新连接为Xi
        for i1 = 1:Lx
            for i2 = 1:xx
                if i2 == 1
                    tA = Have(:,:,(i1-1)*xx+i2);
                else
                   tAtemp = Have(:,:,(i1-1)*xx+i2);
                   tA = [tA,tAtemp];   %子矩阵行连接 
                end
            end
            if i1==1
                tB = tA;
            else
                tB = [tB;tA];           %子矩阵列连接，tB即为Xi
            end
        end  
        %每个H的平均值求取之后，再求取外层的矩阵级别的平均值
        nlen = Nx;
        k = xx;
        reX = zeros(Ly,yy,nlen);
        for j1 = 1:nlen
            sum = 0;
            %****************************
            if j1>=1 && j1<=Lx
                for n = 1:j1
                    sum = sum + tB((n-1)*Ly+1:n*Ly,(j1-n+1-1)*yy+1:(j1-n+1)*yy);                    
                end
                reX(:,:,j1)  = sum./j1;                
            end
            %****************************
            if j1>L && j1<k
                for n = 1:Lx
                    sum = sum + tB((n-1)*Ly+1:n*Ly,(j1-n+1-1)*yy+1:(j1-n+1)*yy);
                end
                reX(:,:,j1) = sum./Lx;
            end
            %****************************
            if j1>=k && j1<=nlen
                for n = j1-k+1:Lx
                    sum = sum + tB((n-1)*Ly+1:n*Ly,(j1-n+1-1)*yy+1:(j1-n+1)*yy);
                end
               reX(:,:,j1)  = sum./(nlen-j1+1);
            end
        end
       
        for kk = 1:nlen
            temp = reX(:,:,kk);
            vector = zeros(1,Ny);
            for k1 = 1:Ny
                if k1<=yy
                   vector(1,k1) = temp(1,k1);
                else
                    for mm = 2:Ly
                        vector(1,k1) = temp(mm,yy);
                    end
                end
            end
            if kk == 1
                map = vector;
            else
                map = [map;vector];
            end
        end
        re2DX(:,:,i) = map; 
    end  
end