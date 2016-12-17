function re2DX = SSA_2D(X,Lx,Ly)
%    对二维数据X进行奇异谱分析
%    input         X---原始二维数据,大小为Nx*Ny
%                  Lx，Ly----embedding窗口大小
%    output        re2DX-----Nx*Ny*L大小的矩阵，L为Lx*Ly，为奇异值的数量
%    原理参考论文Jaime Zabalza 2015年论文
    [Nx,Ny] = size(X);
    if Lx>Nx/2 || Ly>Ny/2
        error('Lx或者Ly大于原始数组维度的一半，请修改！');
    end
    xx = Nx - Lx + 1;
    yy = Ny - Ly + 1;
    K = xx * yy;
    L = Lx * Ly;
    %******************获取X的轨迹矩阵***************************
    HH = zeros(Ly, yy, Nx);
    H = zeros(Ly, yy);
    for i = 1:Nx  %Nx个子矩阵
%         for j = 1:Ly  %子矩阵的行数
%             for k = 1:yy %子矩阵的列数
%                 H(j,k) = X(i,k+j-1);
%             end
%         end
        xrow = X(i,yy:-1:1);
        ycol = X(i,yy:Ny);
        H = toeplitz(ycol,xrow);        
        HH(:,:,i) = fliplr(H);
    end
    trX_cell = cell(Lx,xx);
    for i = 1:Lx
        for j = 1:xx
            trX_cell{i,j} = HH(:,:,i+j-1);
        end
    end      
    trX = cell2mat(trX_cell);
    %*****************奇异值分解******************************
    %L = 10;
    [U,S,V] = svds(trX,L);
    %计算每个奇异值对应的Xi
    recon = zeros(L,K,L);
    for i  = 1:L
        xi = U(:,i)*V(:,i)'.*S(i,i);
        recon(:,:,i) = xi;
    end

    % 斜线处平均，获取每个Xi对应的原始的数据，分两步计算，首先需要计算每个子HbH矩阵的斜线平均值，在计算矩阵间的平均
    % 取出每个H，并对每个H求平均
    re2DX = zeros(Nx,Ny,L);
    for i = 1:L%对每个特征值对应的xi处理
        xi = recon(:,:,i);
        nrow = ones(1,Lx) .* Ly;
        ncol = ones(1,xx) .* yy;
        xi_cell = mat2cell(xi,nrow,ncol);
        for ii = 1:Lx
            for jj = 1:xx
                H = xi_cell{ii,jj};
                H_ave = calave(H,Lx,Ly,Ny);
                xi_cell{ii,jj} = H_ave;
            end
        end 
        %每个H的平均值求取之后，再求取外层的矩阵级别的平均值
        nlen = Nx;
        reX = zeros(Ly,yy,nlen);
        k = xx;
        for j1 = 1:nlen
            sum = 0;
            %****************************
            if j1>=1 && j1<=Lx
                for n = 1:j1
                    sum = sum + xi_cell{n,j1-n+1};                    
                end
                reX(:,:,j1)  = sum./j1;                
            end
            %****************************
            if j1>Lx && j1<k
                for n = 1:Lx
                    sum = sum + xi_cell{n,j1-n+1};
                end
                reX(:,:,j1) = sum./Lx;
            end
            %****************************
            if j1>=k && j1<=nlen
                for n = j1-k+1:Lx
                    sum = sum +  xi_cell{n,j1-n+1};
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