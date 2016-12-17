function re2DX = SSA_2D(X,Lx,Ly)
%    �Զ�ά����X���������׷���
%    input         X---ԭʼ��ά����,��СΪNx*Ny
%                  Lx��Ly----embedding���ڴ�С
%    output        re2DX-----Nx*Ny*L��С�ľ���LΪLx*Ly��Ϊ����ֵ������
%    ԭ��ο�����Jaime Zabalza 2015������
    [Nx,Ny] = size(X);
    if Lx>Nx/2 || Ly>Ny/2
        error('Lx����Ly����ԭʼ����ά�ȵ�һ�룬���޸ģ�');
    end
    xx = Nx - Lx + 1;
    yy = Ny - Ly + 1;
    K = xx * yy;
    L = Lx * Ly;
    %******************��ȡX�Ĺ켣����***************************
    HH = zeros(Ly, yy, Nx);
    H = zeros(Ly, yy);
    for i = 1:Nx  %Nx���Ӿ���
%         for j = 1:Ly  %�Ӿ��������
%             for k = 1:yy %�Ӿ��������
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
    %*****************����ֵ�ֽ�******************************
    %L = 10;
    [U,S,V] = svds(trX,L);
    %����ÿ������ֵ��Ӧ��Xi
    recon = zeros(L,K,L);
    for i  = 1:L
        xi = U(:,i)*V(:,i)'.*S(i,i);
        recon(:,:,i) = xi;
    end

    % б�ߴ�ƽ������ȡÿ��Xi��Ӧ��ԭʼ�����ݣ����������㣬������Ҫ����ÿ����HbH�����б��ƽ��ֵ���ڼ��������ƽ��
    % ȡ��ÿ��H������ÿ��H��ƽ��
    re2DX = zeros(Nx,Ny,L);
    for i = 1:L%��ÿ������ֵ��Ӧ��xi����
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
        %ÿ��H��ƽ��ֵ��ȡ֮������ȡ���ľ��󼶱��ƽ��ֵ
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