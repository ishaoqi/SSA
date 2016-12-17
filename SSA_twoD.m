function re2DX = SSA_twoD(X,Lx,Ly)
%    �Զ�ά����X���������׷���
%    X---ԭʼ��ά����
%    Lx��Ly----embedding���ڴ�С
    [Nx,Ny] = size(X);
    xx = Nx - Lx + 1;
    yy = Ny - Ly + 1;
    K = xx * yy;
    L = Lx * Ly;
    %******************��ȡX�Ĺ켣����***************************
    HH = zeros(Ly, yy, Nx);
    H = zeros(Ly, yy);
    for i = 1:Nx  %Nx���Ӿ���
        for j = 1:Ly  %�Ӿ��������
            for k = 1:yy %�Ӿ��������
                H(j,k) = X(i,k+j-1);
            end
        end
        HH(:,:,i) = H;
    end
    %**********���Ӿ������X�Ĺ켣����*************************
    for i = 1:Lx
        for j = 1:xx
            if j == 1
                A = HH(:,:,i+j-1);
            else
                Atemp = HH(:,:,i+j-1);
                A = [A,Atemp];   %�Ӿ���������
            end
        end
        if i==1
            B = A;
        else
            B = [B;A];           %�Ӿ��������ӣ�B��ΪX�Ĺ켣����
        end
    end
    %*****************����ֵ�ֽ�******************************
    [U,S,V] = svds(B,L);
    %����ÿ������ֵ��Ӧ��Xi
    recon = zeros(L,K,L);
    for i  = 1:L
        xi = U(:,i)*V(:,i)'.*S(i,i);
        recon(:,:,i) = xi;
    end
    % б�ߴ�ƽ������ȡÿ��Xi��Ӧ��ԭʼ�����ݣ����������㣬������Ҫ����ÿ����HbH�����б��ƽ��ֵ���ڼ��������ƽ��
    % ȡ��ÿ��H������ÿ��H��ƽ��
    Hnum = Lx*xx;
    Htotal = zeros(Ly,yy,Hnum,Ny);
    Have = Htotal;
    re2DX = zeros(Nx,Ny,L);
    for i = 1:L%��ÿ������ֵ��Ӧ��xi����
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
        %��ȡÿ��H��ƽ��ֵ֮�󣬽�H��������ΪXi
        for i1 = 1:Lx
            for i2 = 1:xx
                if i2 == 1
                    tA = Have(:,:,(i1-1)*xx+i2);
                else
                   tAtemp = Have(:,:,(i1-1)*xx+i2);
                   tA = [tA,tAtemp];   %�Ӿ��������� 
                end
            end
            if i1==1
                tB = tA;
            else
                tB = [tB;tA];           %�Ӿ��������ӣ�tB��ΪXi
            end
        end  
        %ÿ��H��ƽ��ֵ��ȡ֮������ȡ���ľ��󼶱��ƽ��ֵ
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