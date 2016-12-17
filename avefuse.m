function [data1 info1] = avefuse(data, info, n)
    mm = ceil(info.bands/n);
    info1 = info;
    info1.bands = mm;
    data_after = zeros(info.lines,info.samples,mm);
    for i = 1:mm
        sumdata = zeros(info.lines,info.samples);
        if i == mm
            num = info.bands - (mm-1)*n;
            for j = (mm-1)*n+1:info.bands
                sumdata = sumdata + data(:,:,j);
            end
            avedata = sumdata./num;
            data_after(:,:,i) = avedata;
        else
            for j = (i-1)*n+1:i*n
                sumdata = sumdata + data(:,:,j);
            end
            avedata = sumdata./n;
            data_after(:,:,i) = avedata;
        end        
    end
    data1 = data_after;
end
            