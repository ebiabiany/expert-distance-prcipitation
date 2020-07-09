function masks = getMasks(domain,Vedges,Hedges)

masks = zeros((length(Vedges)-1)*(length(Hedges)-1),(length(domain.lat))*(length(domain.lon)));
k = 1;

for i=1:length(Vedges)-1
    for j=1:length(Hedges)-1
        data = getMask(domain,Vedges(i),Vedges(i+1),Hedges(j),Hedges(j+1));
        masks(k,:) = data(:);
        imwrite(data,['mask_EXP1',num2str(k),'.bmp']);
        k = k + 1;
    endfor
endfor

endfunction