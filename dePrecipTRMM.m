function dist = dePrecipTRMM(j1,j2)

m = 4;

if ismatrix(j1) && ismatrix(j2)

    if sum(isnan(j1))>1
        dist = NaN;
        return;
    elseif sum(isnan(j2))>1
        dist = NaN;
        return;      
    end

    %# CONFIG
    load('domainBx0,25.mat');
    input.masks = getMasks(domainB,[30.125,20.125,5.125],[-66.625,-50.152,-20.125]);

    %#TRAITEMENT
    input.data = [j1;j2];
    output = ptDistance(input,'decoupage');

    clear input;

    [obs,~] = size(output);
    input.cube = struct([]);

    for j=1:m
        eff = output(1,:,j);
        eff(isnan(eff)) = [];
        input.cube(j).data = zeros(obs,length(eff));
        input.cube(j).data(1,:) = eff;
        for i=2:obs
            eff = output(i,:,j);
            eff(isnan(eff)) = [];
            input.cube(j).data(i,:) = eff;
        end
    end
    clear output eff obs comp;

    load('edges.mat');
    
    s = 0;
    for i=1:m
        j1 = histcounts(input.cube(i).data(1,:),edges,'Normalization','probability');
        j2 = histcounts(input.cube(i).data(2,:),edges,'Normalization','probability');
        s = s + getDistance(j1,j2,'kl');
    end
    dist = s / m;
    
elseif ndims(j1)== 3 && ndims(j2)== 3
    s = 0;
    for i=1:m
        s = s + getDistance(j1(:,:,i),j2(:,:,i),'kl');
    end
    
    dist = s / m;
else
    error('(!)dePrecipTRMM : Problème de dimension.. ');
end
end
