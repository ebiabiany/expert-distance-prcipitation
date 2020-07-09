function output = getDistance(input1,input2,type)

switch (type)
    case 'l1'
        output = norm(input2-input1,1);
    case 'l2'
        output = norm(input2-input1,2);
    case 'mahalanobis'
        S = cov(input2);
        mu = mean(input2,1);
        X = sqrt(((input1-mu)/S)*(input1-mu)');
        S = cov(input1);
        mu = mean(input1,1);
        Y = sqrt(((input2-mu)/S)*(input2-mu)');
        output = (X+Y)./2;
    case 'spearman'
        output = pdist([input1;input2],'spearman');
    case 'kullback'
        if (size(input1,2)~=size(input2,2))
            error('the number of columns in input1 and input2 should be the same');
        endif
        if (sum(~isfinite(input1(:))) + sum(~isfinite(input2(:))))
           error('the inputs contain non-finite values!') 
        endif
        if size(input2,1)==1
            input2 = input2 ./sum(input2);
            input1 = input1 ./repmat(sum(input1,2),[1 size(input1,2)]);
            temp =  input1.*log2(input1./repmat(input2,[size(input1,1) 1]));
            temp(isnan(temp))=0;
            X = sum(temp,2);
        elseif (size(input2,1)==size(input1,1))
            input2 = input2 ./repmat(sum(input2,2),[1 size(input2,2)]);
            input1 = input1 ./repmat(sum(input1,2),[1 size(input1,2)]);
            temp =  input1.*log2(input1./input2);
            temp(isnan(temp))=0;
            X = sum(temp,2);
        endif
        C = input1;
        input1 = input2;
        input2 = C;
        if (size(input1,2)~=size(input2,2))
            error('the number of columns in input1 and input2 should be the same');
        endif
        if (sum(~isfinite(input1(:))) + sum(~isfinite(input2(:))))
           error('the inputs contain non-finite values!') 
        endif
        if (size(input2,1)==1)
            input2 = input2 ./sum(input2);
            input1 = input1 ./repmat(sum(input1,2),[1 size(input1,2)]);
            temp =  input1.*log2(input1./repmat(input2,[size(input1,1) 1]));
            temp(isnan(temp))=0;
            Y = sum(temp,2);
        elseif (size(input2,1)==size(input1,1))
            input2 = input2 ./repmat(sum(input2,2),[1 size(input2,2)]);
            input1 = input1 ./repmat(sum(input1,2),[1 size(input1,2)]);
            temp =  input1.*log2(input1./input2);
            temp(isnan(temp))=0;
            Y = sum(temp,2);
        endif
        output = X+Y;
    case 'kl'
        ep = 1/(2176*2);
        n = length(input1);
        for i=1:n
            input1(i) = (input1(i) + ep)/(1+n*ep);
            input2(i) = (input2(i) + ep)/(1+n*ep);
        endfor
        X = 0;
        Y = 0;
        for i=1:length(input1)
                X = X + input1(i)*log2(input1(i)/input2(i));
                Y = Y + input2(i)*log2(input2(i)/input1(i));
        endfor
        output = X + Y;
    case 'd-kl'
        input = input1;
        J1 = input2(1,1);
        J2 = input2(1,2);
        n = length(input.cube);
        load('edges.mat');
        s = 0;
        for i=1:n
            j1 = histc(input.cube(i).data(J1,:),edges,'Normalization','probability');
            j2 = histc(input.cube(i).data(J2,:),edges,'Normalization','probability');
            s = s + getDistance(j1,j2,'kl');
        endfor
        output = s / n;
    otherwise
        %Log
        error(['/!\ GET-DISTANCE: ''',type,''' is unknown..']);
endswitch

endfunction
