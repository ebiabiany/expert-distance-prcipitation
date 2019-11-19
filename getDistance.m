function output = getDistance(input1,input2,type)
% [FUNCTION] Pré-traitemnt pour calcule de la distance particulière
% output = ptDistance(input,pt)
%  input1 : donnée d'entrée
%        matrix
%  input2 : donnée d'entrée
%        matrix
%   type : option en texte
%        'motif-local', 'distribution', 'fluctuation-jour', etc.
% output : données en sortie
%        matrix

switch type
    case 'l1'
        %Traitement
        output = norm(input2-input1,1);
    case 'l2'
        %Traitement
        output = norm(input2-input1,2);
    case 'mahalanobis'
        %Traitement
        S = cov(input2);
        mu = mean(input2,1);
        X = sqrt(((input1-mu)/S)*(input1-mu)');
        S = cov(input1);
        mu = mean(input1,1);
        Y = sqrt(((input2-mu)/S)*(input2-mu)');
        output = (X+Y)./2;
        %output = mahal(input1,input2);
    case 'spearman'
        output = pdist([input1;input2],'spearman');
    case 'kullback'
        if size(input1,2)~=size(input2,2)
            error('the number of columns in input1 and input2 should be the same');
        end
        if sum(~isfinite(input1(:))) + sum(~isfinite(input2(:)))
           error('the inputs contain non-finite values!') 
        end
        % normalizing the input1 and input2
        if size(input2,1)==1
            input2 = input2 ./sum(input2);
            input1 = input1 ./repmat(sum(input1,2),[1 size(input1,2)]);
            temp =  input1.*log2(input1./repmat(input2,[size(input1,1) 1]));
            temp(isnan(temp))=0;% resolving the case when input1(i)==0
            X = sum(temp,2);
        elseif size(input2,1)==size(input1,1)
            input2 = input2 ./repmat(sum(input2,2),[1 size(input2,2)]);
            input1 = input1 ./repmat(sum(input1,2),[1 size(input1,2)]);
            temp =  input1.*log2(input1./input2);
            temp(isnan(temp))=0; % resolving the case when input1(i)==0
            X = sum(temp,2);
        end
        C = input1;
        input1 = input2;
        input2 = C;
        if size(input1,2)~=size(input2,2)
            error('the number of columns in input1 and input2 should be the same');
        end
        if sum(~isfinite(input1(:))) + sum(~isfinite(input2(:)))
           error('the inputs contain non-finite values!') 
        end
        % normalizing the input1 and input2
        if size(input2,1)==1
            input2 = input2 ./sum(input2);
            input1 = input1 ./repmat(sum(input1,2),[1 size(input1,2)]);
            temp =  input1.*log2(input1./repmat(input2,[size(input1,1) 1]));
            temp(isnan(temp))=0;% resolving the case when input1(i)==0
            Y = sum(temp,2);
        elseif size(input2,1)==size(input1,1)
            input2 = input2 ./repmat(sum(input2,2),[1 size(input2,2)]);
            input1 = input1 ./repmat(sum(input1,2),[1 size(input1,2)]);
            temp =  input1.*log2(input1./input2);
            temp(isnan(temp))=0; % resolving the case when input1(i)==0
            Y = sum(temp,2);
        end
        %disp(['>> X=',num2str(X),' Y=',num2str(Y)]);
        output = X+Y;
    case 'kl'
        ep = 1/(2176*2);
        n = length(input1);
        for i=1:n
            input1(i) = (input1(i) + ep)/(1+n*ep);
            input2(i) = (input2(i) + ep)/(1+n*ep);
        end
        X = 0;
        Y = 0;
        for i=1:length(input1)
                X = X + input1(i)*log2(input1(i)/input2(i));
                Y = Y + input2(i)*log2(input2(i)/input1(i));
        end
        output = X + Y;
    case 'd-kl'
        %INFORMATION : input1 = input & input2 = [J1,J2]
        input = input1;
        J1 = input2(1,1);
        J2 = input2(1,2);
        n = length(input.cube);
        edges = load('save_work_09_01_2019_1.mat', 'edges');
        edges = edges.edges;
        s = 0;
        for i=1:n
            j1 = histcounts(input.cube(i).data(J1,:),edges,'Normalization','probability');
            j2 = histcounts(input.cube(i).data(J2,:),edges,'Normalization','probability');
            s = s + getDistance(j1,j2,'kl');
        end
        output = s / n;
    otherwise
        %Log
        error(['/!\ GET-DISTANCE: ''',type,''' est inconu..']);
end

%Log
% disp(['getDistance(',inputname(1),',',inputname(2),',''',type,''') -> ',num2str(output),'.']);

end