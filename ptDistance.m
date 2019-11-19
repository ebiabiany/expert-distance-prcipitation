function output = ptDistance(input,pt)
% [FUNCTION] Pré-traitemnt pour calcule de la distance particulière
% output = ptDistance(input,pt)
%  input : donnée d'entrée
%        matrix ou struct
%     pt : option en texte
%        'motif-local', 'distribution', 'fluctuation-jour', etc.
% output : données en sortie
%        matrix ou struct

switch(pt)
    case 'decoupage'
        %Verification
        if isstruct(input)
            %Traitement
            [m,~] = size(input.masks);
            [obs,comp] = size(input.data);
            output = ones(obs,comp,m)*NaN;
            for i=1:obs
                for j=1:m
                    mask = input.masks(j,:);
                    mask(mask==0) = NaN;
                    output(i,:,j) = input.data(i,:).*mask;
                end
            end
        else
            %Log
            error(['(!) PRE-TRAITEMENT: ''',inputname(1),''' n''est pas une structure, il  faut définir ''',inputname(1),'.data'' ainsi ''',inputname(1),'.masks''.']);
        end
        %Nettoyage
        clear m obs comp i j;
    case 'distribution'
        %Verification
        if isstruct(input)
            %Traitement
            output = hist(input.data(:),input.edges);
        else
            %Log
            error(['(!) PRE-TRAITEMENT: ''',inputname(1),''' n''est pas une structure, il  faut définir ''',inputname(1),'.data'' ainsi ''',inputname(1),'.edges''.']);
        end
    case 'motif-local'
        %Traitement
        m = mean(input(:));
        output = input-m;
        %Nettoyage
        clear m;
    case 'fluctuation-jour'
        %Log
        error(['(!) PRE-TRAITEMENT: ''',pt,''' est indisponible pour le moment..']);
    case 'smooth'
        data = input.data;
        edges = input.edges;
        [l,c] = size(data);
        for i=1:l
            for j=1:c
                for k=1:length(edges)-1
                    if(data(i,j)>edges(k) && data(i,j)<edges(k+1))
                        data(i,j) = (edges(k+1)+edges(k))/2;
                    end
                end
            end
        end
    otherwise
        %Log
        error(['/!\ PRE-TRAITEMENT: ''',pt,''' est inconu..']);
end

%Log
%disp(['ptDistance(',inputname(1),',''',pt,''') -> done.']);

end