function output = ptDistance(input,pt)

switch(pt)
    case 'decoupage'
        if isstruct(input)
            [m,~] = size(input.masks);
            [obs,comp] = size(input.data);
            output = ones(obs,comp,m)*NaN;
            for i=1:obs
                for j=1:m
                    mask = input.masks(j,:);
                    mask(mask==0) = NaN;
                    output(i,:,j) = input.data(i,:).*mask;
                endfor
            endfor
        else
            %Log
            error(['(!) PRE-TRAITEMENT: ''',inputname(1),''' n''est pas une structure, il  faut définir ''',inputname(1),'.data'' ainsi ''',inputname(1),'.masks''.']);
        endif
        clear m obs comp i j;
    case 'distribution'
        if isstruct(input)
            output = hist(input.data(:),input.edges);
        else
            %Log
            error(['(!) PRE-TRAITEMENT: ''',inputname(1),''' n''est pas une structure, il  faut définir ''',inputname(1),'.data'' ainsi ''',inputname(1),'.edges''.']);
        endif
    case 'motif-local'
        m = mean(input(:));
        output = input-m;
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
                    endif
                endfor
            endfor
        endfor
    otherwise
        %Log
        error(['/!\ PRE-TRAITEMENT: ''',pt,''' est inconu..']);
endswitch

endfunction
