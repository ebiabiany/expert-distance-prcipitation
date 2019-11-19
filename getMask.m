function mask = getMask(domain,dLat,fLat,dLon,fLon)

mask = false(length(domain.lat),length(domain.lon));

[~,dLat] = min(abs(domain.lat-dLat));
[~,fLat] = min(abs(domain.lat-fLat));

[~,dLon] = min(abs(domain.lon-dLon));
[~,fLon] = min(abs(domain.lon-fLon));

% disp(['dLat=',num2str(dLat),' fLat=',num2str(fLat),' dLon=',num2str(dLon),' fLon=',num2str(fLon)]);

mask(dLat:fLat,dLon:fLon) = 1;

end