function mask = getMask(domain,dLat,fLat,dLon,fLon)

mask = false(length(domain.lat),length(domain.lon));

[~,dLat] = min(abs(domain.lat-dLat));
[~,fLat] = min(abs(domain.lat-fLat));

[~,dLon] = min(abs(domain.lon-dLon));
[~,fLon] = min(abs(domain.lon-fLon));

mask(dLat:fLat,dLon:fLon) = 1;

endfunction