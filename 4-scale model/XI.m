% cos of included angle between SZA & theta, 
% phi the zimuthal differences
function diff=XI(SZA,theta,phi)
cosdif=cos(SZA)*cos(theta)+sin(SZA)*sin(theta)*cos(phi);
diff=acos(cosdif);


