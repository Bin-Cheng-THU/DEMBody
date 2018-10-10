function Dy = contact_model(t,y)

global E nu r m

%% Model
Kn = -2*E*sqrt(r)/(3*(1-nu^2));

Dy = zeros(2,1);
Dy(1) = y(2);
Dy(2) = Kn/m*y(1)^(3/2);
end