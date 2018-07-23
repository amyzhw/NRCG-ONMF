function ht = transpVect(prob, x, h, d, type)


%% z is the foot of ht
if nargin == 4
    type = 0;
end
if type == 0    
    z = moveEIG(prob,x,d,1);
else % type =1
    z = d;
end
%% the given vector is (x,h) and is transported onto d

ip_old = ip(x,h,h);


M1 = (z.U'*x.U)*h.M*(x.V'*z.V);


Up2 = h.Up*(x.V'*z.V);

Vp3 = h.Vp*(x.U'*z.U);


ht.M = M1;
ht.Up = Up2;
ht.Up = ht.Up - z.U*(z.U'*ht.Up);
ht.Vp = Vp3;
ht.Vp = ht.Vp - z.V*(z.V'*ht.Vp);

ip_new = ip(z,ht,ht);
ht = scaleTxM(ht, ip_old/ip_new);
