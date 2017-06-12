function [p3,h] = line2(k, sinPhi, cosPhi, s, d, plotCurve)
r = 1/k;
pM = [r * cosPhi; r * sinPhi; d];
p1 = [0;0;d];
a = s/r;
cosa = cos(a);
sina = sin(a);
% p2 per rotation matrix
% x = p1-pM;
% n = cross(pM, [0;0;1]);
% n = -n/norm(n);
% n1 = n(1); n2 = n(2); n3 = n(3);
%R = [n1*n1*(1-cosa)+cosa n1*n2*(1-cosa)-n3*sina n1*n3*(1-cosa)+n2*sina;
%     n2*n1*(1-cosa)+n3*sina n2*n2*(1-cosa)+cosa n2*n3*(1-cosa)-n1*sina;
%      n3*n1*(1-cosa)-n2*sina n3*n2*(1-cosa)+n1*sina n3*n3*(1-cosa)+cosa];

%p2 = R*x+pM;
%or shorter
% p2 = n*dot(n,x)+cosa*cross(cross(n,x),n)+sina*cross(n,x)+pM;

p2 = p1 + (pM-p1)*(1-cosa) + [0;0;r*sina];
v1 = p1 - pM;
v2 = p2 - pM;
p3 = -cross(v2, cross(v1,v2));
p3 = p3/norm(p3)*d + p2;
h = {};
if plotCurve
    segments = 100;
    v3 = cross(cross(v1,v2),v1);
    v3 = r*v3/norm(v3);
    t = linspace(0,atan2(norm(cross(v1,v2)),dot(v1,v2)),segments);
    v = v1*cos(t)+v3*sin(t)+pM;
    h = horzcat([0; 0; 0], v, p3);  
end
end