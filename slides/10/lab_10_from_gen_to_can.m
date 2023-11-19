syms x y x_n y_n real
% Ax^2 + Bxy + Cy^2 + Dx+ Ey + F = 0
A = 0;
B = 1;
C = 0;
D = 0;
E = 0;
F = 1;
% A = 34;
% B = 24;
% C = 41;
% D = -44;
% E = 58;
% F = 1;


eq_orig = A*x^2 + B*x*y + C*y^2 + D*x+ E*y + F==0;

Sb = [A B/2 D/2; B/2 C E/2; D/2 E/2 F]
Sm = [A B/2; B/2 C]

Sb_det = det(Sb)
Sm_eig = eig(Sm)

eq = x_n^2/(-Sb_det/(Sm_eig(1)^2*Sm_eig(2))) + y_n^2/(-Sb_det/(Sm_eig(2)^2*Sm_eig(1))) ==1 

ang = rad2deg(atan(B/(A-C)))/2

df_X = diff(eq_orig,x)
df_Y = diff(eq_orig,y)

[x_c,y_c] = solve([df_X df_Y],[x y])


shift_mat = [cosd(-ang) -sind(-ang) 0; sind(-ang) cosd(-ang) 0; 0 0 1] * ...
    [1 0 -x_c; 0 1 -y_c; 0 0 1]

shift_mat_one = [cosd(ang) sind(ang) -x_c*cosd(ang) + y_c*sind(ang);
                 -sind(ang) cosd(ang) x_c*sind(ang)-y_c*cosd(ang);
                 0,0,1]

x_n=0
y_n=-2

X_new_in_old = inv(shift_mat_one)*[x_n;y_n;1]
double(X_new_in_old)
