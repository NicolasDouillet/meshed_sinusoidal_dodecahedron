function [V, T] = meshed_sinusoidal_dodecahedron(sampling, w, option_display)
%% meshed_sinusoidal_dodecahedron : function to compute, mesh, display
% and save a meshed sinusoidal dodecaahedron.
%
% Author & support : nicolas.douillet (at) free.fr 2020.
%
%
% Syntax
%
% meshed_sinusoidal_dodecahedron;
% meshed_sinusoidal_dodecahedron(sampling);
% meshed_sinusoidal_dodecahedron(sampling, w);
% meshed_sinusoidal_dodecahedron(sampling, w, option_display);
% [V, T] = meshed_sinusoidal_dodecahedron(sampling, w, option_display);
%
%
% Description
%
% meshed_sinusoidal_dodecahedron computes and displays a meshed sinusoidal
% icosahedron with parameters sampling = 60, w = 1, and option_display =
% true, by default.
%
% meshed_sinusoidal_dodecahedron(sampling) samples at the value sampling
% each icosahedron basis triangle (20).
%
% meshed_sinusoidal_dodecahedron(sampling, w) uses the given shape parameter
% w.
%
% meshed_sinusoidal_dodecahedron(sampling, w, option_display) displays the
% result when option_display = *true/*1 and doesn't when option_display =
% false/0.
%
% [V, T] = meshed_sinusoidal_dodecahedron(sampling, w, option_display)
% stores the point set coordinates and its corresponding triangulation
% respectively in V and T.
%
%
% See also MESH, TRIMESH, SPHERE
%
%
% Input arguments
%
% - sampling : positive integer scalar double, sampling > 2. Remarkable value : sampling = 3 gives an icosahedron
%   
% - w : real scalar double, the shape parameter. Remarkable value : w = 0 gives a geoid.
%
% - option_display : logical *true (1) / false (0).
%
%
% Output arguments
%
%       [|  |  | ]
% - V = [Vx Vy Vz], real matrix double, the point set. Size(V) = [nb_vertices,3].
%       [|  |  | ]
%
%       [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%       [ |  |  |]
%
%
% Example #1 : default parameters values
% meshed_sinusoidal_dodecahedron;
%
%
% Example #2 : minimum sampling step
% meshed_sinusoidal_dodecahedron(6);
%
%
% Example #3 : negative shape parameter value
% meshed_sinusoidal_dodecahedron(60,-1);
%
%
% Example #4 : large shape parameter value
% meshed_sinusoidal_dodecahedron(60,3);


%% Input parsing
assert(nargin < 4, 'Too many input arguments.');

if nargin < 3
    option_display = true;
    if nargin < 2
        w = 1;
        if nargin < 1
            sampling = 60;
        else
            assert(sampling == floor(sampling) && sampling > 2,'sampling parameter must be an integer with a value greater or equal to 3.');
        end
    else
        assert(isnumeric(w),'w parameter must be of class numeric.'); 
    end
else
    assert(islogical(option_display) || isnumeric(option_display),'option_display parameter class must be either logical or numeric.');
end


%% Body
phi_n = 0.5*(1+sqrt(5));

Mrz = [cos(0.4*pi) -sin(0.4*pi) 0;...
       sin(0.4*pi) cos(0.4*pi) 0;...
       0 0 1];

centre_angle = 2*asin(1/sqrt(phi_n*sqrt(5)));
a = 4/3/sqrt(phi_n*sqrt(5)); % edge length

% Icosahedron vertices coordinates
% 1st equilateral triangle
V0 = [0 0 1]';
V1 = [sin(centre_angle) 0 cos(centre_angle)]';
V2 = Mrz*V1;

% Lower base triangle with /O symetry
V3 = -V0;
V4 = -V1;
V5 = -V2;

% (12) vertices set coordinates vector
U0 = Mrz*V2;
U1 = Mrz^2*V2;
U2 = Mrz^3*V2;
U3 = Mrz*V5;
U4 = Mrz^2*V5;
U5 = Mrz^3*V5;

% Icosahedron face centres
C0 = mean(cat(2,V0,V1,V2),2)';
C1 = mean(cat(2,V0,V2,U0),2)';
C2 = mean(cat(2,V0,U0,U1),2)';
C3 = mean(cat(2,V0,U1,U2),2)';
C4 = mean(cat(2,V0,U2,V1),2)';

C5 = mean(cat(2,V3,V4,V5),2)';
C6 = mean(cat(2,V3,V5,U3),2)';
C7 = mean(cat(2,V3,U3,U4),2)';
C8 = mean(cat(2,V3,U4,U5),2)';
C9 = mean(cat(2,V3,U5,V4),2)';

C10 = mean(cat(2,V1,V2,U4),2)';
C11 = mean(cat(2,V2,U4,U5),2)';
C12 = mean(cat(2,V2,U0,U5),2)';
C13 = mean(cat(2,U0,U5,V4),2)';
C14 = mean(cat(2,U0,U1,V4),2)';

C15 = mean(cat(2,U1,V4,V5),2)';
C16 = mean(cat(2,U1,U2,V5),2)';
C17 = mean(cat(2,U2,V5,U3),2)';
C18 = mean(cat(2,U2,V1,U3),2)';
C19 = mean(cat(2,V1,U3,U4),2)';

C = cat(1,C0,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19)';

% Bidirectional (u,v) sampling + compute corresponding squared distances vector
[T0, I0] = sample_triangle(V0, V1, V2, floor(sampling/3));
T0 = T0';

% Replicate / rotation -> upper crown
Tu = T0;
Iu = I0;

for k = 1:4    
    Tu = cat(2,Tu, Mrz^k*T0);
    Iu = cat(1,Iu,I0+max(max(Iu)));
end

% Lower base triangle with /O symetry
V3 = -V0;
V4 = -V1;
V5 = -V2;

[T1, I1] = sample_triangle(V3, V4, V5, floor(sampling/3));
T1 = T1';

% Lower crown
Tl = T1;
Il = I1;

for k = 1:4
    Tl = cat(2,Tl,Mrz^k*T1);
    Il = cat(1,Il,I1+max(max(Il)));
end

T = [Tu Tl];
I = [Iu; Il+max(max(Iu))];

% 1st belt triangle
V6 = V1;
V7 = V2;
V8 = Mrz^2*V5;

[T2, I2] = sample_triangle(V6, V8, V7, floor(sampling/3));
T2 = T2';

% 2nd belt triangle
V9  = -V6;
V10 = -V7;
V11 = -V8;

[T3, I3] = sample_triangle(V9, V10, V11, floor(sampling/3));
T3 = T3';

% Full belt = centre crown
T4 = [T2 T3];
I4 = [I2; I3+max(max(I2))];
Tc = T4;
Ic = I4;

for k = 1:4
    Tc = cat(2,Tc,Mrz^k*T4);
    Ic = cat(1,Ic,I4+max(max(Ic)));
end

Sinico = [T Tc]; 
every_idx = [I; Ic+max(max(I))];

% Spherical coordinates
X = Sinico(1,:);
Y = Sinico(2,:);
Z = Sinico(3,:);

X = X(:);
Y = Y(:);
Z = Z(:);
V = [X Y Z];

[V, every_idx] = remove_duplicated_vertices(V, every_idx);

for i = 1:size(V,1)
    
    radius = sqrt(sum(V(i,:).^2));
    
    if radius ~= 1
        V(i,:) = V(i,:) / radius;
    end
end

f = [];
coeff = [];

for i = 1:size(V,1)
    
    % Closest vertex
    [~, min_dst] = closest_vertex(V(i,:)',C);        
    
    if min_dst < 0.5*a
        f = cat(1,f,i);
        coeff = cat(1,coeff,2*pi*min_dst/a);
    end
    
end

N = V(f,:);

X_s = N(:,1);
Y_s = N(:,2);
Z_s = N(:,3);

radius = sqrt(sum(X_s.^2+Y_s.^2+Z_s.^2, 2));

theta_c = acos(Z_s ./ radius);
phi_c = atan2(Y_s, X_s);

Rho_s = w*0.5*(1+cos(coeff));

% Arcos option
X_s = X_s + Rho_s.*sin(theta_c).*cos(phi_c);
Y_s = Y_s + Rho_s.*sin(theta_c).*sin(phi_c);
Z_s = Z_s + Rho_s.*cos(theta_c);

N = [X_s Y_s Z_s];
V(f,:) = N;

[V, T] = remove_duplicated_vertices(V, every_idx);
TRI = triangulation(T, V(:,1), V(:,2), V(:,3));


% Display
if option_display 
    
    figure;
    set(gcf,'Color',[1 1 1]), set(gca,'Color',[1 1 1]);    
    trisurf(TRI), hold on;    
    colormap([0 1 1]);
    axis square, axis equal;
    grid on;
    view(3);
    
end


end % meshed_sinusoidal_dodecahedron


%% Closest_vertex subfunction
function [V0, min_dst] = closest_vertex(Point,V)


dst_vect = sqrt(sum((V-repmat(Point, [1,size(V,2)])).^2,1));

f = find(dst_vect == min(dst_vect));
min_dst = dst_vect(f(1));
V0 = V(:,f(1));


end % closest_vertex


%% Remove_duplicated_vertices subfunction
function [V_out, T_out] = remove_duplicated_vertices(V_in, T_in)


tol = 1e4*eps;
[V_out,~,n] = uniquetol(V_in,tol,'ByRows',true);
T_out = n(T_in);


end


%% sample_triangle subfunction
function [V, T] = sample_triangle(V1, V2, V3, nbstep)


% Create sampling grid
global Ndim;

Ndim = size(V1,1);

% (V1V2, V1V3) base
u = (V2 - V1);
v = (V3 - V1);

V = zeros(sum(1:nbstep+1),Ndim);

nu = u / norm(u);
nv = v / norm(v);
stepu = norm(u) / nbstep;
stepv = norm(v) / nbstep;
k = 1;

% Sampling & vertices generation
for m = 0:nbstep
    
    for n = 0:nbstep
        
        if (m+n) <= nbstep % in (V1,V2,V3) triangle conditions ; indices # nb segments
            
            % translation vector
            tv = m*stepu*nu + n*stepv*nv;
            V(k,:) = (V1 + tv)';
            k = k+1;
            
        end
        
    end
    
end


% Index triplets list construction
T = zeros(nbstep^2,3);
row_length = 1 + nbstep;
cum_row_length = row_length;
row_idx = 1;
p = 1;

while (p <= nbstep^2) && (row_length > 1)
    
     i = p;
    
    if p < 2 % "right" triangle serie only
        
        while i < cum_row_length
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;
            i = i +1;
            
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    else
        
        % Since p >= 2
        while i < cum_row_length % both triangle series
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;            
            T(row_idx,:) = [i i-row_length i+1]; % + upside-down triangles serie
            row_idx = row_idx + 1;            
            i = i +1;
            
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    end
    
end

T = sort(T,2);
T = unique(T,'rows','stable');


end % sample_triangle