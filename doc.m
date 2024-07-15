%% meshed_sinusoidal_dodecahedron
%
% Function to compute, mesh, display and
% save a meshed sinusoidal dodecahedron.
%
% Author : nicolas.douillet9 (at) gmail.com, 2016-2024.
%
%% Syntax
%
% meshed_sinusoidal_dodecahedron;
%
% meshed_sinusoidal_dodecahedron(sampling);
%
% meshed_sinusoidal_dodecahedron(sampling, w);
%
% meshed_sinusoidal_dodecahedron(sampling, w, option_display);
%
% [V, T] = meshed_sinusoidal_dodecahedron(sampling, w, option_display);
%
%% Description
%
% meshed_sinusoidal_dodecahedron computes and displays a meshed sinusoidal
% dodecahedron with parameters sampling = 60, w = 1, and option_display =
% true, by default.
%
% meshed_sinusoidal_dodecahedron(sampling) samples at the value sampling
% each dodecahedron basis triangle (20).
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
%% See also
%
% <https://fr.mathworks.com/help/matlab/ref/sphere.html sphere> |
% <https://fr.mathworks.com/help/matlab/ref/mesh.html mesh> |
% <https://fr.mathworks.com/help/matlab/ref/trimesh.html trimesh> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73634-meshed-sinusoidal-icosahedron meshed_sinusoidal_icosahedron> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73159-meshed-reuleaux-tetrahedron meshed_reuleaux_tetrahedron> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/69212-n-level-geoid n_level_geoid>
%
%% Input arguments
%
% - sampling : positive integer scalar double, sampling > 2. Remarkable value : sampling = 3 gives an dodecahedron
%   
% - w : real scalar double, the shape parameter. Remarkable value : w = 0 gives a geoid.
%
% - option_display : logical *true (1) / false (0).
%
%% Output arguments
%
%        [|  |  | ]
% - V = [Vx Vy Vz], real matrix double, the point set. Size(V) = [nb_vertices,3].
%        [|  |  | ]
%
%        [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%        [ |  |  |]
%
%% Example #1
% Default parameters values
meshed_sinusoidal_dodecahedron;

%% Example #2
% Minimum sampling step
meshed_sinusoidal_dodecahedron(6);

%% Example #3
% Negative shape parameter value
meshed_sinusoidal_dodecahedron(60,-1);

%% Example #4
% Large shape parameter value
meshed_sinusoidal_dodecahedron(60,3);