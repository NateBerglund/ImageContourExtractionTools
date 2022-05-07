clear all
close all
clc
%pkg load geometry
%pkg load msh

extrusion_height = 3.0;

polygon = csvread('border_polygon.csv');
polygon(:,2) = 500-polygon(:,2); % compensate for the original being in image coords

fid = fopen('extruded_polygon.stl','wt');
fprintf(fid, 'solid extruded_polygon\n');

n = size(polygon,1);
vertices = zeros(6*n,3);
vertices(1:3:3*n,1:2) = polygon;
vertices(1:3:3*n,3) = 0;
vertices(2:3:3*n,1:2) = [polygon(end,:); polygon(1:end-1,:)];
vertices(2:3:3*n,3) = 0;
vertices(3:3:3*n,1:2) = [polygon(end,:); polygon(1:end-1,:)];
vertices(3:3:3*n,3) = extrusion_height;
vertices(3*n+1:3:6*n,1:2) = polygon;
vertices(3*n+1:3:6*n,3) = extrusion_height;
vertices(3*n+2:3:6*n,1:2) = polygon;
vertices(3*n+2:3:6*n,3) = 0;
vertices(3*n+3:3:6*n,1:2) = [polygon(end,:); polygon(1:end-1,:)];
vertices(3*n+3:3:6*n,3) = extrusion_height;
    
for f = 1:(size(vertices,1)/3)
  fprintf(fid, 'facet normal 0.0 0.0 1.0\n');
  fprintf(fid, '    outer loop\n');
  fprintf(fid, '        vertex %.3f %.3f %.3f\n', vertices(3*f-2,:));
  fprintf(fid, '        vertex %.3f %.3f %.3f\n', vertices(3*f-1,:));
  fprintf(fid, '        vertex %.3f %.3f %.3f\n', vertices(3*f,:));
  fprintf(fid, '    endloop\n');
  fprintf(fid, 'endfacet\n');
endfor

% TODO: Figure out the hard part, which is to triangulate the polygon's interior

% Keep around until you implement, as this may be needed for prettier visualization:
%% Orientation flip
%temp = vertices(1:3:end,:);
%vertices(1:3:end,:) = vertices(2:3:end,:);
%vertices(2:3:end,:) = temp;

fprintf(fid, 'endsolid extruded_polygon\n');
fclose(fid)

% The following is a temporary workaround using the gmsh application:

%filename = "mesh";
%meshsize = sqrt(mean(sumsq(diff(polygon, 1, 1), 2)))/2;
%data2geo(polygon, meshsize, "output", [filename ".geo"]);
% Note: The msh package's 'msh2m_gmsh' function appears to be broken, so instead
% you will need to run gmsh manually with a command similar to:
% gmsh -format msh -2 -o mesh.msh mesh.geo 2>&1
%
% Then you'll need to open mesh.msh in gmsh and export it as an stl.
%
% Finally, make a copy of the stl and edit it in Notepad++, doing a find and
% replace of " 0\r\n" with " 3\r\n" (Extended mode) -- make the '3' in the
% example whatever extrusion_height is. Be patient as the file may be millions
% of lines of text long and it may take Notepad++ a long time to do this, but
% it will eventually finish.
%
% Finally, manually combine extruded_polygon.stl, mesh.stl, mesh_copy.stl in
% a text editor such as Notepad++, to produce the final full stil file.
%
% Note: The file size produced this way is much larger than it needs to be, so
% I really want to implement a better polygon triangulation algorithm in
% Matlab/Octave (or find one I can modify/re-use).

% Display the plot
hold on
axis equal
plot([polygon(:,1);polygon(1,1)], [polygon(:,2);polygon(1,2)], 'r-', 'linewidth', 2);
