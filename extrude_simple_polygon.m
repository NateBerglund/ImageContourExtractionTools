clear all
close all
clc

extrusion_height = 1.0;

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

%% Orientation flip
%temp = vertices(1:3:end,:);
%vertices(1:3:end,:) = vertices(2:3:end,:);
%vertices(2:3:end,:) = temp;

fclose(fid)

% Display the plot
hold on
axis equal
plot([polygon(:,1);polygon(1,1)], [polygon(:,2);polygon(1,2)], 'r-', 'linewidth', 2);
