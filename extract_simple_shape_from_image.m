clear all
close all
clc

Image = imread('mod6.bmp');
dy1 = [diff(Image, 1, 1); zeros(1, size(Image,2))];
dy2 = [dy1(:,2:end) zeros(size(Image,1), 1)];
dx1 = [diff(Image, 1, 2) zeros(size(Image,1), 1)];
dx2 = [dx1(2:end,:); zeros(1, size(Image,2))];
edgeMap = (dx1 ~= 0) | (dy1 ~= 0) | (dx2 ~= 0) | (dy2 ~= 0);

[edgeYs, edgeXs] = find(edgeMap ~= 0);
ex = edgeXs(1);
ey = edgeYs(1);
edgeXs = ex;
edgeYs = ey;
% Classify the pixels incident on this vertex as in the object or not
incidentObjPixels = zeros(5,1);
if (Image(ey, ex) == 0)
  incidentObjPixels(1) = 1;
endif
if (Image(ey+1, ex) == 0)
  incidentObjPixels(2) = 1;
endif
if (Image(ey+1, ex+1) == 0)
  incidentObjPixels(3) = 1;
endif
if (Image(ey, ex+1) == 0)
  incidentObjPixels(4) = 1;
endif
incidentObjPixels(5) = incidentObjPixels(1); % wrap around
idx = find(diff(incidentObjPixels) == -1, 1, 'first');
nextEx = ex;
nextEy = ey;
if (idx == 1)
  nextEx = ex - 1;
elseif (idx == 2)
  nextEy = ey + 1;
elseif (idx == 3)
  nextEx = ex + 1;
else % idx == 4
  nextEy = ey - 1;
endif

startX = ex;
startY = ey;
while (nextEx ~= startX  || nextEy ~= startY)
  edgeXs = [edgeXs; nextEx];
  edgeYs = [edgeYs; nextEy];
  edgeMap(nextEy, nextEx) = 0;
  if (edgeMap(nextEy, nextEx - 1) ~= 0 && ...
    Image(nextEy, nextEx) ~= Image(nextEy + 1, nextEx) && ...
    (nextEx - 1 ~= ex || nextEy != ey))
    ex = nextEx;
    ey = nextEy;
    nextEx = nextEx - 1;
    nextEy = nextEy;
  elseif (edgeMap(nextEy + 1, nextEx) ~= 0 && ...
    Image(nextEy + 1, nextEx) ~= Image(nextEy + 1, nextEx + 1) && ...
    (nextEx ~= ex || nextEy + 1 != ey))
    ex = nextEx;
    ey = nextEy;
    nextEx = nextEx;
    nextEy = nextEy + 1;
  elseif (edgeMap(nextEy, nextEx + 1) ~= 0 && ...
    Image(nextEy, nextEx + 1) ~= Image(nextEy + 1, nextEx + 1) && ...
    (nextEx + 1 ~= ex || nextEy != ey))
    ex = nextEx;
    ey = nextEy;
    nextEx = nextEx + 1;
    nextEy = nextEy;
  elseif (edgeMap(nextEy - 1, nextEx) ~= 0 && ...
    Image(nextEy, nextEx) ~= Image(nextEy, nextEx + 1) && ...
    (nextEx ~= ex || nextEy - 1 != ey))
    ex = nextEx;
    ey = nextEy;
    nextEx = nextEx;
    nextEy = nextEy - 1;
  else
    break;
  endif
end

edgeYs = edgeYs + 0.5;
edgeXs = edgeXs + 0.5;

imshow(Image);
hold on
plot(edgeXs, edgeYs, 'r-', 'linewidth', 2);