clear all
close all
clc

epsilon = 0.05;

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

% Simplify runs of vertices that are all in a straight line
inALine = false;
startIdx = 0;
endIdx = 0;
idx = 3;
while (idx <= numel(edgeXs) || inALine) % we check inALine to allow "wrap around"
  
  % adjust idx as necessary so it is within bounds
  wrappedIdx = mod(idx - 1, numel(edgeXs)) + 1;
  
  % will this current point be in a line of at least 3 points?
  willBeInALine = ...
    (edgeXs(wrappedIdx) - edgeXs(wrappedIdx-1) == edgeXs(wrappedIdx-1) - edgeXs(wrappedIdx-2) &&...
     edgeYs(wrappedIdx) - edgeYs(wrappedIdx-1) == edgeYs(wrappedIdx-1) - edgeYs(wrappedIdx-2));
     
  if (~inALine && willBeInALine) % starting a line
    startIdx = idx - 2;
    inALine = true;
  endif
  
  if (inALine && ~willBeInALine) % ending a line
    endIdx = idx - 1;
    inALine = false;
    
    % How many points we need to remove
    countToRemove = endIdx - startIdx - 1;
    
    % adjust endIdx so it is within bounds (do this after computing countToRemove)
    endIdx = mod(endIdx - 1, numel(edgeXs)) + 1;
    
    % Circularly rotate the arrays so the final point of the line will be at index 1.
    edgeXs = [edgeXs(endIdx:end); edgeXs(1:endIdx-1)];
    edgeYs = [edgeYs(endIdx:end); edgeYs(1:endIdx-1)];
    
    % Now we just have to remove the last countToRemove points
    edgeXs = edgeXs(1:end-countToRemove);
    edgeYs = edgeYs(1:end-countToRemove);
    
    idx = 2; % adjust the index (the point we just examined is now at index 2)
  endif
  
  idx = idx + 1; % increment the index
end

% ----- Remove "jaggies" by collapsing any edges of length 1 -----

% First, create doubly-sized arrays with an extra element in between each
% pair of vertices to (potentially) store the new vertex along the edge between
% them.
n = numel(edgeXs);
newEdgeXs = zeros(2 * n, 1);
newEdgeYs = zeros(2 * n, 1);
includeMask = zeros(2 * n, 1, 'logical');
includeMask(1:2:end) = true;
newEdgeXs(1:2:end) = edgeXs;
newEdgeYs(1:2:end) = edgeYs;
% Next, initialize some variables that need to be set before the first
% iteration of the loop
isPriorEdgeExtremal = false;
isCurrentEdgeExtremal = ...
  [edgeXs(end-1)-edgeXs(end-2) edgeYs(end-1)-edgeYs(end-2)] * ...
  [edgeXs(end)-edgeXs(1); edgeYs(end)-edgeYs(1)] > epsilon; 
isNextEdgeExtremal = ...
  [edgeXs(end)-edgeXs(end-1) edgeYs(end)-edgeYs(end-1)] * ...
  [edgeXs(1)-edgeXs(2); edgeYs(1)-edgeYs(2)] > epsilon;
% Main loop to identify edges to collapse
for idx = 1:n
  isPriorEdgeExtremal = isCurrentEdgeExtremal;
  isCurrentEdgeExtremal = isNextEdgeExtremal; 
  isNextEdgeExtremal = ...
      [edgeXs(idx)-edgeXs(mod(idx-2,n)+1) ...
       edgeYs(idx)-edgeYs(mod(idx-2,n)+1)] * ...
      [edgeXs(mod(idx,n)+1)-edgeXs(mod(idx+1,n)+1); ...
       edgeYs(mod(idx,n)+1)-edgeYs(mod(idx+1,n)+1)] > epsilon;
  % Is the edge between the previous and current vertex length 1?
  if (abs(sqrt((edgeXs(idx)-edgeXs(mod(idx-2,n)+1))^2 +
               (edgeYs(idx)-edgeYs(mod(idx-2,n)+1))^2) - 1) ...
               < epsilon && ~isCurrentEdgeExtremal)
    middleIdx = 2*(idx-1); % index of the "middle of the edge" in the doubly-sized arrays
    if (isPriorEdgeExtremal && ~isNextEdgeExtremal) % favor the prior vertex if only the prior edge was extremal
      newEdgeXs(middleIdx) = edgeXs(idx-1);
      newEdgeYs(middleIdx) = edgeYs(idx-1);
    elseif (~isPriorEdgeExtremal && isNextEdgeExtremal) % favor the current vertex if only the next edge is extremal
      newEdgeXs(middleIdx) = edgeXs(idx);
      newEdgeYs(middleIdx) = edgeYs(idx);
    else % Take the average of the prior and current vertices in all other cases
      newEdgeXs(middleIdx) = 0.5 * (edgeXs(idx-1) + edgeXs(idx));
      newEdgeYs(middleIdx) = 0.5 * (edgeYs(idx-1) + edgeYs(idx));
    endif
    % Include the new vertex in the middle but not the vertices at the endpoints
    includeMask(middleIdx-1) = false;
    includeMask(middleIdx) = true;
    includeMask(middleIdx+1) = false;
  endif
end
edgeXs = newEdgeXs(includeMask);
edgeYs = newEdgeYs(includeMask);

% Display the plot
imshow(Image);
hold on
plot([edgeXs; edgeXs(1)] + 0.5, [edgeYs; edgeYs(1)] + 0.5, 'r-', 'linewidth', 2);