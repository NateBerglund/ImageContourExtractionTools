clear all
close all
clc

ppmm = 20;%39.8630; % pixels per mm
epsilon = 0.05;

Image = imread('mug_profile2.bmp');
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
endwhile

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
endwhile

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
endfor
edgeXs = newEdgeXs(includeMask);
edgeYs = newEdgeYs(includeMask);

% ----- Resample the curve about once per pixel of length (will make our Taubin smoothing work better) -----
% Do the resampling in a way that will respect "corners", however.

newEdges = [edgeXs'; edgeYs'];
fixedPoints = ...
  sum(([newEdges(:,2:end) newEdges(:,1)] - newEdges) .* ...
      ([newEdges(:,end) newEdges(:,1:end-1)] - newEdges), 1) > ...
      sqrt(sum(([newEdges(:,2:end) newEdges(:,1)] - newEdges).^2, 1)) .* ...
      sqrt(sum(([newEdges(:,end) newEdges(:,1:end-1)] - newEdges).^2, 1)) * ...
      (-0.8); % angle must be more acute than arccos(-0.8) ~= 143 degrees

n = numel(edgeXs);
polygonLength = 0;
for idx = 1:n
  % Add the length of the edge between the current and next vertex
  polygonLength = polygonLength + ...
    sqrt((edgeXs(mod(idx,n)+1)-edgeXs(idx))^2 + ...
         (edgeYs(mod(idx,n)+1)-edgeYs(idx))^2);
endfor
polygonLength = polygonLength + epsilon;
edgUnit = polygonLength / round(polygonLength); % approximately 1 pixel

fixedPointIndices = find(fixedPoints);
if numel(fixedPointIndices) > 0
  % Circularly rotate the points so that they start at a fixed point
  edgeXs = [edgeXs(fixedPointIndices(1):end); edgeXs(1:fixedPointIndices(1)-1)];
  edgeYs = [edgeYs(fixedPointIndices(1):end); edgeYs(1:fixedPointIndices(1)-1)];
  fixedPointIndices = fixedPointIndices - fixedPointIndices(1) + 1;
else
  fixedPointIndices = 1; % treat the first point as if it were a fixed point
endif
fixedPointIndices = [fixedPointIndices n+1]; % "wrap around"

newEdges = zeros(2,0);
fixedPoints = zeros(1,0,'logical');
for fpIdx = 1:numel(fixedPointIndices)-1
  % Add the fixed-point vertex as-is
  newEdges = [newEdges [edgeXs(fixedPointIndices(fpIdx)); ...
                        edgeYs(fixedPointIndices(fpIdx))]];
  fixedPoints = [fixedPoints; true];
  
  % Delete the previous vertex if this will produce an edge of length < 0.5 * edgUnit
  if (size(newEdges,2) > 1 && ...
    sqrt((newEdges(1,end)-newEdges(1,end-1))^2 + ...
         (newEdges(2,end)-newEdges(2,end-1))^2) < 0.5 * edgUnit)
    newEdges(:,end-1) = [];
    fixedPoints(end-1) = [];
  endif
  
  distanceRemaining = edgUnit;
  for idx = fixedPointIndices(fpIdx):fixedPointIndices(fpIdx+1)-1
    % Add the length of the edge between the current and next vertex
    edgeLen = sqrt((edgeXs(mod(idx,n)+1)-edgeXs(idx))^2 + ...
                   (edgeYs(mod(idx,n)+1)-edgeYs(idx))^2);
    edgeUnitVec = [edgeXs(mod(idx,n)+1)-edgeXs(idx); ...
                   edgeYs(mod(idx,n)+1)-edgeYs(idx)] / edgeLen;
    
    currentVertex = [edgeXs(idx); edgeYs(idx)];
    edgeLenOffset = 0; % offset along the original edge
    while (edgeLenOffset < edgeLen)
      if (distanceRemaining > 0)
        if (distanceRemaining < edgeLen)
          edgeLenOffset = distanceRemaining;
          newEdges = [newEdges currentVertex + edgeLenOffset * edgeUnitVec];
          fixedPoints = [fixedPoints; false];
          distanceRemaining = 0;
        else
          edgeLenOffset = edgeLen;
          distanceRemaining = distanceRemaining - edgeLen;
          break;
        endif
      endif
      if (edgeLenOffset + edgUnit < edgeLen)
        edgeLenOffset = edgeLenOffset + edgUnit;   
        newEdges = [newEdges currentVertex + edgeLenOffset * edgeUnitVec];
        fixedPoints = [fixedPoints; false];
      else
        distanceRemaining = edgeLenOffset + edgUnit - edgeLen;
        break;
      endif
    endwhile
  endfor
endfor

%debugging
edgUnit
maxDeviation = max(abs(sqrt(sum(diff(newEdges,1,2).^2,1))-edgUnit))
maxEdgeLen = max(sqrt(sum(diff(newEdges,1,2).^2,1)))
minEdgeLen = min(sqrt(sum(diff(newEdges,1,2).^2,1)))

% ----- Run Taubin smoothing, but fixing vertices likely to be corners -----

lambda = 0.5;
mu = -0.51;
n = size(newEdges, 2);
coeffs = [lambda; mu];
for iter = 1:100
  convolvedEdges = (1/3)*newEdges + ...
                   (1/3)*[newEdges(:,2:end) newEdges(:,1)] + ...
                   (1/3)*[newEdges(:,end) newEdges(:,1:end-1)];
  displacements = convolvedEdges - newEdges;
  coeff = coeffs(mod(iter-1,2)+1);
  newEdges(:,~fixedPoints) = newEdges(:,~fixedPoints) + coeff * displacements(:,~fixedPoints);
endfor

%debugging
disp('after smoothing')
edgUnit
maxDeviation = max(abs(sqrt(sum(diff(newEdges,1,2).^2,1))-edgUnit))
maxEdgeLen = max(sqrt(sum(diff(newEdges,1,2).^2,1)))
minEdgeLen = min(sqrt(sum(diff(newEdges,1,2).^2,1)))

% Simplify edges likely to be straight lines
run_distance = 30;
straight_line_squared_tolerance = 0.05;
forward_diffs = (circshift(newEdges, [0,-run_distance]) - newEdges) / run_distance;
newEdgesExt = [newEdges newEdges(:,1:run_distance)];
isLineMask = zeros(1, size(newEdges,2), 'logical');
for i=1:size(newEdges,2)
  predicted = repmat(newEdges(:,i), 1, run_distance-1) + forward_diffs(:,i)*(1:1:run_distance-1);
  actual = newEdgesExt(:,i+1:i+run_distance-1);
  diffssq = sum((actual-predicted).^2,1);
  if (all(diffssq < straight_line_squared_tolerance))
    %isLineMask(i+1:i+run_distance-1) = true;
    isLineMask(i+1) = true;   
  endif
endfor
for i=size(newEdges,2)+1:-1:run_distance
  if (isLineMask(i+1-run_distance) && ~isLineMask(i+2-run_distance))
    isLineMask(i+1-run_distance:i-1) = true;
    if (i <= size(newEdges,2))
      isLineMask(i) = false;
    endif
  endif
endfor
for i=1:size(newEdges,2)
  if (fixedPoints(i))
    isLineMask(mod(i+1,numel(isLineMask))+1) = true;
    isLineMask(mod(i,numel(isLineMask))+1) = true;
    isLineMask(i) = false;
    isLineMask(mod(i+numel(isLineMask)-2,numel(isLineMask))+1) = true;
    isLineMask(mod(i+numel(isLineMask)-3,numel(isLineMask))+1) = true;
  endif
endfor
isLineMask

detected_line_starts = [];
detected_line_ends = [];
transitions = diff([isLineMask isLineMask(1)]);
transitions = [transitions transitions(1:run_distance+2)];
for i=1:numel(isLineMask)-1
  if (transitions(i) > 0 && all(transitions(i+1:i+run_distance-1)==0))
    detected_line_starts = [detected_line_starts; i+1];
  endif
  if (transitions(i+run_distance) < 0 && all(transitions(i+1:i+run_distance-1)==0))
    detected_line_ends = [detected_line_ends; i+run_distance];
  endif
endfor
if (numel(detected_line_ends) < numel(detected_line_starts))
  detected_line_ends = [detected_line_ends; numel(isLineMask)];
endif

detected_line_starts
detected_line_ends

isLineMask = zeros(1, size(newEdges,2), 'logical'); 
for i=1:numel(detected_line_starts)
  isLineMask(detected_line_starts(i):detected_line_ends(i)) = true;
endfor

newEdges(:,isLineMask) = [];
fixedPoints(isLineMask) = [];

fixedPointIndices = find(fixedPoints)


% Display the plot
imshow(Image);
hold on
plot([newEdges(1,:) newEdges(1,1)] + 0.5, [newEdges(2,:) newEdges(2,1)] + 0.5, 'r-', 'linewidth', 2);
plot(newEdges(1,fixedPoints) + 0.5, newEdges(2,fixedPoints) + 0.5, 'bo', 'markersize', 10);
plot(newEdges(1,:) + 0.5, newEdges(2,:) + 0.5, 'g+', 'markersize', 2);

csvwrite('mug_border_polygon.csv', newEdges' / ppmm);
