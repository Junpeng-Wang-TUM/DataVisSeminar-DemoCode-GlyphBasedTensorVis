%% This code is created for visualizing the 2D stress tensor field via the Elliptic Tensor glyphs
%% Author: Junpeng Wang (junpeng.wang@tum.de)
%% Date: 2021-09-14
%% Refer to Fig.2a in paper "Stress Topology Analysis for Porous Infill Optimization" by 
%%	Wang, J., Wu, J. and Westermann, R. (2021, arXiv:2108.09675) for detailed explanations.

function EllipticTensorGlyph2D(stressfileName,samplingSpace,scalingGlyph)
	%%1. Read Data
	LoadStressField(stressfileName);
	
	%%2. Setup Sampling Points
	SetupSamplingPoints(samplingSpace);

	%%3. Create Glyphs
	figure(1); CreateGlyphics(scalingGlyph*50);
end

function SetupSamplingPoints(samplingSpace)
	global nelx_;
	global nely_;
	global carEleMapBack_;
	global eNodMat_;
	global nodeCoords_;
	global cartesianStressField_;
	global samplingPoints_;
	global principalStressList_;
	
	q = samplingSpace;
	validElementsVolume = zeros(nelx_*nely_,1);
	validElementsVolume(carEleMapBack_) = (1:length(carEleMapBack_))';
	validElementsVolume = reshape(validElementsVolume, nely_, nelx_);
	sampledElements = validElementsVolume(q:q:nely_-q, q:q:nelx_-q);	
	%sampledElements = validElementsVolume([1 q:q:nely_], [1 q:q:nelx_]);
	sampledElements = reshape(sampledElements, numel(sampledElements), 1);
	sampledElements(0==sampledElements) = [];
	
	numSampledElements = size(sampledElements,1);
	samplingPoints_ = zeros(numSampledElements,2);
	stressTensorElementCentroid = zeros(numSampledElements,3);
	shapeFuncs = ShapeFunction(0, 0);
	for ii=1:numSampledElements
		iEle = sampledElements(ii);
		iNodes = eNodMat_(iEle,:);
		iCoords = nodeCoords_(iNodes, :);
		cartesianStress = cartesianStressField_(iNodes, :);
		samplingPoints_(ii,:) = shapeFuncs*iCoords;
		stressTensorElementCentroid(ii,:) = shapeFuncs*cartesianStress;
	end
	principalStressList_ = ComputePrincipalStress(stressTensorElementCentroid);		
end

function CreateGlyphics(scalingFactor)
	global numEles_;
	global eNodMat_;
	global nodeCoords_;
	global nodesOnBoundary_;
	global samplingPoints_;
	global principalStressList_;
	
    lw = 2.0;
	tmp = abs(principalStressList_(:,[1 4])); maxEigen = max(tmp(:));
	numTensors = size(principalStressList_,1);
	for ii=1:numTensors
		iTensor = principalStressList_(ii,:);
		iCoord = samplingPoints_(ii,:);
		D = [iTensor(1) 0; 0 iTensor(4)];
		V = iTensor([2 3 5 6])'; V = reshape(V,2,2);	
		eigen1 = D(1,1);
		eigen2 = D(2,2);
		if abs(eigen1) < abs(eigen2)
			tmp = [V(1,2) V(2,2)];
			V(1,2) = V(1,1); V(2,2) = V(2,1);
			V(1,1) = tmp(1);
			V(2,1) = tmp(2);
			D(2,2) = eigen1;
			D(1,1) = eigen2;
		end
		hold on;
		dScale = 1.6*abs(D(1,1))/maxEigen; dScale = scalingFactor*dScale;
		%dScale = sqrt(dScale);
		xx = [iCoord(1)-dScale*0.5*V(1,1), iCoord(1)+dScale*0.5*V(1,1)];
		yy = [iCoord(2)-dScale*0.5*V(2,1), iCoord(2)+dScale*0.5*V(2,1)];		
		if D(1,1) > 0
			plot(xx, yy, 'Color', [252 141 98]/255, 'linewidth', lw);
		else
			plot(xx, yy, 'Color', [102 194 165]/255, 'linewidth', lw);
		end
		dScale = 1.6*abs(D(2,2))/maxEigen; dScale = scalingFactor*dScale;
		%dScale = sqrt(dScale);
		xx = [iCoord(1)-dScale*0.5*V(1,2), iCoord(1)+dScale*0.5*V(1,2)];
		yy = [iCoord(2)-dScale*0.5*V(2,2), iCoord(2)+dScale*0.5*V(2,2)];
		if D(2,2) > 0
			plot(xx, yy, 'Color', [252 141 98]/255, 'linewidth', lw);
		else
			plot(xx, yy, 'Color', [102 194 165]/255, 'linewidth', lw);
		end

		dScale = 1.6*abs(D(1,1))/maxEigen; dScale = scalingFactor*dScale;
		%dScale = sqrt(dScale);
		x1 = iCoord(1)-dScale*0.5*V(1,1);
		x2 = iCoord(1)+dScale*0.5*V(1,1);
		y1 = iCoord(2)-dScale*0.5*V(2,1);
		y2 = iCoord(2)+dScale*0.5*V(2,1);
		e = sqrt(1-D(2,2)*D(2,2) / D(1,1)/D(1,1));
		[x,y] = ellipseShape(x1,x2,y1,y2,e);
		plot(x,y,'k', 'linewidth', lw); 		
	end

	%%Show silhouette
	edgeIndices = eNodMat_(:, [1 2 2 1  2 3 3 2  3 4 4 3  4 1 1 4])';
	edgeIndices = reshape(edgeIndices(:), 4, 4*numEles_);	
	tmp = zeros(size(nodeCoords_,1),1); tmp(nodesOnBoundary_) = 1;
	tmp = tmp(edgeIndices); tmp = sum(tmp,1);
	boundaryEleEdges = edgeIndices(:,find(4==tmp));
	xPatchs = nodeCoords_(:,1); xPatchs = xPatchs(boundaryEleEdges);
	yPatchs = nodeCoords_(:,2); yPatchs = yPatchs(boundaryEleEdges);		
	cPatchs = zeros(size(yPatchs));
	hold on; hd = patch(xPatchs, yPatchs, cPatchs);
	set(hd, 'FaceColor', 'None', 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 2);
	
	axis equal; axis tight; axis off;
end

function [x,y] = ellipseShape(x1,x2,y1,y2,e)
	a = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
	b = a*sqrt(1-e^2);
	t = linspace(0,2*pi);
	X = a*cos(t);
	Y = b*sin(t);
	w = atan2(y2-y1,x2-x1);
	x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
	y = (y1+y2)/2 + X*sin(w) + Y*cos(w);
end

function principalStress = ComputePrincipalStress(cartesianStress)
	principalStress = zeros(size(cartesianStress,1), 1+2+1+2);
	iPS = zeros(1, 6);
	for ii=1:size(cartesianStress,1)
		iCartesianStress = cartesianStress(ii,:);
		A = iCartesianStress([1 3; 3 2]);
		[eigenVec, eigenVal] = eig(A);
		iPS([1 4]) = diag(eigenVal);
		iPS([2 3 5 6]) = reshape(eigenVec,1,4);
		principalStress(ii,:) = iPS;
	end		
end


function N = ShapeFunction(s, t)
	%				   	   __s (parametric coordinate system)
	%				  	 /-t
	%				*4			*3
	%			*1			*2
	%
	%				nodes
	s = s(:);
	t = t(:);
	N = zeros(size(s,1), 4);
	N(:,1) = 0.25*(1-s).*(1-t);
	N(:,2) = 0.25*(1+s).*(1-t);
	N(:,3) = 0.25*(1+s).*(1+t);
	N(:,4) = 0.25*(1-s).*(1+t);	
end

function LoadStressField(fileName)
	global boundingBox_;
	global numNodes_;
	global nodeCoords_;
	global numEles_;
	global eNodMat_;
	global eleState_;
	global nodState_;
	global elementsOnBoundary_;
	global nodesOnBoundary_;
	global cartesianStressField_;
	global loadingCond_; 
	global fixingCond_;
	global eleCentroidList_;

	global nelx_;
	global nely_;
	global carNodMapForward_;
	global voxelizedVolume_;	
	
	%%read mesh and cartesian stress field
	fid = fopen(fileName, 'r');
	fgetl(fid); fgetl(fid); fgetl(fid); 
	tmp = fscanf(fid, '%s', 1);	
	meshType = fscanf(fid, '%s', 1);
	
	tmp = fscanf(fid, '%s', 1);
	tmp = fscanf(fid, '%d %d', [1 2]);
	nelx_ = tmp(1); nely_ = tmp(2); 
	tmp = fscanf(fid, '%s', 1);
	lBound = fscanf(fid, '%f %f', [1 2]);
	tmp = fscanf(fid, '%s', 1);
	uBound = fscanf(fid, '%f %f', [1 2]);
	boundingBox_ = [lBound; uBound];
	tmp = fscanf(fid, '%s', 1); 
	numValidEles = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%s', 1);
	validElements = fscanf(fid, '%d', [1, numValidEles])';
	validElements = validElements + 1;		

	%%read cartesian stress field
	tmp = fscanf(fid, '%s %s %s %s', 4);
	numStressFields = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%s %s', 2); numLoadedNodes = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%d %f %f', [3, numLoadedNodes]); 
	tmp(1,:) = tmp(1,:)+1; 
	loadingCond_ = tmp';
	tmp = fscanf(fid, '%s %s', 2); numFixedNodes = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%d', [1, numFixedNodes]); 
	fixingCond_ = tmp'+1;
	tmp = fscanf(fid, '%s %s', 2); numValidNods = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%f %f %f', [3, numValidNods]);
	cartesianStressField_ = tmp';
	fclose(fid);

	%%recover cartesian mesh
	voxelizedVolume_ = zeros(nelx_*nely_,1);
	voxelizedVolume_(validElements) = 1;
	voxelizedVolume_ = reshape(voxelizedVolume_, nely_, nelx_);
	RecoverCartesianMesh();
	numNod2ElesVec = zeros(numNodes_,1);
	for ii=1:numEles_
		iNodes = eNodMat_(ii,:);
		numNod2ElesVec(iNodes,1) = numNod2ElesVec(iNodes,1)+1;
	end
	nodesOnBoundary_ = find(numNod2ElesVec<4);	
	nodState_ = zeros(numNodes_,1); nodState_(nodesOnBoundary_) = 1;
	eleState_ = 4*ones(1, numEles_);	
	%%%%
	allNodes = zeros(numNodes_,1);
	allNodes(nodesOnBoundary_) = 1;	
	tmp = zeros(numEles_,1);
	for ii=1:4
		tmp = tmp + allNodes(eNodMat_(:,ii));
	end
	elementsOnBoundary_ = find(tmp>0);				
	%%%%
	%% element centroids
	eleNodCoordListX = nodeCoords_(:,1); eleNodCoordListX = eleNodCoordListX(eNodMat_);
	eleNodCoordListY = nodeCoords_(:,2); eleNodCoordListY = eleNodCoordListY(eNodMat_);
	eleCentroidList_ = [sum(eleNodCoordListX,2) sum(eleNodCoordListY,2)]/4;
	
	loadingCond_(:,1) = carNodMapForward_(loadingCond_(:,1));
	fixingCond_ = carNodMapForward_(fixingCond_);
end

function RecoverCartesianMesh()	
	global nelx_; global nely_;
	global voxelizedVolume_;
	global boundingBox_;
	global numEles_; global numNodes_; global eleSize_;
	global nodeCoords_; global eNodMat_; 
	global carEleMapBack_; global carEleMapForward_;
	global carNodMapBack_; global carNodMapForward_;
	%    __ x
	%   / 
	%  -y         
	%		4--------3
	%	    |		 |		
	%		|		 |
	%		1--------2
	%	rectangular element		
	eleSize_ = min((boundingBox_(2,:)-boundingBox_(1,:))./[nelx_ nely_]);
	carEleMapBack_ = find(1==voxelizedVolume_);	
	numEles_ = length(carEleMapBack_);
	carEleMapForward_ = zeros(nelx_*nely_,1);	
	carEleMapForward_(carEleMapBack_) = (1:numEles_)';
	nodenrs = reshape(1:(nelx_+1)*(nely_+1), 1+nely_, 1+nelx_);
	eNodVec = reshape(nodenrs(1:end-1,1:end-1)+1, nelx_*nely_, 1);
	eNodMat_ = repmat(eNodVec(carEleMapBack_),1,4);
	tmp = [0 nely_+[1 0] -1];
	for ii=1:4
		eNodMat_(:,ii) = eNodMat_(:,ii) + repmat(tmp(ii), numEles_,1);
	end
	carNodMapBack_ = unique(eNodMat_);
	numNodes_ = length(carNodMapBack_);
	carNodMapForward_ = zeros((nelx_+1)*(nely_+1),1);
	carNodMapForward_(carNodMapBack_) = (1:numNodes_)';	
	for ii=1:4
		eNodMat_(:,ii) = carNodMapForward_(eNodMat_(:,ii));
	end
	nodeCoords_ = zeros((nelx_+1)*(nely_+1),2);

	xSeed = boundingBox_(1,1):(boundingBox_(2,1)-boundingBox_(1,1))/nelx_:boundingBox_(2,1);
	ySeed = boundingBox_(2,2):(boundingBox_(1,2)-boundingBox_(2,2))/nely_:boundingBox_(1,2);		
	nodeCoords_(:,1) = reshape(repmat(xSeed, nely_+1, 1), (nelx_+1)*(nely_+1), 1);
	nodeCoords_(:,2) = repmat(ySeed, 1, nelx_+1)';	
	
	nodeCoords_ = nodeCoords_(carNodMapBack_,:);
end
