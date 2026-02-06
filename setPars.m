% -------------------------------------------------------------------------
% setPars.m -- Set default parameters for solving the PDE problem. 
%              setPars; lists all parameters and options. 
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function pars = setPars(varargin)
%
% Listing of parameters and default values
% 
if nargin == 0 && nargout == 0
    disp("Parameters required by RBFsolver, options are given between [] and default parameters given between {}:")
    disp("        dim: Dimension of problem [ positive integer in range [1 | {2} | 3]")
    disp("       geom: Geometry to solve the problem in [ 'cube' of 'ball' {'ball'}]")
    disp("       prob: Problem to solve [ 'poisson' {'poisson'} ]")
    disp("     method: Localised RBF numerical method to use [ 'PUM' | 'FD' | {FD} ]")
    disp("       mode: Extra classification of method used [ {'collocation'} | 'fitted' | 'unfitted']")
    disp("     bcMode: Way to impose boundary conditions only relevant for fitted method [ 'strong' | {'weak'} ]")
    disp("    scaling: Flag to scale the LS problem [ 0 | {1} ]")
    disp("    display: Display flag [ 0 | {1} ]")
    disp("  mvCentres: Flag to move evaluation points on centre points, irrelevant for collocation [ 0 | {1} ]")
    disp("        psi: Compactly supported weight function for PUM [ 'bmp' | {'w2'} ]")
    disp("        phi: Local RBF [ {'phs'} | 'mq' | 'gs' | 'rbfqr' | 'iq']")
    disp("         ep: Shape parameter (smooth phi) [ positive scalar {0.1} ], order for 'phs' [ positive integer {3} ] ")
    disp("       pdeg: Degree of polynomial added to basis, p [ positive integer | 0 | {-1} ]")
    disp("     memTol: Tolerance used to return error when system is too large to solve with cpp library [ positive scalae {0.1} ]")
    disp("        del: Relative overlap between patches [ positive scalar {0.4} ]")
    disp("          q: Oversampling factor between nodes and evaluation points [ positive scalar > 1 {3} ]")
    disp("     rbfdeg: Polynomial degree localy supported by either patches or stencils in case of smooth basis [ positive integer ]")
    disp("          P: Number of generated patches for PU methods [ positive integer {50} ]")
    disp("          N: Global number of centre/node points for FD methods [ positive integer {200} ] ")
    disp("   extCoeff: Fraction of stencil used for the extension in unfitted RBF-FD method case [ positive scalar {0.5} ]")
    return
end
%
%
%
fieldName = ["dim", 'dimPositiveInt';
             "geom", 'positiveDoubleZero';
             "aspectRatio", 'positiveDouble';
             "nLoc", 'positiveInt';
             "numPatches", 'positiveInt';
             "overlap", 'positiveDouble';
             "ep", 'positiveDoubleZero';
             "numPtch0", 'positiveInt';
             "phi", 'phiType';
             "pdeg", 'positiveMinusOneInt';
             "nodeGen", 'nodeGenType';
             "debug", 'flag';
             "display", 'flag';
             "suppressOutput", 'flag';
             "noRef", 'flag',
             "ptchMargin", 'positiveDouble';
             "psi", 'phiType';
             "nrmValWght", 'positiveDouble';
             "seed", 'positiveIntZero';
             "smallWdthTol", "positiveDouble";
             "neiNumFlip", "positiveInt";
             "tolBnd", 'positiveDouble';
             "maxitBnd", 'positiveInt';
             "maxNewtonMvTol", 'positiveDouble'; 
             "neiTol", 'numeric';
             "nClosePtch", 'positiveInt';
             "nPCAmin", 'positiveInt';
             "maxPtchRef", 'positiveDouble';
             "minPtchMerge", 'positiveDouble';
             "scaleData", 'flag';
             "dataScale", 'positiveDouble']; % The numeric test works for intended doubles

%
% For all other call signatures, dim is a required first argument
% 
if (nargin > 0 & ispositiveInt(varargin{1}))
    dim = varargin{1};

    % Set default values according to dim
    if (dim==2)
        pars.extraTol = 1.5;
        pars.dataTol = 0.5;
        pars.aspectRatio = 2;
        pars.nLoc = 21;
        pars.numPatches = 24;
        pars.overlap = 0.5;
        pars.ep = 0.5;
        pars.numPtch0 = 10;
        pars.maxitBnd = 30;
        pars.nClosePtch = 5;
        pars.nPCAmin = 12;
        pars.smallWdthTol = 0.2;
        
    elseif(dim==3)
        pars.extraTol = 2;
        pars.dataTol = 0.1;
        pars.aspectRatio = 0.8;
        pars.nLoc = 165;
        pars.numPatches = 250;
        pars.overlap = 0.7;
        pars.ep = 0.95;
        pars.numPtch0 = 20;
        pars.maxitBnd = 50;        
        pars.nClosePtch = 12;
        pars.nPCAmin = 25;
        pars.smallWdthTol = 0.15;
    else
        error('Only dim=2 and dim=3 are currently supported')
    end

    % Set dimension agnostic parameters
    pars.phi = "mq";
    pars.pdeg = -1;
    pars.nodeGen = "halt";
    pars.debug = 0;
    pars.display = 1;
    pars.suppressOutput = 0;
    pars.noRef = 0;
    pars.ptchMargin = 0.01;
    pars.psi = "bump";
    pars.nrmValWght = 0.05;
    pars.seed = 0;
    pars.neiNumFlip = 1;
    pars.tolBnd = 5e-7;
    pars.maxNewtonMvTol = 0.01;
    pars.neiTol = 0.3;
    pars.maxPtchRef = 0.9;
    pars.minPtchMerge = 0.9;
    pars.scaleData = 0;
    pars.dataScale = 1;
else
    error('The parameter dim is required and must take the values 2 or 3')
end    
  
pin = 0; % Indicates that the second parameter pars is not present

% Check if there is a pars parameter
argnum = 2; % The place where we expect pars
if (mod(nargin-1,2)==1  & isstruct(varargin{argnum}))
    % First argument is a pars struct that we want to update
    for i = 1:size(fieldName,1)
        if isfield(varargin{argnum},fieldName(i,1))

          typeTest = str2func(strcat('is',fieldName(i,2)));
          istype = typeTest(varargin{argnum}.(fieldName(i,1)));

          if ~istype
            error("""%s"" should be of the expected data type.",fieldName(i,1));
          end
	  
          pars.(fieldName(i,1)) = varargin{argnum}.(fieldName(i,1));
        end
     end
     pin = 1; % There was a second argument pars
end
  
  for i = (pin+2):2:nargin

    if ~(isa(varargin{i},'string') | isa(varargin{i},'char'))
        error("Argument %i should be a string equal to the name of a specific field.",i);
    end
    
    if ismember(varargin{i},fieldName) && nargin > i
        loc = find(fieldName(:,1) == varargin{i});
	typeName = fieldName(loc, 2);
	
	isfunc = strcat('is',typeName);
	typeTest = str2func(isfunc);

	istype = typeTest(varargin{i+1});
	if ~istype
          error("""%s"" should be of the expected data type.",varargin{i});
	end

        pars.(fieldName(loc,1)) = varargin{i+1};
    else
        error("Field name ""%s"" not available.",num2str(varargin{i}));
    end

  end  
end

function isf = isflag(data)
  isf = (data==0 | data==1);
end

function ispi = ispositiveInt(data)
  % Note that floor of a string will work without the first check
  ispi = isnumeric(data) & all(data==floor(data) & data>0);
end

function ispi = ispositiveIntZero(data)
  ispi = isnumeric(data) & all(data==floor(data) & data>=0);
end

function ispmi = ispositiveMinusOneInt(data)
  ispmi = isnumeric(data) & all(data==floor(data) & data>=-1);
end  

function ispd = ispositiveDouble(data)
  ispd = all(isnumeric(data) & data>0);
end

function ispd = ispositiveDoubleZero(data)
  ispd = all(isnumeric(data) & data>=0);
end

function ispd = ispositiveDoublePlusOne(data)
  ispd = all(isnumeric(data) & data>=1);
end

function ispt = isphiType(data)
  names = {'mq','bump','wendland_c2','phs','tps','iq','r3','gs','Bmq','rbfqr'};
  ispt = (isa(data,'string') | isa(data,'char'));
  if ispt
    ispt = ispt & ismember(data,names);
  end  
end

function isngt = isnodeGenType(data)
  names = {'halt','uni'};
  isngt = (isa(data,'string') | isa(data,'char'));
  if isngt
    isngt = isngt & ismember(data,names);
  end  
end
