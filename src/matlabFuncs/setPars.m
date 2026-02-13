% -------------------------------------------------------------------------
% setPars.m -- Set default parameters for solving the PDE problem.
%              Outputs structure with relevant parameters.
% Syntax    -- setPars;                                            -> Lists all parameters and options.
%              pars = setPars;                                     -> Sets default parameters. 
%              pars = setPars("geom",'ball',"method",'PUM',"P",10) -> Sets some parameters to non-default values.
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function pars = setPars(varargin)
    %
    % Listing of parameters and default values
    % 
    if nargout == 0
        disp("Parameters required by RBFsolver, options are given between [] and default parameters given between {}:")
        disp("        dim: Dimension of problem [ positive integer in range [ 1 | {2} | 3 ]")
        disp("       geom: Geometry to solve the problem in [ 'cube' of 'ball' {'ball'} ]")
        disp("       prob: Problem to solve [ {'Poisson'} ]")
        disp("     method: Localised RBF numerical method to use [ 'PUM' | {'FD'} ]")
        disp("       mode: Extra classification of method used [ {'collocation'} | 'fitted' | 'unfitted']")
        disp("     bcMode: Way to impose boundary conditions only relevant for fitted method [ 'strong' | {'weak'} ]")
        disp("    scaling: Flag to scale the LS problem [ 0 | {1} ]")
        disp("    display: Display flag [ 0 | {1} ]")
        disp("  mvCentres: Flag to move evaluation points on centre points, irrelevant for collocation [ 0 | {1} ]")
        disp("        psi: Compactly supported weight function for PUM [ 'bmp' | {'w2'} ]")
        disp("        phi: Local RBF [ {'phs'} | 'mq' | 'gs' | 'rbfqr' | 'iq' ]")
        disp("         ep: Shape parameter (smooth phi) [ positive scalar {0.1} ], order for 'phs' [ positive integer {3} ] ")
        disp("       pdeg: Degree of polynomial added to basis, p [ positive integer | 0 | {-1} ]")
        disp("     memTol: Tolerance used to return error when system is too large to solve with cpp library [ positive scalar {0.1} ]")
        disp("        del: Relative overlap between patches [ positive scalar <= 1 {0.4} ]")
        disp("          q: Oversampling factor between nodes and evaluation points [ positive scalar >= 1 {3} ]")
        disp("     rbfdeg: Polynomial degree localy supported by either patches or stencils in case of smooth basis [ positive integer ]")
        disp("          P: Number of generated patches for PU methods [ positive integer {50} ]")
        disp("          N: Global number of centre/node points for FD methods [ positive integer {200} ] ")
        disp("   extCoeff: Fraction of stencil used for the extension in unfitted RBF-FD method case [ positive scalar {0.5} ]")
        disp("      cppOn: Flag to choose whether to use generated c++ library [ 0 | {1} ]")
        disp("      debug: Debug flag specifically for cpp functions [ {0} | 1 ]")
        disp("      parOn: Include parallelisation in implementation [ {0} | 1 ]")
        return
    end
    %
    % Defining types for all inputs
    %
    fieldName = ["dim", 'dimType';
                 "geom", 'geomType';
                 "prob", 'probType';
                 "method", 'methodType';
                 "mode", 'modeType';
                 "bcMode", 'bcModeType';
                 "scaling", 'flag';
                 "display", 'flag';
                 "mvCentres", 'flag';
                 "psi", 'psiType';
                 "phi", 'phiType';
                 "ep", 'epType';
                 "pdeg", 'positiveMinusOneInt';
                 "memTol", 'positiveDouble';
                 "del", 'positiveDoubleLessOne';
                 "q", 'positiveDoublePlusOne';
                 "rbfdeg", 'positiveInt';
                 "P", 'positiveInt';
                 "N", 'positiveInt';
                 "extCoeff", "positiveDouble";
                 "cppOn", "flag";
                 "debug", "flag";
                 "parOn", "flag"]; 
    
    %
    % Set default values
    % 
    pars.dim = 2;                        
    pars.geom = 'ball';   
    pars.prob = 'Poisson';    
    pars.method = 'PUM';          
    pars.mode = 'unfitted';     
    pars.bcMode = 'weak';                     
    pars.scaling = 1;                           
    pars.display = 1;                        
    pars.mvCentres = 1;                 
    pars.psi = 'w2';                    
    pars.phi = 'rbfqr';                             
    pars.pdeg = -1;     
    pars.ep = .1;
    pars.memTol = 0.1;                        
    pars.del = 0.4;       
    pars.rbfdeg = 4;     
    pars.P = 50; 
    pars.N = 200;
    pars.extCoeff = 0.5;
    pars.q = 3;                     
    pars.ep = 0.1;     
    pars.cppOn = 1;
    pars.debug = 1;
    pars.parOn = 0;

    pin = 0; % Indicates that the second parameter pars is not present
    
    for i = (pin+1):2:nargin
        
        if ~ismember(varargin{i},fieldName(:,1))
            error("setPars:IncorrectArg","setPars: Argument %i should be a string equal to the name of an available parameter.",i);
        end
        
        fld = find(fieldName(:,1) == varargin{i});
        typeName = fieldName(fld, 2);
        
        isfunc = strcat('is',typeName);
        typeTest = str2func(isfunc);
        
        if ~strcmp(fieldName(fld, 1),'ep')
            istype = typeTest(varargin{i+1});
        else
            istype = typeTest(varargin{i+1},pars.phi);
        end

        if ~strcmp(fieldName(fld, 1),'ep')
        else
        end
        if ~istype
              error("setPars:IncorrectType","setPars: ""%s"" should be of the expected data type.",varargin{i});
        end
    
        pars.(fieldName(fld,1)) = varargin{i+1};
    end  
end

%
% Functions used to check type of input parameter types
%
function isf = isdimType(data)
  isf = (data==1 | data==2 | data==3);
end

function isf = isgeomType(data)
  typeOpt = {'cube','ball'};
  isf = (isa(data,'string') | isa(data,'char'));
  if isf
    isf = isf & ismember(data,typeOpt);
  end  
end

function isf = isprobType(data)
  typeOpt = {'Poisson'};
  isf = (isa(data,'string') | isa(data,'char'));
  if isf
    isf = isf & ismember(data,typeOpt);
  end  
end

function isf = ismethodType(data)
  typeOpt = {'PUM','FD'};
  isf = (isa(data,'string') | isa(data,'char'));
  if isf
    isf = isf & ismember(data,typeOpt);
  end  
end

function isf = ismodeType(data)
  typeOpt = {'collocation','fitted','unfitted'};
  isf = (isa(data,'string') | isa(data,'char'));
  if isf
    isf = isf & ismember(data,typeOpt);
  end  
end

function isf = isbcModeType(data)
  typeOpt = {'strong','weak'};
  isf = (isa(data,'string') | isa(data,'char'));
  if isf
    isf = isf & ismember(data,typeOpt);
  end  
end

function isf = isflag(data)
  isf = (data==0 | data==1);
end

function isf = ispsiType(data)
  names = {'bmp','w2'};
  isf = (isa(data,'string') | isa(data,'char'));
  if isf
    isf = isf & ismember(data,names);
  end  
end

function isf = isphiType(data)
  names = {'mq','phs','iq','gs','rbfqr'};
  isf = (isa(data,'string') | isa(data,'char'));
  if isf
    isf = isf & ismember(data,names);
  end  
end

function isf = isepType(data,phiPar)
    if strcmp(phiPar,'phs')
        isf  = ispositiveInt(data);
    else
        isf = ispositiveDouble(data);
    end
end

function isf = ispositiveInt(data)
  isf = isnumeric(data) & all(data==floor(data) & data>0);
end

function isf = ispositiveDouble(data)
  isf = all(isnumeric(data) & data>0);
end

function isf = ispositiveMinusOneInt(data)
  isf = isnumeric(data) & all(data==floor(data) & data>=-1);
end  

function isf = ispositiveDoublePlusOne(data)
  isf = all(isnumeric(data) & data>=1);
end

function isf = ispositiveDoubleLessOne(data)
  isf = all(isnumeric(data) & data<=1);
end
