function modelout = changeMinimalMedia(model, limEx, lowerBound)
%changeMedia -  Changes the lower and upper bounds of the model's exchange
%reactions to alter the in silico minimal media with choice of carbon source.

%INPUTS
% model - COBRA model structure
% limEX - cell array containing exchange reactions to turn on in minimal
% media
%
%OUTPUT
% modelout - COBRA model structure with modified exchange reaction bounds
%
% Adapted from the following:
% Matthew Oberhardt, 1-7-2010
% Jennifer Bartell, 3-27-2013
% Anna Blazier, 9-18-2012 
% Laura Dunphy 9-22-2016

% Finds the indices of the exchange reactions in the model
if nargin < 3
  lowerBound = -10;
end 
 
exchangerxns = [];
for rxn = 1:length(model.rxns)
    exchangerxns(rxn) = strncmp('EX_',model.rxns{rxn},3);
end
exchangeindices = find(exchangerxns); 

exchangeindices = find(findExcRxns(model));


modelout = model;

modelout.lb(exchangeindices) = zeros(size(exchangeindices));
modelout.ub(exchangeindices) = 1000*ones(size(exchangeindices));
%Limit aerobic growth to 20 mmol/gDW/hr of O2 (EX_O2 removed from openexchanges)
modelout.lb(findRxnIDs(modelout,{'EX_cpd00007(e)'})) = -20; 


%% Minimal Media
%Nutrients such as ions that are freely exchanged:    
openexchanges = {     
 'EX_cpd00001(e)' %H2O
 'EX_cpd00009(e)' %Phosphate
 'EX_cpd00011(e)' %CO2
 'EX_cpd00021(e)' %Fe2+ LJD 5/26/17
 'EX_cpd00030(e)' %Mn2+
 'EX_cpd00034(e)' %Zn2+
 'EX_cpd00048(e)' %Sulfate
 'EX_cpd00058(e)' %Cu2+
 'EX_cpd00067(e)' %H+
 'EX_cpd00149(e)' %Co2+
 'EX_cpd00205(e)' %K+
 'EX_cpd00254(e)' %Mg
 'EX_cpd00528(e)' % Nitrogen LJD 5/26/17
 'EX_cpd00971(e)' %Na+
 'EX_cpd00013(e)' %NH3
 'EX_cpd01012(e)' %Cd2+
 'EX_cpd10516(e)' %fe3
 'EX_cpd00244(e)' %Ni2+ ADDED BY PHIL
};

limitedexchanges = limEx;

%changes the lower bound of openexchanges to -1000 and the lower bound of
%the limitedexchanges to -10.  Also changes the upper bound of
%limitedexchanges to 0.
modelout.lb(find(ismember(model.rxns,openexchanges))) = -1000*ones(size(openexchanges));
modelout.lb(find(ismember(model.rxns,limitedexchanges))) = lowerBound*ones(size(limitedexchanges));


modelout;

end
