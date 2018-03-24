function addToOpenSet( node , predecessor , evaluation ) ;
global nodesInOpenSet numberOfNodes openSet
    % Insert node in openSet. Here openSet is not ordered
    if nodesInOpenSet == numberOfNodes ,
        error('Impossible to add a node in openSet: full.');
    end
    nodesInOpenSet = nodesInOpenSet + 1 ;
    openSet(nodesInOpenSet,:) = [ node predecessor evaluation ] ;
end
