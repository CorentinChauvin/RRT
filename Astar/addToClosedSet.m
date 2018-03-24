function addToClosedSet( node , predecessor , evaluation ) 
global nodesInClosedSet closedSet numberOfNodes

    % Insert node in closedSet. Here closedSet is not ordered
    if nodesInClosedSet == numberOfNodes ,
        error('Impossible to add a node in closedSet: full.');
    end
    nodesInClosedSet = nodesInClosedSet + 1 ;
    closedSet(nodesInClosedSet,:) = [ node predecessor evaluation ] ;
end

