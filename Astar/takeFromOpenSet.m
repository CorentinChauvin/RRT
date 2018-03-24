function [node,predecessor,eval] = takeFromOpenSet
global nodesInOpenSet openSet 

    if nodesInOpenSet == 0 ,
        error('Impossible to take node in openSet: empty.');
    end
    % The node with lowest evaluation must be taken from open
    iMin = 1 ;
    evalMin = openSet(1,3) ;
    for i = 2 : nodesInOpenSet ,
        if openSet(i,3) < evalMin ,
            evalMin = openSet(i,3) ;
            iMin = i ;
        end
    end
    node        = openSet(iMin,1) ;
    predecessor = openSet(iMin,2) ;
    eval        = evalMin         ;
    % Remove the node and push nodes to its right one spot to the left.
    for i = iMin : nodesInOpenSet-1 ,
        openSet(i,:) = openSet(i+1,:) ;
    end
    openSet(nodesInOpenSet,:) = [ 0 0 0 ] ; 
    nodesInOpenSet = nodesInOpenSet-1 ;
end

