function removeFromClosedSet(indexInClosedSet)
global nodesInClosedSet closedSet 

    if (indexInClosedSet>nodesInClosedSet) || (indexInClosedSet<1) ,
        error( ['Illegal index in closedSet: ',int2str(indexInClosedSet)] ) ;
    end
    % Push nodes one spot to the left from indexInClosedSet+1
    for i = indexInClosedSet : nodesInClosedSet-1 ,
        closedSet(i,:) = closedSet(i+1,:) ;
    end
    closedSet(nodesInClosedSet,:) = [ 0 0 0 ] ;
    nodesInClosedSet = nodesInClosedSet-1 ;
end

