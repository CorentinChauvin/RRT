function updateInOpenSet( indexInOpenSet , evaluation ) ;
global nodesInOpenSet openSet 

    if indexInOpenSet < 1 ,
        error( ['Illegal index in openSet: ',int2str(indexInOpenSet)] ) ;
    end
    if indexInOpenSet > nodesInOpenSet ,
        error( ['Impossible to update openSet at index ', ...
                int2str(indexInOpenSet),': out of bounds'] );
    end
    openSet( indexInOpenSet , 3 ) = evaluation ;
end

