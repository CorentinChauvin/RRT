function [present,index] = isInOpenSet( node )
    global nodesInOpenSet openSet 
    present = false ;
    index = 0 ;
    i = 1 ;
    while (i<=nodesInOpenSet) && (openSet(i,1)~=node) ,
        i = i+1 ;
    end
    if i<=nodesInOpenSet ,
        present = true ;
        index = i ;
    end
end

