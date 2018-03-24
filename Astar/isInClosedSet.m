function [present,index] = isInClosedSet( node )
    global nodesInClosedSet closedSet 
    present = false ;
    index = 0 ;
    i = 1 ;
    while (i<=nodesInClosedSet) && (closedSet(i,1)~=node) ,
        i = i+1 ;
    end
    if i<=nodesInClosedSet ,
        present = true ;
        index = i ;
    end
end

