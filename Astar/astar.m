% Find a path in a graph using the A* algorithm
%
% Inputs:
%   - start,goal: start and goal vertices of the graph (integer number).
%   - graph: the sparse NxN matrix of the graph, N number of vertices.
%   - nodeInfo: a N-cell array with the characteristics of each node.
%               Node info depends on the problem at hand.
%   - hStarHandle: the handle to the h* function. Typically, the function
%                  uses goal and nodeInfo to calculate the estimated cost
%                  from a given node to the goal. The function will 
%                  be called by h*(goal,nodeInfo) so h* should know what
%                  nodeInfo is all about.
% Outputs:
%   - path: a list of node numbers, from start to goal.
%   - cost: the cost of the path
% Note:
%   The program displays a progress bar. The progress restarts every 
%   1000 loops: there is no way to know beforehand how many loops will be
%   performed...
%
% Author: 
%   Gaetan Garcia - Ecole Centrale de Nantes

function [path,cost] = astar( start , goal , graph , nodeInfo , hStarHandle )
global numberOfNodes nodesInOpenSet nodesInClosedSet openSet closedSet

    numberOfNodes = size(graph,1) ;
    initializeOpenSet    ;
    initializeClosedSet  ;
    
    goalReached = false ;
    % Put start in open. No predecessor (denoted as 0).
    addToOpenSet( start , 0 , hStarHandle(nodeInfo(start,:)-nodeInfo(goal,:)) ) ;
    
    wbHandle = waitbar(0,'Searching the graph...') ;
    nbLoops = 0 ;
    while ~goalReached && (nodesInOpenSet>0) ,

        nbLoops = nbLoops + 1 ;
        waitbar( mod(nbLoops,1000)/1000 ) ;
        [node,predecessor,eval] = takeFromOpenSet  ;

        [isInClosed,index] = isInClosedSet(node) ;
        if isInClosed ,
            closedSet(index,3) = min( closedSet(index,3) , eval ) ;
        else
            addToClosedSet( node , predecessor , eval ) ;
        end
        if node == goal ,
            goalReached = true ;
        else
            edges = find( graph(node,:) > 0 ) ;
            for i = 1 : length(edges) ,
                successor = edges(i) ;
                successorEval = eval - ...
                    hStarHandle(nodeInfo(node,:)-nodeInfo(goal,:)) + ...
                    hStarHandle(nodeInfo(successor,:)-nodeInfo(goal,:)) + ...
                    graph(node,successor) ;
                [nodeIsInOpen  ,indexInOpen]   = isInOpenSet(successor)   ;
                [nodeIsInClosed,indexInClosed] = isInClosedSet(successor) ;
                if nodeIsInOpen ,
                    if successorEval < openSet(indexInOpen,3) ,
                        updateInOpenSet(indexInOpen,successorEval) ;
                    end
                elseif nodeIsInClosed ,
                    if successorEval < closedSet(indexInClosed,3) ,
                        addToOpenSet(successor,node,successorEval) ;
                        removeFromClosedSet(indexInClosed) ;
                    end
                else
                    addToOpenSet( successor , node , successorEval ) ;
                end
            end
        end
    end
    
    if ~goalReached ,
        path = [] ;
    else
        node = goal ;
        indexInPath = 1 ;
        path(indexInPath) = goal ;
        while node ~= start ,
            [isPresent,index] = isInClosedSet( node ) ;
            if ~isPresent ,
                error(['Can''t find ',num2str(node),' in closedSet!']);
            else
                node = closedSet(index,2) ;
                indexInPath = indexInPath+1 ;
                path(indexInPath) = node ;
            end
        end
    end
    
    cost = 0 ;
    for i = 1 : length(path)-1 ,
        cost = cost + graph(path(i),path(i+1)) ;
    end
    close(wbHandle) ;
end