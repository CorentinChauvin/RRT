function initializeOpenSet 
% An information in open is:
%    - node number
%    - node predecessor
%    - node evaluation
     global numberOfNodes openSet nodesInOpenSet ;
     openSet = zeros( numberOfNodes , 3 ) ;
     nodesInOpenSet = 0 ;
end