exampleName = 'cqn_multiserver';
evalc('cqn_multiserver');
TOL=1e-2; % 1-percent tolerance

spaceRunningsaved = [
     4     3     2     1     0     0     0
     4     3     2     0     1     0     0
     3     4     2     1     0     0     0
     3     4     2     0     1     0     0
     ];
 
% One row, not two: State.fromMarginalAndStarted returns the single canonical
% (descending-sorted) buffer ordering instead of enumerating every permutation.
spaceStartedsaved = [
     4     3     2     1     0     0     0
    ];

spacesaved = [
     4     3     2     1     0     0     0
     4     3     2     0     1     0     0
     4     2     2     0     0     1     0
     4     1     1     1     0     1     0
     4     1     1     0     1     1     0
     3     4     2     1     0     0     0
     3     4     2     0     1     0     0
     3     2     2     0     0     0     1
     3     1     1     1     0     0     1
     3     1     1     0     1     0     1
     2     4     2     0     0     1     0
     2     3     2     0     0     0     1
     2     1     1     0     0     1     1
     1     4     1     1     0     1     0
     1     4     1     0     1     1     0
     1     3     1     1     0     0     1
     1     3     1     0     1     0     1
     1     2     1     0     0     1     1
     1     1     0     1     0     1     1
     1     1     0     0     1     1     1
    ];

try
    assert(max(max(abs(spaceRunning-spaceRunningsaved)))<TOL,sprintf('%s changed on spaceRunning.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(spaceRunning),mat2str(spaceRunningsaved))
end

try
    assert(max(max(abs(spaceStarted-spaceStartedsaved)))<TOL,sprintf('%s changed on spaceStarted.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(spaceStarted),mat2str(spaceStartedsaved))
end

try
    assert(max(max(abs(space-spacesaved)))<TOL,sprintf('%s changed on space.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(space),mat2str(spacesaved))
end

