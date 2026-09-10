% Tandem of two caches sharing one item set.
%
% A request reads Cache1; on a miss it reads Cache2 for THE SAME item; on a
% second miss it is fetched from the origin. Item identity is carried by giving
% each item its OWN job class, so no class switching happens on an arc into a
% cache and the miss of item i at Cache1 IS the read class of item i at Cache2.
%
% Item popularity is expressed by the per-class request rates rather than by a
% popularity distribution the cache draws from, which is also how the TTL
% cache-network analysis is parameterised.
%
% Note on the capacities: with LRU(1) at BOTH levels, Cache2 can never hit. A
% miss at Cache1 inserts the item there, so the next Cache1 miss is necessarily
% for a different item and consecutive arrivals at Cache2 are always distinct.
% Cache2 therefore needs room for more than one item to show any hit at all.
clc; clear solver AvgTable;

n       = 3;                          % items, shared by both caches
pAccess = [0.5 0.3 0.2];              % item popularity
lambda  = 1.0;                        % request rate scale

model = Network('CacheTandem');

think  = Delay(model, 'Think');
cache1 = Cache(model, 'Cache1', n, 1, ReplacementStrategy.LRU);
cache2 = Cache(model, 'Cache2', n, 2, ReplacementStrategy.LRU);

% One dedicated class per item for each role.
read = cell(1,n); hit1 = cell(1,n); hit2 = cell(1,n); miss = cell(1,n);
for f = 1:n
    read{f} = ClosedClass(model, sprintf('Read%d', f), 1, think, 0);
    hit1{f} = ClosedClass(model, sprintf('Hit1_%d', f), 0, think, 0);
    hit2{f} = ClosedClass(model, sprintf('Hit2_%d', f), 0, think, 0);
    miss{f} = ClosedClass(model, sprintf('Miss_%d', f), 0, think, 0);
    think.setService(read{f}, Exp(lambda * pAccess(f)));
end

cache1.setItemReadClasses(read, hit1);          % read{f} reads item f at Cache1
cache1.setMissCache(read{1}, cache2, hit2);     % Cache1 miss of item f = Cache2 read of item f
cache2.setItemMissClass(read{1}, miss);         % Cache2 miss leaves the network

P = model.initRoutingMatrix();
for f = 1:n
    P{read{f}, read{f}}(think, cache1) = 1.0;
    P{hit1{f}, read{f}}(cache1, think) = 1.0;   % hit at Cache1
    P{hit2{f}, read{f}}(cache2, think) = 1.0;   % hit at Cache2
    P{miss{f}, read{f}}(cache2, think) = 1.0;   % miss at both
end
% Cache1 -> Cache2 for each item is registered by setMissCache.
model.link(P);

CTMC(model).getAvgCacheTable
SSA(model, 'samples', 5e4, 'seed', 1).getAvgCacheTable
LDES(model, 'samples', 2e5, 'seed', 1).getAvgCacheTable
