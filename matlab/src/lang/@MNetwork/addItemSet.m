function addItemSet(self, itemSet)
% ADDITEMSET(ITEMSET)

if sum(cellfun(@(x) strcmp(x.name, itemSet.name), self.items))>0
    line_error(mfilename,sprintf('An item type with name %s already exists.\n', itemSet.name));
end
nItemSet = size(self.items,1);
itemSet.index = nItemSet+1;
self.items{end+1,1} = itemSet;
self.setUsedLangFeature(class(itemSet));
end
