function saveh5struct(fname,s,groupName)
% save heirarchical structure data an h5 file or subgroup
% by David Frey
%
% inputs:
% fname - name of h5 file to load in
% groupName - name of subgroup to load
% s - structure with same heirarchy as h5 file
%

    import util.*

    if nargin < 3 || isempty(groupName)
        groupName = '/';
    end

    % save datasets in this group
    fn = fieldnames(s);
    for i = 1:length(fn)
        name = fn{i};
        if isstruct(s.(name)) % recursively save subgroups
            saveh5struct(fname,s.(name),[groupName '/' name]);
        else % save the value
            h5create(fname, [groupName '/' name], size(s.(name)), 'Datatype', class(val));
            h5write(fname, [groupName '/' name], s.(name));
        end
    end

end

