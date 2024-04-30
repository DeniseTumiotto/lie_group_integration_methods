function new_sols = delete_duplicate(sols)

saved = 0;
j = 1;

for i = 1:max(size(sols))
    if mod(sols{i}.problemset, 10) == 1 && all(saved(:)~=sols{i}.problemset)
        new_sols{j} = sols{i};
        j = j+1;
        saved = [saved sols{i}.problemset];
    elseif mod(sols{i}.problemset, 10) ~= 1
        new_sols{j} = sols{i};
        j = j+1;
    end
end

end