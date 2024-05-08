function outputy = embed(entriesx, basexvector, entriesy)

outputyinit = zeros(1,length(basexvector));

for i = 1:length(outputyinit)
    for j = 1:length(entriesx)
        if basexvector(i) == entriesx(j)
            % store entriesy(j) in outputyinit(i)
            outputyinit(i) = entriesy(j);
        end
    end
end

outputy = outputyinit;

