# ArbOvl
ArbOvl supports computing pairwise and multiple genomic interval sets intersection under various conditions and standards, the main output is a
    list containing two objects: the first one is a binary data frame representing intersection sizes in each combination, the second one is a list
    containing intersecting intervals of each combination. Users can use the binary data frame as input to visualize the result, both venn plot(up to 5
    sets) and upset plot are supported. Users can also output the intersecting intervals into a bed-format data frame.
