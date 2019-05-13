function pdist = pointMassMixture(points, masses, N)
% This function produces an 1-d array of size N with entries chosen from
% mixture of point masses
% Input: 
% points: points at which the pmf has non-zero mass
% masses: mass/weight at the points
% N length of the output array

pdist = -ones(1, N);

cdfmasses = cumsum(masses);


    cdfmasses = [0, cdfmasses];
    cointosses = rand(1, N);
    lenpoints = length(points);
    for i = 1 : lenpoints
        ind = find( cointosses > cdfmasses(lenpoints - i + 1) );
        pdist(ind) = points(lenpoints - i + 1);
        cointosses(ind) = -1;
        clear ind
    end

