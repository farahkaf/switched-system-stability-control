function ok = isSubsetScaled(S1,S2,gamma)
    S2g = gamma*S2;
    V = S1.V';
    ok = all(S2g.contains(V));
end