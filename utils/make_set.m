function P = make_set(H,w)
    P = Polyhedron(H,w);
    P = P.minHRep();   
end
