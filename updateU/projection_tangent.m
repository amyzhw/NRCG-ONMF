  function Up = projection_tangent(X, U)
        
        XtU = multiprod(multitransp(X), U);
        symXtU = multisym(XtU);
        Up = U - multiprod(X, symXtU);
        
    end