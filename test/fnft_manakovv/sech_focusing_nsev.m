A = sym(3.2);
        hlf = sym(1)/sym(2);
        im = sym(1j);
        syms lam;
        a(lam) = gamma(-im*lam + hlf).^2 ./ ( gamma(-im*lam + A + hlf) ...
            * gamma(-im*lam - A + hlf) );
        da = diff(a, lam);
        b(lam) = im*sin(sym(pi)*A) ./ cosh(sym(pi) * lam);

        bound_states = sym(1j)*(A - (floor(A):-1:1) + hlf)
        normconsts = b(bound_states)
        syms h;
        da_vals = limit(da(bound_states + h), h, 0);
        residues = normconsts ./ da_vals

        xi = (-sym(7):sym(8))/5;
        XI = [xi(1) xi(end)]
        digits(40);
        contspec = vpa(b(xi) ./ a(xi)).'
        ab = vpa([a(xi) b(xi)]).'
        
        