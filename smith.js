/* Find g = gcd(a, b) and find x, y with x*a + y*b == g. */
function xgcd(a, b) {
    // We have the following equations:
    // 1a + 0b == a
    // 0a + 1b == b
    // Do a Euclidean algorithm to the right-hand side,
    // Carrying the left-hand side along for the ride.
    let [x, next_x] = [1n, 0n];
    let [y, next_y] = [0n, 1n];
    let [g, next_g] = [BigInt(a), BigInt(b)];
    while (next_g) {
        q = g / next_g;
        [x, next_x] = [next_x, x - q * next_x];
        [y, next_y] = [next_y, y - q * next_y];
        [g, next_g] = [next_g, g - q * next_g];
    }
    if (g < 0n) {
        [x, y, g] = [-x, -y, -g];
    }
    return [x, y, g];
}

function check_rectangular(A, num_cols) {
    A.forEach((row) => {
        if (row.length != num_cols) {
            throw new Error("num_cols does not match given matrix");
        }
    })
}

function identity_matrix(n) {
    let M = [];
    for (i = 0; i < n; i++) {
        let row = [];
        for (j = 0; j < n; j++) {
            row.push((i == j) ? 1n : 0n);
        }
        M.push(row);
    }
    return M;
}

function smithify(A, num_cols) {
    // Find matrices S and T such that SAT is diagonal.
    check_rectangular(A, num_cols);
    let D = A.map((row) => row.map((x) => BigInt(x)));
    let m = D.length;
    let n = num_cols;

    // Initialize S, T, and their inverses to identity matrices.
    let S = identity_matrix(m);
    let Sinv = identity_matrix(m);
    let T = identity_matrix(n);
    let Tinv = identity_matrix(n);

    // The Algorithm
    // -------------
    // Do row and column operations to mutate D
    // until D in Smith Normal Form.
    // Do the corresponding row ops to S and col ops to T
    // to maintain the invariant the SAT == D.
    // Also keep track of the matrix inverses Sinv and Tinv.

    function generalized_row_op(i1, i2, x, y, z, w) {
        // Mutate D by left-multiplying
        // the two-row matrix [D[i1,:] ; D[i2, :]]
        // by [x y ; z w]
        for (jj = 0; jj < n; jj++) {
            let aa = D[i1][jj];
            let bb = D[i2][jj]
            D[i1][jj] = x*aa + y*bb;
            D[i2][jj] = z*aa + w*bb;
        }
        // The corresponding effects on S and Sinv
        for (jj = 0; jj < m; jj++) {
            let aa = S[i1][jj];
            let bb = S[i2][jj];
            S[i1][jj] = x*aa + y*bb;
            S[i2][jj] = z*aa + w*bb;
        }
        for (jj = 0; jj < m; jj++) {
            let aa = Sinv[jj][i1];
            let bb = Sinv[jj][i2];
            Sinv[jj][i1] = w*aa - z*bb;
            Sinv[jj][i2] = x*bb - y*aa;
        }
    }

    function improve_with_row_ops(i1, i2, j) {
        let a = D[i1][j];
        let b = D[i2][j];
        // Pretend that [a; b] is a column matrix.
        // Left-multiply it by some other matrix
        // to turn it into [gcd(a,b); 0].
        // Carry the rest of the respective rows along for the ride.

        // We could do this using only row-operations back and forth,
        // doing the entire Euclidean algorithm to both rows,
        // but with the generalized row op approach,
        // all the heavy lifting can be done ahead of time
        // so we only need one sweep across the rows.
        if (b == 0) {
            return;
        }
        else if (a != 0 && b % a == 0) {
            // This case is important for termination.
            // [ 1  0] [ a]  ==  [a]
            // [-q  1] [qa]      [0]
            q = b / a
            generalized_row_op(i1, i2, 1n, 0n, -q, 1n)
            console.assert(D[i1][j] == a);
            console.assert(D[i2][j] == 0);
        }
        else {
            let [x, y, g] = xgcd(a, b);
            // [x       y] [a]  ==  [g]
            // [-b/g  a/g] [b]      [0]
            generalized_row_op(i1, i2, x, y, -b/g, a/g);
            console.assert(D[i1][j] == g);
            console.assert(D[i2][j] == 0);
        }
    }

    function generalized_col_op(j1, j2, x, y, z, w) {
        // Mutate D by right-multiplying
        // the two-column matrix D[:, j1] | D[:, j2]
        // by [x y ; z w]
        for (ii = 0; ii < m; ii++) {
            let aa = D[ii][j1];
            let bb = D[ii][j2];
            D[ii][j1] = x*aa + y*bb;
            D[ii][j2] = z*aa + w*bb;
        }
        for (ii = 0; ii < n; ii++) {
            let aa = T[ii][j1];
            let bb = T[ii][j2];
            T[ii][j1] = x*aa + y*bb;
            T[ii][j2] = z*aa + w*bb;
        }
        for (ii = 0; ii < n; ii++) {
            let aa = Tinv[j1][ii];
            let bb = Tinv[j2][ii];
            Tinv[j1][ii] = w*aa - z*bb;
            Tinv[j2][ii] = x*bb - y*aa;
        }
    }

    function improve_with_col_ops(j1, j2, i) {
        // Pretend that [a; b] is a column matrix.
        // Right-multiply it by some other matrix
        // to turn it into [gcd(a,b); 0].
        // Carry the rest of the respective columns along for the ride.
        let a = D[i][j1];
        let b = D[i][j2];
        if (b == 0) {
            return;
        }
        else if (a != 0 && b % a == 0) {
            // [a  qa]  [1  -q]  ==  [a  0]
            //          [0   1]
            q = b / a
            generalized_col_op(j1, j2, 1n, 0n, -q, 1n)
            console.assert(D[i][j1] == a);
            console.assert(D[i][j2] == 0n);
        }
        else {
            let [x, y, g] = xgcd(a, b);
            // [a  b]  [x  -b/g]  ==  [g  0]
            //         [y   a/g]
            generalized_col_op(j1, j2, x, y, -b/g, a/g)
            console.assert(D[i][j1] == g);
            console.assert(D[i][j2] == 0n);
        }
    }


    // ===== Phase 1: diagonalize =====
    for (k = 0; k < m && k < n; k++) {
        while (1) {
            for (i = k + 1; i < m; i++) {
                improve_with_row_ops(k, i, k);
            }
            let done_in_row = true;
            for (j = k + 1; j < n; j++) {
                if (D[k][j]) {
                    done_in_row = false;
                }
            }
            if (done_in_row) {
                // The row ops fixing this column didn't mess up this row.
                break;
            }

            for (j = k + 1; j < n; j++) {
                improve_with_col_ops(k, j, k);
            }
            let done_in_col = true;
            for (i = k + 1; i < m; i++) {
                if (D[i][k]) {
                    done_in_col = false;
                }
            }
            if (done_in_col) {
                // The column ops fixing this row didn't mess up this column.
                break;
            }
        }
    }

    // ===== Phase 2: fix the divisibility =====
    while (1) {
        // Bubble the most divisible numbers toward the end.
        let done = true;
        for (k = 0; k < m - 1 && k < n - 1; k++) {
            d1 = D[k][k];
            d2 = D[k+1][k+1];
            if (d1 == 0) {
                if (d2 == 0) {
                    continue;
                }
            }
            else {
                if (d2 % d1 == 0) {
                    continue;
                }
            }
            // Handle a non-divisibility

            // [A  0] [1  0]  ==  [A  0]
            // [0  B] [1  1]  ==  [B  B]
            generalized_col_op(k, k+1, 1n, 1n, 0n, 1n);

            // [A  0]  -->  [gcd(A,B)  X]
            // [B  B]       [       0  Y]
            improve_with_row_ops(k, k+1, k);

            // Because we used row operations, X and Y are multiples of B,
            // which is a multiple of gcd(A,B). To finish off, do
            // [gcd(A,B), X]  -->  [gcd(A,B), 0]
            // [       0, Y]       [       0, Y]
            improve_with_col_ops(k, k+1, k);
            done = false;
        }
        if (done) {
            break;
        }
    }

    return {D:D, S:S, Sinv:Sinv, T:T, Tinv:Tinv};
}

function cokernel(A, num_cols) {
    let m = A.length;
    let n = num_cols;
    check_rectangular(A, n);

    // A = UDV
    smith_A = smithify(A, n);
    U = smith_A.Sinv
    D = smith_A.D

    // coker A = Z^m / im(UDV)
    //         = Z^m / im(UD)
    //         = Ubar * (Z^m / im(D)),
    // where Ubar * [x + (im D)] := [Ux + (im UD)],
    // Note that U is a morphism of pairs:
    //     (X^m, im D) --> (X^m, im UD)
    // The morphism of pairs has an inverse, so when descending
    // to the quotients, Ubar also has an inverse.
    // The generators of Z^m/im(D) are the standard basis, and they
    // generate Z/kZ where k is the corresponding diagonal entry.
    // Apply Ubar to the standard basis.

    // Z/1Z is trivial, Z/0Z is free, and Z/kZ is torsion for k > 1.
    trivialities = []
    torsion_generators = []
    free_generators = []

    for (j = 0; j < m; j++) {
        let order = j < n ? D[j][j] : 0;
        order = order < 0 ? -order : order;
        column = [];
        for (i = 0; i < m; i++) {
            column[i] = U[i][j];
        }
        if (order == 0) {
            free_generators.push(column);
        }
        else if (order == 1) {
            trivialities.push(column);
        }
        else {
            torsion_generators.push([column, order]);
        }
    }

    return {trivialities:trivialities,
            torsion_generators:torsion_generators,
            free_generators:free_generators}
}

function kernel_basis(A, num_cols) {
    let m = A.length;
    let n = num_cols;
    check_rectangular(A, n);
    let smith_A = smithify(A, n);
    let D = smith_A.D
    let T = smith_A.T
    // If D = SAT then {Ax=0} = {SAx=0} = T*{SATx=0} = T*{Dx=0}.
    let ker_D_indices = []
    for (j = 0; j < n; j++) {
        if (j >= m || D[j][j] == 0) {
            ker_D_indices.push(j);
        }
    }
    let columns = [];
    ker_D_indices.forEach((j) => {
        let col = [];
        for (i = 0; i < n; i++) {
            col.push(T[i][j]);
        }
        columns.push(col);
    });
    // TODO: do some operations on these to simplify?
    return columns;
}

function homology(A, num_A_cols, B, num_B_cols) {
    check_rectangular(A, num_A_cols);
    check_rectangular(B, num_B_cols);

    // Compute the homology at
    //        A           B
    // Z^n ------> Z^m ------> Z^k
    let n = num_A_cols;
    let m = A.length;
    let k = B.length;
    if (num_B_cols != m) {
        throw new Error(`Matrices not composable A has ${m} rows, but B has ${num_B_cols} columns.`);
    }

    // Assert BA=0.
    for (i1 = 0; i1 < k; i1++) {
        for (i3 = 0; i3 < n; i3++) {
            s = 0n;
            for (i2 = 0; i2 < m; i2++) {
                s += B[i1][i2] * A[i2][i3];
            }
            if (s != 0n) {
                throw new Error("Matrices do not compose to zero");
            }
        }
    }

    // Since BA=0, we can write Bbar: coker(A) --> Z^k.
    // defined by Bbar([x]) = Bx.
    // Then homology is (ker B)/(im A) = ker Bbar,
    // since these are both the same subquotient.
    let coker_A = cokernel(A, n);

    // compute ker(Bbar: coker(A) --> Z^k)
    // Bbar kills all of the torsion in coker(A)
    // because there's no torsion left in Z^k.
    // The entire torsion part of coker(A)
    // is therefore the torsion part of ker(Bbar).
    let torsion_generators = coker_A.torsion_generators;

    // Write G: Z^(m-r) --> Z^m
    // Such that [] o G: Z^(m-r) --> coker(A)
    // is embedding of free summand
    let G = [];
    let freegen = coker_A.free_generators;
    let m_r = freegen.length;
    // transpose: freegen has m_r lists of m each
    // we want G to have m lists of m_r each
    for (i = 0; i < m; i++) {
        let row = [];
        for (j = 0; j < m_r; j++) {
            row[j] = freegen[j][i];
        }
        G.push(row);
    }

    // H = (ker B) / (im A)
    //   = ker(Bbar: coker(A) --> Z^k)
    //   = ker(Bbar: Tors(coker(A)) (+) Free(coker(A)) --> Z^k)
    //   = Tors(coker(A)) (+) ker(Bbar: [G(Z^(m-r))] --> Z^k)
    // The second summand is
    // ker(Bbar: [G(Z^(m-r))] --> Z^k)
    //   = {[G(v)] in [G(Z^(m-r))] : Bbar([G(v)]) = 0}
    //   = [G({v in Z^(m-r) : Bbar([G(v)]) = 0})]
    //   = [G({v in Z^(m-r) : BG(v) = 0})]
    //   = [G(ker(B o G))]

    let BG = []
    for (i = 0; i < k; i++) {
        BG[i] = [];
        let Bi = B[i];
        for (j = 0; j < m_r; j++) {
            let s = 0n;
            for (q = 0; q < m; q++) {
                s += Bi[q] * G[q][j];
            }
            BG[i][j] = s;
        }
    }
    let ker_BG = kernel_basis(BG, m_r);
    let free_generators = [];
    for (_i = 0; _i < ker_BG.length; _i++) {
        free_generators[_i] = [];
        let gen = ker_BG[_i];
        for (i = 0; i < m; i++) {
            let s = 0n;
            let Gi = G[i];
            for (q = 0; q < Gi.length; q++) {
                s += Gi[q] * gen[q];
            }
            free_generators[_i][i] = s;
        }
    }

    return {torsion_generators: torsion_generators,
            free_generators: free_generators}
}

function chain_complex_from_names(dimension_face_names, boundary) {
    let name_to_dimension = new Map();
    let name_to_index = new Map();
    let dimension_to_size = new Map([[0, 0]]);
    let max_dim = 0;
    let min_dim = 0;
    let matrices = new Map();
    dimension_face_names.forEach((name_list, dim) => {
        dimension_to_size.set(dim, name_list.length);
        if (dim > max_dim) {
            max_dim = dim;
        }
        if (dim < min_dim) {
            min_dim = dim;
        }
    });
    dimension_to_size.set(min_dim - 1, 0);
    dimension_to_size.set(max_dim + 1, 0);

    // Make sure everything lines up
    dimension_face_names.forEach((name_list, dim) => {
        for (index = 0; index < name_list.length; index++) {
            let name = name_list[index];
            if (name_to_dimension.has(name)) {
                throw new Error(`duplicate name "${name}"`);
            }
            name_to_dimension.set(name, dim);
            name_to_index.set(name, index);
        }
    });
    for (dim = min_dim; dim <= max_dim; dim++) {
        if (!dimension_to_size.has(dim)) {
            dimension_to_size.set(dim, 0);
        }
    }
    boundary.forEach((dF, F) => {
        if (!name_to_dimension.has(F)) {
            throw new Error(`Unknown face "${F}"`);
        }
        let dim_F = name_to_dimension.get(F);
        dF.forEach(([coeff, face]) => {
            if (!name_to_dimension.has(face)) {
                throw new Error(`Unknown face "${face}"`);
            }
            dim_face = name_to_dimension.get(face);
            if (name_to_dimension.get(face) != dim_F - 1) {
                throw new Error(`Boundary of ${dim_F}-dimensional ${F} includes ${dim_face}-dimensional ${face}`);
            }
        });
    });

    // Make the matrices
    for (dim = max_dim + 1; dim >= min_dim; dim--) {
        let m = dimension_to_size.get(dim - 1);
        let n = dimension_to_size.get(dim);
        let M = [];
        for (i = 0; i < m; i++) {
            M[i] = [];
            for (j = 0; j < n; j++) {
                M[i][j] = 0n;
            }
        }
        matrices.set(dim, M);
    }
    boundary.forEach((dF, F) => {
        let A = matrices.get(name_to_dimension.get(F));
        dF.forEach(([coeff, face]) => {
            let [i, j] = [name_to_index.get(face), name_to_index.get(F)];
            A[i][j] += BigInt(coeff);
        });
    });

    // Assert chain complex
    boundary.forEach((dF, F) => {
        let dim = name_to_dimension.get(F);
        ddF = new Map();
        if (!dimension_face_names.has(dim-2)) {
            return;
        }
        dimension_face_names.get(dim-2).forEach((name) => {
            ddF.set(name, 0);
        });
        dF.forEach(([coeff, face]) => {
            if (boundary.has(face)) {
                boundary.get(face).forEach(([coeff2, face2]) => {
                    ddF.set(face2, ddF.get(face2) + coeff * coeff2);
                });
            }
        });
        ddF.forEach((coeff, name) => {
            if (coeff != 0) {
                throw new Error(`Boundary of boundary of ${F} was nonzero; included ${coeff} copies of ${name}`);
            }
        });
    });

    return {min_dim:min_dim, max_dim:max_dim, matrices:matrices,
            dimension_to_size:dimension_to_size}
}

function homology_from_names(dimension_face_names, boundary) {
    let {min_dim, max_dim, matrices, dimension_to_size}
        = chain_complex_from_names(dimension_face_names, boundary);
    result = [];
    for (dim = min_dim; dim <= max_dim; dim++) {
        let A = matrices.get(dim + 1);
        let B = matrices.get(dim);
        let n = dimension_to_size.get(dim + 1);
        let m = dimension_to_size.get(dim);
        let H = homology(A, n, B, m);
        let name_list = dimension_face_names.get(dim);
        function to_names(gen) {
            let namegen = [];
            for (i = 0; i < gen.length; i++) {
                let coeff = gen[i];
                if (coeff != 0n) {
                    namegen.push([coeff, name_list[i]]);
                }
            }
            return namegen;
        }
        let free_generators = H.free_generators.map(to_names);
        let torsion_generators = H.torsion_generators.map(
            ([gen, order]) => [to_names(gen), order]
        );
        result.push([dim, {free_generators: free_generators,
                           torsion_generators:torsion_generators}]);
    }
    return result;
}

function coeffs_to_string(gen) {
    let parts = [];
    gen.forEach(([coeff, name]) => {
        parts.push(coeff < 0n ? "-" : "+");
        if (parts.length == 1 && parts[0] == "+") {
            parts.pop();
        }
        let c = coeff < 0n ? -coeff : coeff;
        if (c > 1) {
            parts.push(c.toString());
        }
        parts.push(name);
    });
    return parts.join('');
}

