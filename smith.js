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
        let q = g / next_g;
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
        row.forEach((x) => {
            if (typeof x != "bigint") {
                throw new Error("not matrix every entry was a bigint.");
            }
        })
    })
}

function identity_matrix(n) {
    let M = [];
    for (let i = 0; i < n; i++) {
        let row = [];
        for (let j = 0; j < n; j++) {
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
    // until D is in Smith Normal Form.
    // Do the corresponding row ops to S and col ops to T
    // to maintain the invariant that SAT == D.
    // Also keep track of the matrix inverses Sinv and Tinv.

    function generalized_row_op(i1, i2, x, y, z, w) {
        // Mutate D by left-multiplying
        // the two-row matrix [D[i1,:] ; D[i2, :]]
        // by [x y ; z w]
        for (let jj = 0; jj < n; jj++) {
            let aa = D[i1][jj];
            let bb = D[i2][jj]
            D[i1][jj] = x*aa + y*bb;
            D[i2][jj] = z*aa + w*bb;
        }
        // The corresponding effects on S and Sinv
        for (let jj = 0; jj < m; jj++) {
            let aa = S[i1][jj];
            let bb = S[i2][jj];
            S[i1][jj] = x*aa + y*bb;
            S[i2][jj] = z*aa + w*bb;
        }
        for (let jj = 0; jj < m; jj++) {
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
            let q = b / a;
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
        for (let ii = 0; ii < m; ii++) {
            let aa = D[ii][j1];
            let bb = D[ii][j2];
            D[ii][j1] = x*aa + y*bb;
            D[ii][j2] = z*aa + w*bb;
        }
        for (let ii = 0; ii < n; ii++) {
            let aa = T[ii][j1];
            let bb = T[ii][j2];
            T[ii][j1] = x*aa + y*bb;
            T[ii][j2] = z*aa + w*bb;
        }
        for (let ii = 0; ii < n; ii++) {
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
            let q = b / a
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


    // ===== Phase 1: Make D diagonal =====
    for (let k = 0; k < m && k < n; k++) {
        for (;;) {
            for (let i = k + 1; i < m; i++) {
                improve_with_row_ops(k, i, k);
            }
            let done_in_row = true;
            for (let j = k + 1; j < n; j++) {
                if (D[k][j]) {
                    done_in_row = false;
                }
            }
            if (done_in_row) {
                // The row ops fixing this column didn't mess up this row.
                break;
            }

            for (let j = k + 1; j < n; j++) {
                improve_with_col_ops(k, j, k);
            }
            let done_in_col = true;
            for (let i = k + 1; i < m; i++) {
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

    // ===== Phase 2: Fix the divisibility =====
    for (;;) {
        // Bubble the most divisible numbers toward the end.
        let done = true;
        for (let k = 0; k < m - 1 && k < n - 1; k++) {
            let d1 = D[k][k];
            let d2 = D[k+1][k+1];
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
    let smith_A = smithify(A, n);
    let U = smith_A.Sinv;
    let D = smith_A.D;

    // coker A = Z^m / im(UDV)
    //         = Z^m / im(UD)
    //         = Ubar * (Z^m / im(D)),
    // where Ubar * [x + (im D)] := [Ux + (im UD)],
    // Note that U is a morphism of pairs:
    //     (Z^m, im D) --> (Z^m, im UD)
    // The morphism of pairs has an inverse, so when descending
    // to the quotients, Ubar also has an inverse.
    // The generators of Z^m/im(D) are the standard basis, and they
    // generate Z/kZ where k is the corresponding diagonal entry.
    // Apply Ubar to the standard basis.

    // Z/1Z is trivial, Z/0Z is free, and Z/kZ is torsion for k > 1.
    let trivialities = [];
    let torsion_generators = [];
    let free_generators = [];

    for (let j = 0; j < m; j++) {
        let order = j < n ? D[j][j] : 0;
        order = order < 0 ? -order : order;
        let column = [];
        for (let i = 0; i < m; i++) {
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

function matrix_multiply(A, num_A_cols, B, num_B_cols) {
    check_rectangular(A, num_A_cols);
    check_rectangular(B, num_B_cols);
    if (num_A_cols != B.length) {
        throw new Error(`Matrices not composable: left has ${num_A_cols} but right has ${B.length} columns.`);
    }
    let result = [];
    A.forEach((row) => {result.push([])});
    for (let i = 0; i < A.length; i++) {
        for (let j = 0; j < num_B_cols; j++) {
            let s = 0n;
            for (let k=0; k < num_A_cols; k++) {
                s += A[i][k]*B[k][j];
            }
            result[i][j] = s;
        }
    }
    return result;
}

function homology(A, num_A_cols, B, num_B_cols, coeff) {
    check_rectangular(A, num_A_cols);
    check_rectangular(B, num_B_cols);
    if (typeof coeff != "bigint") {
        throw new Error(`modulus ${coeff} is not a BigInt.`);
    }
    matrix_multiply(B, num_B_cols, A, num_A_cols).forEach(
        (row) => row.forEach((x) => {
            if (x != 0n) {
                throw new Error("Matrices do not compose to zero");
            }
        })
    );

    /*
     * Compute the homology at
     *        A'          B'
     * R^n ------> R^m ------> R^k
     * where R is the cyclic group Z/coeff.
     */

    /*
     * We'll operate with the integer matrices A and B
     * instead of the R-matrices A' and B'.
     *            A           B
     *     Z^n ------> Z^m ------> Z^k
     *      |p         |p           |p
     *      v    A'    v     B'     v
     *     R^n ------> R^m ------> R^k
     * 
     * Take Smith normal form SBT = D.
     * Now:
     * 
     *     p(v) in ker B'
     *     iff B'(p(v)) = 0
     *     iff p(B(v)) = 0
     *     iff B(v) in coeff*Z^k
     *     iff Sinv D Tinv v in coeff*Z^k
     *     iff D Tinv v in coeff*Z^k
     *     iff Tv has  Dv in coeff*Z^k
     * 
     * The set {v | D v in coeff*Z^k}
     * has a basis of vectors [0 ... 0 coeff/Dii 0 ... 0].
     * Write E for the matrix with these columns (omitting zero columns),
     * so that {v | Dv in coeff*Z^k} = E Z^l
     * where E is an m-by-l matrix.
     * 
     * Now p(v) in ker(B') iff v in (T E Z^l), so we know the cycles.
     * Now the homology group is
     * 
     *     H = ker(B') / im(A')
     *       = p(T E Z^l) / p(A Z^n)
     * 
     * We can quotient out by the coefficient modulus and by the
     * boundaries simultaneously:
     * 
     *     H = (T E Z^l)/(A Z^n + coeff*Z^m)
     *       = (T E Z^l)/([A|coeff*id] Z^(n+m))
     *       = (T E Z^l)/(T Tinv [A|coeff*id] Z^(n+m))
     *       = (T E Z^l)/(T [Tinv A|coeff*id] Z^(n+m))
     *       = Tbar ( (E Z^l)/(Tinv [A|coeff*id] Z^(n+m)) )
     * 
     * where Tbar([v]):=[Tv] is the map T descended to the quotient.
     * Because E is diagonal and injective,
     * we get can undo multiplication by E
     * to find a matrix F satisfying: E F = Tinv [A|coeff*id]
     * We then have:
     * 
     *     H = Tbar((E Z^l)/(E F Z^(n+m)))
     *       = Tbar Ebar (Z^l / (F Z^(n+m)))
     *       = Tbar(Ebar(coker(F)))
     */

    let smith_B = smithify(B, num_B_cols);
    let T = smith_B.T;
    let Tinv = smith_B.Tinv;
    let D = smith_B.D;

    function cyclic_kernel(x0, m) {
        // Kernel of the multiplication by x0 map Z --x0--> Z/mZ
        // Examples: ker(Z --1--> Z) = 0Z
        //           ker(Z --2--> Z) = 0Z
        //           ker(Z --0--> Z) = 1Z
        //           ker(Z --1--> Z/4) = 4Z
        //           ker(Z --2--> Z/4) = 2Z
        //           ker(Z --0--> Z/4) = 1Z
        let [_, __, x] = xgcd(x0, m);
        if (x == 0n) {
            return 1n;
        }
        else {
            if (m % x !== 0n) {
                throw new Error("bad divisibility");
            }
            return m / x;
        }
    }

    // Consider the m-by-m diagonal matrix E0
    // where the nonzero columns are the basis for ker(D).
    // Because D has the least-divisible entries first,
    // E0 will have the most-divisible entries first,
    // so E0 will have all of the zeros first.
    // The conceptual matrix E will consist of only the nonzero
    // columns of E0.

    // Store the diagonal entries of E0, but note that
    // the entries of E omit the initial string of zeros.
    let E0_entries = [];
    for (let j = 0; j < num_B_cols; j++) {
        let D_entry = j >= B.length ? 0n : D[j][j];
        E0_entries.push(cyclic_kernel(D_entry, coeff));
    }
    let num_zeros;
    for (num_zeros = 0; num_zeros < E0_entries.length; num_zeros++) {
        if (E0_entries[num_zeros] != 0n) {
            break;
        }
    }
    for (let i = num_zeros; i < E0_entries.length; i++) {
        if (E0_entries[i] == 0n) {
            throw new Error("Bad SNF");
        }
    }

    if (num_zeros > 0 && coeff != 0n) {
        throw new Error("cyclic_kernel returned 0 unexpectedly");
    }

    let num_Tinv_A_cols = num_A_cols;
    let Tinv_A = matrix_multiply(Tinv, num_B_cols, A, num_A_cols);
    check_rectangular(Tinv_A, num_Tinv_A_cols);
    console.assert(Tinv_A.length == num_B_cols);
    if (coeff != 0n) {
        // if coefficient modulus is nonzero, kill it
        // by splicing in coeff*(identity matrix)
        for (let i = 0; i < num_B_cols; i++) {
            let row = Tinv_A[i];
            for (let j = 0; j < num_B_cols; j++) {
                row.push(i == j ? coeff : 0n);
            }
        }
        num_Tinv_A_cols = num_A_cols + num_B_cols;
        check_rectangular(Tinv_A, num_Tinv_A_cols);
    }
    for (let i = 0; i < num_zeros; i++) {
        Tinv_A[i].forEach((x) => {
            if (x !== 0n) {
                throw new Error("rows unexpectedly nonzero");
            }
        });
    }
    let F = [];
    for (let i = num_zeros; i < E0_entries.length; i++) {
        let row = Tinv_A[i];
        let E_entry = E0_entries[i];
        row.forEach((x) => {
            if (x % E_entry != 0n) {
                throw new Error("Entry was not divisible");
            }
        });
        F.push(row.map((x) => x / E_entry));
    }
    let cok = cokernel(F, num_Tinv_A_cols);

    function TE(v) {
        console.assert(v.length == E0_entries.length - num_zeros);
        // re-introduce the initial zeros
        let Ev = Array(num_zeros).fill(0n);
        for (let i = num_zeros; i < E0_entries.length; i++) {
            Ev.push(v[i - num_zeros] * E0_entries[i]);
        }
        return T.map((row) => {
            s = 0n;
            for (let i = 0; i < num_B_cols; i++) {
                s += row[i] * Ev[i];
            }
            return s
        });
    }

    let result_torsion_generators = cok.torsion_generators.map(([v, order]) => [TE(v), order]);
    let result_free_generators = cok.free_generators.map(TE)
    return {torsion_generators: result_torsion_generators,
            free_generators: result_free_generators}
}

function chain_complex_from_names(dimension_face_names, boundary, relative) {
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
        for (let index = 0; index < name_list.length; index++) {
            let name = name_list[index];
            if (name_to_dimension.has(name)) {
                throw new Error(`duplicate name "${name}"`);
            }
            name_to_dimension.set(name, dim);
            name_to_index.set(name, index);
        }
    });
    for (let dim = min_dim; dim <= max_dim; dim++) {
        if (!dimension_to_size.has(dim)) {
            dimension_to_size.set(dim, 0);
        }
    }
    boundary.forEach((dF, F) => {
        if (!name_to_dimension.has(F)) {
            throw new Error(`Unknown face "${F}"`);
        }
        let dim_F = name_to_dimension.get(F);
        dF.forEach(([, face]) => {
            if (!name_to_dimension.has(face)) {
                throw new Error(`Unknown face "${face}"`);
            }
            let dim_face = name_to_dimension.get(face);
            if (name_to_dimension.get(face) != dim_F - 1) {
                throw new Error(`Boundary of ${dim_F}-dimensional ${F} includes ${dim_face}-dimensional ${face}`);
            }
        });
    });

    // Make the matrices
    for (let dim = max_dim + 1; dim >= min_dim; dim--) {
        let m = dimension_to_size.get(dim - 1);
        let n = dimension_to_size.get(dim);
        let M = [];
        for (let i = 0; i < m; i++) {
            M[i] = [];
            for (let j = 0; j < n; j++) {
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
        let ddF = new Map();
        if (!dimension_face_names.has(dim-2)) {
            return;
        }
        dimension_face_names.get(dim-2).forEach((name) => {
            ddF.set(name, 0n);
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

    if (relative.length == 0) {
        return {dimension_face_names:dimension_face_names,
                min_dim:min_dim,
                max_dim:max_dim,
                matrices:matrices,
                dimension_to_size:dimension_to_size}
    }

    //////////////////////////////////////////////
    // Now handling relativizing the chain complex
    //////////////////////////////////////////////

    let set_relative = new Set(relative);
    relative.forEach((name) => {
        if (!name_to_dimension.has(name)) {
            throw new Error(`Unknown cell "${name}" to relativize`);
        }
        if (boundary.has(name)) {
            boundary.get(name).forEach(([coeff, face]) => {
                if (!set_relative.has(face)) {
                    throw new Error(`Relative complex not a complex: boundary of ${name} included ${face}`);
                }
            });
        }
    });
    let relative_dimension_face_names = new Map();
    let relative_dimension_to_size = new Map(dimension_to_size);
    dimension_face_names.forEach((name_list, dim) => {
        let filtered = name_list.filter((name)=>!set_relative.has(name));
        relative_dimension_face_names.set(dim, filtered);
        relative_dimension_to_size.set(dim, filtered.length);
    });
    let to_original_index = new Map();
    relative_dimension_face_names.forEach((name_list, dim) => {
        to_original_index.set(dim, name_list.map((name) => name_to_index.get(name)));
    })
    let relative_matrices = new Map();
    // Make the sub-matrices
    for (let dim = max_dim + 1; dim >= min_dim; dim--) {
        let m = relative_dimension_to_size.get(dim - 1);
        let n = relative_dimension_to_size.get(dim);
        let i_to_i0 = to_original_index.get(dim - 1);
        let j_to_j0 = to_original_index.get(dim);
        let M0 = matrices.get(dim);
        let M = [];
        for (let i = 0; i < m; i++) {
            let i0 = i_to_i0[i];
            M[i] = [];
            for (let j = 0; j < n; j++) {
                let j0 = j_to_j0[j];
                M[i][j] = M0[i0][j0];
            }
        }
        relative_matrices.set(dim, M);
    }
    return {dimension_face_names:relative_dimension_face_names,
            min_dim:min_dim,
            max_dim:max_dim,
            matrices:relative_matrices,
            dimension_to_size:relative_dimension_to_size}
}

function transpose(A, num_A_cols) {
    check_rectangular(A, num_A_cols);
    let result = [];
    for (let j = 0; j < num_A_cols; j++) {
        let new_row = [];
        A.forEach((old_row) => {
            new_row.push(old_row[j]);
        });
        result.push(new_row);
    }
    check_rectangular(result, A.length);
    return result;
}

function homology_from_names(dimension_face_names0, boundary, co, coeff, relative) {
    function to_names(name_list, gen) {
        let namegen = [];
        for (let i = 0; i < gen.length; i++) {
            let coeff = gen[i];
            if (coeff != 0n) {
                console.assert(i < name_list.length);
                namegen.push([coeff, name_list[i]]);
            }
        }
        return namegen;
    }
    let {dimension_face_names, min_dim, max_dim, matrices, dimension_to_size}
        = chain_complex_from_names(dimension_face_names0, boundary, relative);
    let result = [];
    for (let dim = min_dim; dim <= max_dim; dim++) {
        let A = matrices.get(dim + 1);
        let B = matrices.get(dim);
        let n = dimension_to_size.get(dim + 1);
        let m = dimension_to_size.get(dim);
        let H;
        if (co) {
            H = homology(transpose(B, m), B.length, transpose(A, n), A.length, coeff);
        }
        else {
            H = homology(A, n, B, m, coeff);
        }
        let name_list = dimension_face_names.get(dim);
        let free_generators = H.free_generators.map(
            (gen) => to_names(name_list, gen)
        );
        let torsion_generators = H.torsion_generators.map(
            ([gen, order]) => [to_names(name_list, gen), order]
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

