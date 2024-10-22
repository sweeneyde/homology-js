const face_name_regex = /^[a-zA-Z_][a-zA-Z0-9_]*$/;

function parse_cell_names(text) {
    let result = new Map();
    let lines = text.split(/\r?\n|;/g);
    lines.forEach((line) => {
        line = line.trim();
        if (line === "") {
            return;
        }
        let colonsplit = line.split(":");
        if (colonsplit.length != 2) {
            throw new Error(`couldn't parse "${line}"`);
        }
        let [dim_str, rest] = colonsplit;
        dim_str = dim_str.trim();
        if (!/^[-+]?\d+$/.test(dim_str)) {
            throw new Error(`couldn't parse dimension ${dim_str}`);
        }
        let dim = parseInt(dim_str);
        rest = rest.trim();
        if (rest.endsWith(",")) {
            rest = rest.slice(0, rest.length - 1).trim();
        }
        let names;
        if (rest === "") {
            names = [];
        }
        else {
            names = rest.split(",").map((s)=>s.trim());
            names.forEach((name) => {
                if (!face_name_regex.test(name)) {
                    throw new Error(`name ${name} is invalid`);
                }
            });
        }
        if (result.has(dim)) {
            throw new Error(`Multiple lines for dimension ${dim}`);
        }
        result.set(dim, names);
    });
    return result;
}

function parse_boundary(text) {
    // Example:
    // return new Map([
    //     ["x", [[1, "v"], [-1, "v"]]],
    //     ["y", [[1, "v"], [-1, "v"]]],
    //     ["F", [[1, "x"], [1, "y"], [1, "x"], [-1, "y"]]],
    // ]);
    let result = new Map();
    let lines = text.split(/\r?\n|;/g);
    lines.forEach((line) => {
        line = line.trim();
        if (line === "") {
            return;
        }
        let colonsplit = line.split(":");
        if (colonsplit.length != 2) {
            throw new Error(`couldn't parse "${line}"`);
        }
        let [name, rest] = colonsplit;
        name = name.trim();
        if (!face_name_regex.test(name)) {
            throw new Error(`name "${name}" is invalid`);
        }
        if (result.has(name)) {
            throw new Error(`duplicate boundary entry for ${name}`)
        }
        rest = rest.trim();
        let chain;
        if (rest === "" || rest === "0") {
            result.set(name, []);
            return;
        }
        if (rest[0] !== "+" && rest[0] !== "-") {
            rest = "+" + rest;
        }
        let splits = [0];
        for (let i = 1; i < rest.length; i++) {
            let c = rest.charAt(i);
            if (c === "+" || c === "-") {
                splits.push(i);
            }
        }
        splits.push(rest.length);
        chain = [];
        for (let j = 1; j < splits.length; j++) {
            let chunk = rest.slice(splits[j - 1], splits[j]);
            if (chunk === "") {
                throw new Error(`Duplicate +- in ${rest}`);
            }
            let sign = chunk.charAt(0) === "+" ? +1 : -1;
            let spaces = chunk.slice(1).search(/\S/);
            if (spaces === -1) {
                throw new Error("Empty face name");
            }
            let digits = chunk.slice(1 + spaces).search(/\D/);
            if (digits == -1) {
                throw new Error(`Integer ${chunk} is not a cell`);
            }
            let coeff = BigInt(sign);
            if (digits > 0) {
                coeff *= BigInt(chunk.slice(1 + spaces, 1 + spaces + digits));
            }
            let boundary_name = chunk.slice(1 + spaces + digits).trim();
            if (!face_name_regex.test(boundary_name)) {
                throw new Error(`name "${boundary_name}" is invalid`);
            }
            chain.push([coeff, boundary_name]);
        }
        result.set(name, chain);
    });
    return result;
}

function parse_relative(relative_text) {
    relative_text = relative_text.trim();
    if (relative_text.endsWith(",")) {
        relative_text = relative_text.slice(0, relative_text.length - 1).trim();
    }
    if (relative_text == "") {
        return [];
    }
    let names = relative_text.split(",").map((s=>s.trim()));
    let set_names = new Set();
    names.forEach((name) => {
        if (!face_name_regex.test(name)) {
            throw new Error(`name "${name}" is invalid`);
        }
        if (set_names.has(name)) {
            throw new Error(`Duplicate cell ${name} in relative subcomplex`)
        }
        set_names.add(name);
    });
    return names;
}

function parse_coeff(coeff_text) {
    if (coeff_text == "") {
        return 0n;
    }
    else if (/^(\d|\s|_)+$/.test(coeff_text)) {
        return BigInt(coeff_text.replace(/\s|_/g, ''));
    }
    else {
        throw new Error(`Couldn't parse modulus "${coeff_text}" as a nonnegative integer`);
    }
}

function get_homology_result_label_and_string(
    cell_names_text,
    boundary_text,
    co,
    coeff_text,
    relative_text,
) {
    let result = null;
    let coeff = null;
    let relative = null;
    try {
        coeff = parse_coeff(coeff_text);
        let cell_names = parse_cell_names(cell_names_text);
        let boundary = parse_boundary(boundary_text);
        relative = parse_relative(relative_text);
        result = homology_from_names(cell_names, boundary, co, coeff, relative);
    }
    catch (err) {
        console.error(err);
        return ["Result", "Error: " + err.message, false];
    }
    let label_text = (
        (relative.length > 0 ? "Relative " : "")
        + (co ? "Cohomology" : "Homology")
        + " with "
        + (coeff === 0n ? "ℤ" : `ℤ/${coeff.toString()}ℤ`)
        + " coefficients"
    );
    let string_result = [];
    result.forEach(([dim, {free_generators, torsion_generators}]) => {
        let head;
        if (co) {
            head = `H^${dim}:`
        }
        else {
            head = `H_${dim}:`
        }
        function print(s) {
            let line = head + " ".repeat(6 - head.length) + s;
            string_result.push(line);
            head = "";
        }
        if (free_generators.length == 0 && torsion_generators.length == 0) {
            print("trivial");
        }
        free_generators.forEach((freegen) => {
            print(`a copy of ℤ generated by ${coeffs_to_string(freegen)}`);
        });
        torsion_generators.forEach(([torgen, order]) => {
            print(`a copy of ℤ/${order.toString()}ℤ generated by ${coeffs_to_string(torgen)}`);
        });
        string_result.push("");
    });
    return [label_text, string_result.join("\r\n"), true];
}

let cell_names_textarea = document.getElementById("cell_names");
let boundary_textarea = document.getElementById("boundary");
let output_textarea = document.getElementById("output_textarea");
let auto_check = document.getElementById("auto_check");
let co_check = document.getElementById("co_check");
let coeff_input = document.getElementById("coeff_input");
let output_label = document.getElementById("output_label");
let relative_input = document.getElementById("relative");

function update() {
    // Update result
    output_textarea.value = "...";
    output_label.innerHTML = "Processing...";
    output_textarea.style.color = "black";
    // output_textarea.style.border = "none"
    output_label.style.color = "black";

    let [title_string, result_string, ok] = get_homology_result_label_and_string(
        cell_names_textarea.value,
        boundary_textarea.value,
        co_check.checked,
        coeff_input.value,
        relative_input.value,
    );
    output_textarea.value = result_string;
    output_label.innerHTML = title_string;

    let color = ok ? "#0000ff" : "#ff0000";
    output_textarea.style.color = color;
    output_label.style.color = color;
}


// https://stackoverflow.com/a/14029861/11461120
[
    cell_names_textarea,
    boundary_textarea,
    coeff_input,
    relative_input,
].forEach((element) => {
    if (element.addEventListener) {
        element.addEventListener("input", update, false);
    } else if (element.attachEvent) {
        element.attachEvent("onpropertychange", update);
    }
    else {
        alert("This browser cannot detect when text fields change?");
    }
})
co_check.addEventListener("change", update);

function do_example(ex_name) {
    let [ex_names, ex_boundary] = EXAMPLES.get(ex_name);
    cell_names_textarea.value = ex_names;
    boundary_textarea.value = ex_boundary;
    update();
}

co_check.checked = false;
coeff_input.value = "";
relative_input.value = "";
do_example("torus_delta");
