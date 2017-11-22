
// ideal gas constant
var R = 0.0019858;

/** Function that count occurrences of a substring in a string;
 * @param {String} string               The string
 * @param {String} subString            The sub string to search for
 * @param {Boolean} [allowOverlapping]  Optional. (Default:false)
 *
 * @author Vitim.us https://gist.github.com/victornpb/7736865
 * @see Unit Test https://jsfiddle.net/Victornpb/5axuh96u/
 * @see http://stackoverflow.com/questions/4009756/how-to-count-string-occurrence-in-string/7924240#7924240
 */
function occurrences(string, subString, allowOverlapping) {

    string += "";
    subString += "";
    if (subString.length <= 0) return (string.length + 1);

    var n = 0,
        pos = 0,
        step = allowOverlapping ? 1 : subString.length;

    while (true) {
        pos = string.indexOf(subString, pos);
        if (pos >= 0) {
            ++n;
            pos += step;
        } else break;
    }
    return n;
}

function count_nuc(sequence, bMixed){
	/* "This function returns a
        matrix of CC, TC/CT, TT composition if bMixed is false, 
		and a matrix of C, U and CC composition if bMixed is set to true.
    */

    if (bMixed == false){
	    A = nj.zeros([1,3]);
        A.set(0, 0, occurrences(sequence, "CC", true));
        A.set(0, 1, occurrences(sequence, "UC", true) + occurrences(sequence, "CU", true));
        A.set(0, 2, occurrences(sequence, "UU", true));
    } else {
	    A = nj.zeros([1,4]);
        A.set(0, 0, occurrences(sequence, "C", true));
        A.set(0, 1, occurrences(sequence, "U", true));
        A.set(0, 2, occurrences(sequence, "CC", true));
        A.set(0, 3, 1);
    }
	return A;
}

function apply_ph_correct(om_dg, p_dg, x_ph, cond){
    /* Apply corrections to DG based on the number of C and CC
        The corrections formula is C*(pH-5.6)*(a+b*CC), where C and CC are the
        number of C and CC in the sequence, entries om_dg[0] and om_dg[2]

        :rtype: numpy matrix with corrected dg
    */

    rdg = p_dg + om_dg.get(0,0) * (cond.get(0) - 5.6)*(x_ph.get(0)+x_ph.get(1)*om_dg.get(0,2));
    return rdg;
}

function compute_tm(pred_dh, pred_dg, cond){
    /*  Compute the Tm based on DG, DH and the concentration
        tm = (298*DH / ( DH-DG - 298*R*log(4/Ct)))
        the Ct (concentration) is the second record in cond_mat

    	:rtype: numpy matrix with Tm estimates
    */
    ref_t = 298
    tm = (ref_t*pred_dh /
            (pred_dh-pred_dg - 298*R*Math.log(4/cond.get(1)/1e-6)));

    return tm;

}

function compute_c50(pred_dg, cond){
    /* Compute concentration of TFO needed to achieve 50% of
        triplex formation. If I am not mistaken, that should go down
        like this, where names of species imply concentration:
        triplex@c50 = duplex_initial / 2  == d/2
        Oo = initial TFO conc
        d/2 / (d/2*Oo-d/2) = K_eq
        Solving for Oo we have
        Oo = (1/K_eq) + d/2 = exp(DG/kT) + d/2
    */

    ref_t = 310
    // we use uM units, hence the 1e6 factor
    oligo_c50 = 1e6*Math.exp(pred_dg/(ref_t*R)) + 0.5 * cond.get(2);

    return oligo_c50;

}

function compute_tfo_tm(tfo_cond){
    /* Apply the nearest-neighbour (NN) model to predict DG and DH for
     * triplex formation ((C,T)-triplex only). The stability has a very
     * strong pH dependence. From the DG and DH we can calculate the Tm as
     * function of the concentration of oligo and also the C_50
     *
     * Read in the sequence and the conditions (pH, oligo and duplex
     * concentration) and apply the NN model to get the DG and DH.  Then apply
     * a pH correction (only to DG), and from there get the Tm and C_50. 
     */

    sequence = tfo_cond.sequence var cond = nj.array([tfo_cond.pH,
            tfo_cond.tfo_conc, tfo_cond.dup_conc]); 
    // this are part of the model
    var nn_dh = nj.array([-10.95, -5.73, -6.44]); 
    var nn_dg = nj.array([-1.891, -0.758, -0.331, 2.646]); 
    var x_ph = nj.array([0.893, -0.005]);

    // bMixed switches between counting nuc. or dinuc. for the DH
    // model or for the DG model.
    var bMixed = false;
    om_dh = count_nuc(sequence, bMixed);
    var bMixed = true;
    om_dg = count_nuc(sequence, bMixed);
    // mulitply coeficients and sequence composition to get the prediction
    pred_dh = om_dh.dot(nn_dh).get(0);
    pred_dg = om_dg.dot(nn_dg).get(0);

    // correct for the pH
    corrected_dg = apply_ph_correct(om_dg, pred_dg, x_ph, cond);
    // predict Tm and C50
    pred_tm = compute_tm(pred_dh, corrected_dg, cond);
    c50 = compute_c50(corrected_dg, cond);
    // output as a JS object, stringified just because.
    var results = {
            "DH": pred_dh.toFixed(1),
            "DG": corrected_dg.toFixed(1),
            "Tm": (pred_tm-273.15).toFixed(1),
            "C50": c50.toFixed(3)
        }

    var obj = JSON.stringify(results);
    return obj;
}


var file_ready = false;
var file_name = "";

function isNumeric(n) {
	return !isNaN(parseFloat(n)) && isFinite(n);
}

// generate a random id name
var rand_id = function() {
	// Math.random should be unique because of its seeding algorithm.
	// Convert it to base 36 (numbers + letters), and grab the first 9 characters
	// after the decimal.
	return '_' + Math.random().toString(36).substr(2, 9) + ".fna";
};

// Event handling
document.addEventListener("DOMContentLoaded",
	function(event) {
		document.getElementById('gen_sample_input')
			.addEventListener("click", function() {
				document.getElementById("frm")
					.elements["seq"].defaultValue = "CUUCUCUCUUUUCCU";
				document.getElementById("frm")
					.elements["pH"].defaultValue = 7.2;
				document.getElementById("frm")
					.elements["tfo_conc"].defaultValue = 5;
				document.getElementById("frm")
					.elements["dup_conc"].defaultValue = 10;
			});
	}
);

window.onload = function() {
	document.getElementById("frm")
		.elements["seq"].placeholder = "string of Cs and Us";
	document.getElementById("frm")
		.elements["pH"].placeholder = "pH range from 4.5 to 8";
	document.getElementById("frm")
		.elements["tfo_conc"].placeholder = "e.g. 5";
	document.getElementById("frm")
		.elements["dup_conc"].placeholder = "e.g. 10";

};

document.addEventListener("DOMContentLoaded",
	function(event) {
		document.body
			.addEventListener("uploadReady", function() {
				document.getElementById("submit_job")
					.style.background = "#cb4b16";
				file_ready = true;
			});
	}
);

function check_input_simple_pred(tfo_cond) {
	var message = "ok";
	if (!isNumeric(tfo_cond.pH) || tfo_cond.pH < 4.5 || tfo_cond.pH > 8) {
		return "wrong_ph";
	}
	if (!isNumeric(tfo_cond.tfo_conc)) {
		return "wrong_conc";
	}
	var re = /[^(U|C)]/
	if (re.test(tfo_cond.sequence)) {
		return "wrong_seq"
	}
	return message;
};

function read_conditions() {

	var seq = document.getElementById("frm")
		.elements["seq"].value;
	var pH = document.getElementById("frm")
		.elements["pH"].value;
	var tfo_conc = document.getElementById("frm")
		.elements["tfo_conc"].value;
	var dup_conc = document.getElementById("frm")
		.elements["dup_conc"].value;
	var conditions = {
		"sequence": seq,
		"pH": pH,
		"tfo_conc": tfo_conc,
		"dup_conc": dup_conc
	};

	return conditions;

}

function print_why_wrong_input(check_input) {
	// print out an error message with reasons
	// why the input is wrong

	if (check_input === "wrong_ph") {
		document.querySelector("#content")
			.innerHTML =
			"<h5 id='results'>Error: pH out of range [4.5,8] or non-numeric </h5>";
	};
	if (check_input === "wrong_conc") {
		document.querySelector("#content")
			.innerHTML =
			"<h5 id='results'>Error: value of concentration is non-numeric </h5>";
	};
	if (check_input === "wrong_seq") {
		document.querySelector("#content")
			.innerHTML =
			"<h5 id='results'>Error: sequence should only " +
			"contain Us or Cs (capital letters)</h5>";
	};
}


function tfo_thermo_one_seq(tfo_conditions){
    result = compute_tfo_tm(tfo_conditions);
    output_pred = JSON.parse(result);
    document.querySelector("#content")
        .innerHTML =
        "<h5 id='results'>Estimated thermodynamic properties </h5>" +
        "<p>" +
        "Sequence: " + tfo_conditions.sequence + "<br>" +
        "pH: " + tfo_conditions.pH + "<br>" +
        "TFO conc.: " + tfo_conditions.tfo_conc + "mM<br>" +
        "Duplex conc.: " + tfo_conditions.dup_conc + "mM<br>" +
        "</p>" +
        "<ul> " +
        "<li>DG: " + output_pred.DG + " kcal/mol </li>" +
        "<li>DH: " + output_pred.DH + " kcal/mol </li>" +
        "<li>Tm: " + output_pred.Tm + " C </li>" +
        "<li>C50: " + output_pred.C50 + " mM </li>" +
        "</ul>";
    EPPZScrollTo.scrollVerticalToElementById('results', 20);

}


// Where the action happens
document.addEventListener("DOMContentLoaded",
	function(event) {

		// Unobtrusive event binding
		document.getElementById('submit_job')
			.addEventListener("click", function() {
				// here you should validate the input and if it is
				// ok then call the calculations. No more error checking
                // going on... TODO: error checking
				var tfo_conditions = read_conditions();
				var check_input = check_input_simple_pred(tfo_conditions);
				if (check_input !== "ok") {
					print_why_wrong_input(check_input);
                } else {
                    tfo_thermo_one_seq(tfo_conditions);
                }
			});
	}
);
