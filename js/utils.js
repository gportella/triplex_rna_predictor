function isNumeric(n) {
  return !isNaN(parseFloat(n)) && isFinite(n);
}

// generate a random id name
var rand_id = function () {
  // Math.random should be unique because of its seeding algorithm.
  // Convert it to base 36 (numbers + letters), and grab the first 9 characters
  // after the decimal.
  return '_' + Math.random().toString(36).substr(2, 9) + ".fna";
};


