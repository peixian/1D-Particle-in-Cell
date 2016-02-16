fn main() {
    println!("Hello, world!");
}

//evaluates the electron number density n(0:J-1) from l
//l is the array from [0, n-1], with 64 bit floats
fn density(l: &[f64]) {
	
}

//solves the 1D poisson equation
fn poisson1d(u: &[f64], v: &[f64], kappa: f64) {

}

//calculates the electric field from potential, where
//E is the potential
fn electricField(phi: &[f64], E: &[f64]) {

}

//gets called whenever a timestep is taken, extrapolates x and v for a particle for a single timestep
fn leapfrog() {

}

//calculates the cloud in cell weight for the force
fn cicWeight() {

}