use std::f64::consts;
use std::vec::Vec;

#[allow(dead_code)]
fn main() {
	let length: f64 = 2.0 * consts::PI; //length of the PIC box
	let np: i32 = 20; //number of particles
	let ng: i32 = 64; //number of zones for mesh
	let m1: i32 = 1; //mass of each particle in the first stream
	let m2: i32 = 1; //mass of each particle in the second stream
	let k: i32 = 1; //used for the peturbation of the sin wave
	let q1: i32 = -1; //charges of each of the particles
	let q2: i32 = -1;
	let vStream: f64 = 4.0; //initial velocity of each stream, this number will turned into positive and negative!
	let dt: f64 = 0.0; //timestep
	let deltaX: f64 = length/(np as f64); //L/np, or the delta x of the lattice
	let tMax: f64 = 16 as f64*consts::PI/((2 as f64*(np as f64)/length).sqrt()); //maximum time
	let mut velocity: Vec<f64> = Vec::new(); //velocities of each particle, index is 
	let mut xPosition: Vec<f64> = Vec::new(); //x-position of each particle
	
	//seed the initial positions and velocities
	for i in 0..(np*2) {
		if i < np {
			xPosition.push(i as f64*length/np as f64);
			velocity.push(-vStream);
		}
		else {
			xPosition.push((i-np) as f64*length/np as f64);
			velocity.push(vStream);
		}
		//rust only accepts usize for indexes
		println!("x: {}, v: {}", xPosition[i as usize], velocity[i as usize]);
	}
	
	//peturb the x positions
	perturb(&mut xPosition, length, k);
	printXPos(xPosition)
}

//evaluates the electron number density n(0:J-1) from l
//l is the array from [0, n-1], with 64 bit floats
fn density(l: &mut[f64]) {
	
}

//solves the 1D poisson equation
fn poisson1d(u: &mut[f64], v: &mut[f64], kappa: f64) {
	
}

//calculates the electric field from potential, where
//E is the potential
fn electricField(phi: &mut[f64], E: &mut[f64]) {

}

//gets called whenever a timestep is taken, extrapolates x and v for a particle for a single timestep
fn leapfrog() {

}

//calculates the cloud in cell weight for the force
fn cicWeight() {

}

//petrubs the function by a sine wave with wavelength = L/k and amplitude of .01K
fn perturb(xPosition: &mut Vec<f64>, length: f64, k: i32) {
	let wavelength: f64 = length/(k as f64);
	let amplitude: f64 = 0.01 * length;
	for i in 0..xPosition.len() {
		xPosition[i as usize] = xPosition[i as usize] + amplitude*(wavelength*xPosition[i as usize]).sin();
	}
}

//prints the list of X positions. THIS IS A COPY, NOT A DIRECT REFERENCE
fn printXPos(xPosition: Vec<f64>) {
	for i in 0..xPosition.len() {
		println!("particle: {}, x: {}", i, xPosition[i as usize])
	}
}