use std::f64::consts;
use std::vec::Vec;

const LENGTH: f64 = 2.0*consts::PI; //length of PIC box
const NG: i32 = 64; //number of zones for mesh
const MASS: i32 = 1; //mass of each particle
const CHARGE: i32 = -1; //charge of each particle
const V_STREAM: f64 = 4.0; //initial velocity of each stream, this number will turned into positive and negative!

fn main() {
	let particle_loading: i32 = 4; //particle loading
	let np: i32 = particle_loading*NG; //number of particles
	let k: i32 = 1; //used for the peturbation of the sin wave
	let dt: f64 = LENGTH/(NG as f64)/(4.0); //timestep
	let delta_x: f64 = LENGTH/(np as f64); //L/np, or the delta x of the lattice
	let t_max: f64 = 16 as f64*consts::PI/((2 as f64*(np as f64)/LENGTH).sqrt()); //maximum time
	let mut velocity: Vec<f64> = Vec::new(); //velocities of each particle, index is 
	let mut x_position: Vec<f64> = Vec::new(); //x-position of each particle
	
	//seed the initial positions and velocities
	for i in 0..(np*2) {
		if i < np {
			x_position.push(i as f64*LENGTH/np as f64);
			velocity.push(-V_STREAM);
		}
		else {
			x_position.push((i-np) as f64*LENGTH/np as f64);
			velocity.push(V_STREAM);
		}
		//rust only accepts usize for indexes
		println!("x: {}, v: {}", x_position[i as usize], velocity[i as usize]);
	}
	
	//peturb the x positions
	perturb(&mut x_position, k);
	let mut rho: Vec<f64> = vec![0.0; np as usize]; //mesh densities!
	let mut phi: Vec<f64> = vec![0.0; NG as usize]; //vector of phi's to solve for the poisson 
	let mut electric_mesh: Vec<f64> = vec![0.0; NG as usize]; //vector of electric field from potential

}

//evaluates the electron number density using CIC weight.
//Each mesh loops through all the particles to add in their contribution, calculates the electron density for that mesh block alone.
fn density(rho: &mut Vec<f64>, x_position: & Vec<f64>, np: i32) {
	for rho_i in 0..rho.len() {
		let mut sum: f64 = 0.0;
		for x_p in 0..np*2 {
			sum += (MASS as f64)*cic_weight(rho_i as f64 - x_position[x_p as usize], np);
		}
		rho[rho_i as usize] = 1.0/(LENGTH/NG as f64)*sum;
	}
}

//Cloud in Cell weighting
fn cic_weight(x: f64, np: i32) -> f64{
	let delta_x: f64 = LENGTH/(np as f64);
	if x > -delta_x {
		return 1.0 + x/delta_x;
	}
	else if x < delta_x {
		return 1.0 - x/delta_x;
	}
	else {
		return 0.0;
	}
}


//solves the 1D poisson equation
fn poisson1d(phi: &mut Vec<f64>, rho: &Vec<f64>) {
	for i in 0..NG {
		phi.push(0.0);
	}
	let dx_g: f64 = LENGTH as f64/NG as f64;
	let error_tolerance: f64 = 0.000001; //10^-6 
	let source = source_l2_norm(&rho, dx_g);
	let residual = residual_l2_norm(&rho, &phi, dx_g);
}

//calculates the L2Norm for the source
//probably don't want to call this outside of the poisson solver?
fn source_l2_norm(rho: &Vec<f64>, dx_g: f64) -> f64 {
	let mut source: f64 = 0.0;
	//sum the entire array
	for i in 0..rho.len() {
		source += rho[i as usize].powi(2);
	}
	//divide it by 1/Ng
	source /= dx_g;
	//then square root it
	source = source.sqrt();
	return source;
}

//calculates the L2Norm for the residual
fn residual_l2_norm(rho: &Vec<f64>, phi: &Vec<f64>, dx_g: f64) -> f64 {
	let mut residual: f64 = 0.0;
	
	for i in 0..rho.len() {
		residual += rho[i as usize];
		if i == 0 {
			residual += phi[i+1 as usize] - 2.0*phi[i as usize] + phi[(rho.len()-1) as usize];
		}
		else if i+1 == LENGTH as usize {
			residual += phi[0] - 2.0*phi[i as usize] + phi[(i-1) as usize];
		}
		else {
			residual += phi[(i+1 as usize)] - 2.0*phi[i as usize] + phi[(i-1) as usize];
			residual /= (LENGTH as f64/NG as f64).powi(2);
			residual = residual.sqrt();
		}
		
	}
	residual /= dx_g;
	residual = residual.sqrt();
	return residual;
}


//calculates the electric field from potential, where
//E is the potential
//uses the second order finite difference
fn electric_field (phi: &mut Vec<f64>, electric_mesh: & mut Vec<f64>) {
	//get dx_g
	let dx: f64 = LENGTH/(NG as f64);
	//set the iteration length for easier reference
	let iter_length = phi.len();
	for i in 1..iter_length-1 {
		//finite difference
		let phi_next = phi[i+1 as usize];
		let phi_prev = phi[i-1 as usize];
		electric_mesh[i] = (phi_next - 2.0*phi[i as usize] + phi_prev)/(dx.powi(2));
	}
	//particles at the edge need to wrap around, so wrap them around
	electric_mesh[0] = (phi[1] - 2.0*phi[0] + phi[iter_length-1])/(dx.powi(2));
	electric_mesh[iter_length-1] = (phi[0] - 2.0*phi[iter_length-1] + phi[iter_length-2])/(dx.powi(2));
}


//gets called whenever a timestep is taken, extrapolates x and v for a particle for a single timestep
fn leapfrog(electric_mesh: &Vec<f64>, x_position: &mut Vec<f64>, velocity: &mut Vec<f64>, dt: f64, delta_x: f64) {
	for i in 0..x_position.len() {
		let (x_new, v_new) = leap(x_position[i as usize], delta_x, velocity[i as usize], dt, electric_mesh);
		x_position[i as usize] = x_new;
		velocity[i as usize] = v_new;
	}
}

//calculates the drift/kick/drift for a specific value
//returns the new X and new V
//THIS SHOULD *ONLY* BE CALLED FROM THE LEAPFROG METHOD
fn leap(x: f64, delta_x: f64, v: f64, dt: f64, electric_mesh: &Vec<f64>) -> (f64, f64) {
	let j: i32 = (x/delta_x).floor() as i32;
	let y: f64 = x/delta_x - j as f64;
	//drift
	let x_half: f64 = x + 0.5*v*dt;
	//kick
	let v_new: f64 = v + accel(j, y, x_half, electric_mesh)*dt;
	//drift
	let x_new: f64 = x_half + 0.5*v_new*dt;
	return (x_new, v_new)
}

//calculates the accleration
fn accel(j: i32, y: f64, x_value: f64, electric_mesh: &Vec<f64>) -> f64{
	let mut dphi_dx: f64 = 0.0;
	if j+1 == NG {
		dphi_dx = electric_mesh[j as usize]*(1.0-y) + electric_mesh[0]*y;
	}
	else {
		dphi_dx = electric_mesh[j as usize]*(1.0-y) + electric_mesh[(j+1) as usize]*y;
	}
	let a: f64 = -CHARGE as f64/MASS as f64*dphi_dx;
	return a;
}


//petrubs the function by a sine wave with waveLENGTH = L/k and amplitude of .01K
fn perturb(x_position: &mut Vec<f64>, k: i32) {
	let wavelength: f64 = LENGTH/(k as f64);
	let amplitude: f64 = 0.01 * LENGTH;
	for i in 0..x_position.len() {
		x_position[i as usize] = x_position[i as usize] + amplitude*(wavelength*x_position[i as usize]*2.0*consts::PI).sin();
	}
}

//prints the list of X positions. THIS IS A COPY, NOT A DIRECT REFERENCE
fn print_vec(vector: &Vec<f64>) {
	for i in 0..vector.len() {
		println!("{}", vector[i as usize])
	}
}