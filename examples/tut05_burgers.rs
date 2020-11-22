// use ndarray::*;
use core::f64::consts::PI;
use gnuplot::*;
fn analytical_u(x:f64,t:f64,nu:f64)->f64{
    -2.*nu*(-(-8.*t + 2.*x)*(-(-4.*t + x).powi(2)/(4.*nu*(t + 1.))).exp()/(4.*nu*(t + 1.)) - (-8.*t + 2.*x - 4.*PI)*(-(-4.*t + x - 2.*PI).powi(2)/(4.*nu*(t + 1.))).exp()/(4.*nu*(t + 1.)))/((-(-4.*t + x - 2.*PI).powi(2)/(4.*nu*(t + 1.))).exp() + (-(-4.*t + x).powi(2)/(4.*nu*(t + 1.))).exp()) + 4.
}
pub fn linspace(start:f64,end:f64,n:usize)->Vec<f64>{
    let mut result = vec![start;n];
    let interval = (end-start)/n as f64;
    for i in 1..n{
        result[i] += interval*i as f64;
    }
    result
}
fn main(){
    let nx = 101;
    let nt = 100;
    let dx = 2.*PI/(nx-1) as f64;
    let sigma = 0.5;
    let nu = 0.07;
    // let dt = sigma*dx.powi(2)/nu;
    let dt = dx * nu;
    let t = 0;
    let mut u = vec![0.;nx];
    for i in 0..nx{
        u[i] = analytical_u(i as f64*dx, 0., nu);
    }
    // println!("{:?}",u);
    let mut un = u.clone();
    let mut fg = Figure::new();

    // fg.axes2d().lines_points(&Array1::linspace(0.,2.*PI,nx),&u,&[]);
    fg.axes2d().lines_points(&linspace(0.,2.*PI,nx),&u,&[]);
    for n in 1..nt{
        for i in 1..(nx-1){
            un[i] = u[i]-u[i]*dt/dx*(u[i]-u[i-1])+nu*dt/dx.powi(2)*(u[i+1]-2.*u[i]+u[i-1]);            
        }
        un[0] = u[0]-u[0]*dt/dx*(u[0]-u[nx-2])+nu*dt/dx.powi(2)*(u[1]-2.*u[0]+u[nx-2]);            
        un[nx-1] = un[0];
        // fg.axes2d().lines_points(&Array1::linspace(0.,2.*PI,nx),&un,&[]);
        u = un.clone();
    }
    fg.axes2d().lines_points(&linspace(0.,2.*PI,nx),&un,&[]);
    fg.show();
}