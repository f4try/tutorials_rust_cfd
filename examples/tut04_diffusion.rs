use ndarray::*;
use gnuplot::*;

fn main(){
    let nx = 41;
    let dx = 2./(nx-1) as f64;
    let nt = 20;
    let sigma=0.2;
    let nu = 0.3;
    let dt = sigma*dx.powi(2)/nu;
    let mut u = Array1::ones(nx);
    Zip::from(u.slice_mut(s![(0.5/dx) as usize..(1.0/dx) as usize])).apply(|ue|*ue=2.);
    // println!("{:?}",u);
    let mut un = Array1::ones(nx);
    let mut fg=Figure::new();
    fg.axes2d().lines_points(&Array1::linspace(0.,2.,nx),&u,&[]);
    for n in 1..nt{
        for i in 1..nx-1{
            un[i]=u[i]+nu*dt/dx.powi(2)*(u[i+1]-2.*u[i]+u[i-1]);
        }
        fg.axes2d().lines_points(&Array1::linspace(0.,2.,nx),&un,&[]);
        u = un.clone();
    }
    fg.show();
}