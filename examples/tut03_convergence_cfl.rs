use ndarray::{Array1,s,Zip};
use gnuplot::*;
macro_rules! printVar {
    ($v:expr) => {
        println!("{:?}",$v);
    };
}
fn linearconv(nx:usize){
    // let nx = 41;
    let dx = 2./(nx-1) as f64;
    let nt = 25;
    let dt = 0.025;
    let c = 1.;

    let mut u = Array1::<f64>::ones(nx);
    let mut uh = u.slice_mut(s![(0.5/dx) as usize..(1./dx +1.) as usize]);
    Zip::from(uh).apply(|ue|*ue=2.);
    printVar!(u);
    let mut fg1=Figure::new();
    fg1.axes2d().lines_points(&Array1::linspace(0.,2.,nx),&u,&[Caption("u0")]);

    let mut un = Array1::<f64>::ones(nx);
    for n in 1..nt{
        for i in 1..nx{
            un[i] = u[i]-c*dt/dx*(u[i]-u[i-1]);
            // un[i] = (u[i]/dt+un[i-1]*c/dx)/(1./dt+c/dx);
        } 
        fg1.axes2d().lines_points(&Array1::linspace(0.,2.,nx),&un,&[Caption("u0")]);
        u = un.clone();
    }
    printVar!(un);
    // fg1.show().unwrap();

    let mut fg2=Figure::new();
    fg2.axes2d().lines_points(&Array1::linspace(0.,2.,nx),&un,&[Caption("u0")]);
    fg2.show().unwrap();
}

fn linearconv_adapt(nx:usize){
    let dx = 2./(nx-1) as f64;
    // let dt = 0.025;
    let c = 1.;
    let sigma = 0.5;
    let dt = sigma*dx;
    let nt = (0.625/dt) as usize;
    let mut u = Array1::<f64>::ones(nx);
    let mut uh = u.slice_mut(s![(0.5/dx) as usize..(1./dx +1.) as usize]);
    Zip::from(uh).apply(|ue|*ue=2.);
    printVar!(u);
    let mut fg1=Figure::new();
    fg1.axes2d().lines_points(&Array1::linspace(0.,2.,nx),&u,&[Caption("u0")]);

    let mut un = Array1::<f64>::ones(nx);
    for n in 1..nt{
        for i in 1..nx{
            un[i] = u[i]-c*dt/dx*(u[i]-u[i-1]);
            // un[i] = (u[i]/dt+un[i-1]*c/dx)/(1./dt+c/dx);
        } 
        fg1.axes2d().lines_points(&Array1::linspace(0.,2.,nx),&un,&[Caption("u0")]);
        u = un.clone();
    }
    printVar!(un);
    // fg1.show().unwrap();

    let mut fg2=Figure::new();
    fg2.axes2d().lines_points(&Array1::linspace(0.,2.,nx),&un,&[Caption("u0")]);
    fg2.show().unwrap();
}

fn main(){
    linearconv(41);
    linearconv(61);
    linearconv(71);
    linearconv(81);
    linearconv(82);

    linearconv_adapt(41);
    linearconv_adapt(61);
    linearconv_adapt(71);
    linearconv_adapt(81);
    linearconv_adapt(82);
    linearconv_adapt(1000);
    linearconv_adapt(2000);
}