use ndarray::*;
use gnuplot::*;
use std::fs;
use std::io::Write;
use std::time::{Duration, SystemTime};
fn build_up_b(rho:f64,dt:f64,dx:f64,dy:f64,u:&Array2<f64>,v:&Array2<f64>)->Array2<f64>{
    let row = u.shape()[0];
    let col = u.shape()[1];
    let mut b = Array2::<f64>::zeros((row,col));
    for i in 1..row-1{
        for j in 1..col-1{
            b[[i,j]] = rho * dx.powi(2) * dy.powi(2)/2./(dx.powi(2)+dy.powi(2))
                        *(1./dt*(u[[i,j+1]]-u[[i,j-1]])/2./dx
                          +1./dt*(v[[i+1,j]]-v[[i-1,j]])/2./dy
                          -(u[[i,j+1]]-u[[i,j-1]])/2./dx*(u[[i,j+1]]-u[[i,j-1]])/2./dx
                          -2.*(u[[i+1,j]]-u[[i-1,j]])/2./dy*(v[[i,j+1]]-v[[i,j-1]])/2./dx
                          -(v[[i+1,j]]-v[[i-1,j]])/2./dy*(v[[i+1,j]]-v[[i-1,j]])/2./dy
                         );
        }
    }
    for i in 1..row-1{
        b[[i,col-1]] = rho * dx.powi(2) * dy.powi(2)/2./(dx.powi(2)+dy.powi(2))
                        *(1./dt*(u[[i,0]]-u[[i,col-2]])/2./dx
                          +1./dt*(v[[i+1,col-1]]-v[[i-1,col-1]])/2./dy
                          -(u[[i,0]]-u[[i,col-2]])/2./dx*(u[[i,0]]-u[[i,col-2]])/2./dx
                          -2.*(u[[i+1,col-1]]-u[[i-1,col-1]])/2./dy*(v[[i,0]]-v[[i,col-2]])/2./dx
                          -(v[[i+1,col-1]]-v[[i-1,col-1]])/2./dy*(v[[i+1,col-1]]-v[[i-1,col-1]])/2./dy
                         );
    }
    for i in 1..row-1{
        b[[i,0]] = rho * dx.powi(2) * dy.powi(2)/2./(dx.powi(2)+dy.powi(2))
                        *(1./dt*(u[[i,1]]-u[[i,col-1]])/2./dx
                          +1./dt*(v[[i+1,0]]-v[[i-1,0]])/2./dy
                          -(u[[i,1]]-u[[i,col-1]])/2./dx*(u[[i,1]]-u[[i,col-1]])/2./dx
                          -2.*(u[[i+1,0]]-u[[i-1,0]])/2./dy*(v[[i,1]]-v[[i,col-1]])/2./dx
                          -(v[[i+1,0]]-v[[i-1,0]])/2./dy*(v[[i+1,0]]-v[[i-1,0]])/2./dy
                         );
    }
    return b;
}

fn pressure_poisson_periodic(p:&Array2<f64>,dx:f64,dy:f64)->Array2<f64>{
    let mut p = p.clone();
    let row = p.shape()[0];
    let col = p.shape()[1];
    let mut pn = Array2::<f64>::zeros((row,col));
    let nit = 50;
    for q in 0..nit{
        pn = p.clone();
        for i in 1..row-1{
            for j in 1..col-1{
                p[[i,j]] = ((pn[[i,j+1]]+pn[[i,j-1]])*dy.powi(2)
                            +(pn[[i+1,j]]+pn[[i-1,j]])*dx.powi(2))
                            /2./(dx.powi(2)+dy.powi(2));
            }
            p[[i,col-1]] = ((pn[[i,0]]+pn[[i,col-2]])*dy.powi(2)
                            +(pn[[i+1,col-1]]+pn[[i-1,col-1]])*dx.powi(2))
                            /2./(dx.powi(2)+dy.powi(2));
            p[[i,0]] = ((pn[[i,1]]+pn[[i,col-1]])*dy.powi(2)
                            +(pn[[i+1,0]]+pn[[i-1,0]])*dx.powi(2))
                            /2./(dx.powi(2)+dy.powi(2));
        }
        for j in 0..col{
            p[[row-1,j]] = p[[row-2,j]];
            p[[0,j]] = p[[1,j]];
        }
    }
    return p;
}

fn main(){
    let now = SystemTime::now();
    let nx = 41;
    let ny = 41;
    let nt = 10;
    // let nit = 50;
    let c = 1;
    let dx = 2./(nx as f64 -1.);
    let dy = 2./(ny as f64 -1.);
    
    let rho = 1.;
    let nu = 0.1;
    let F = 1.;
    let dt = 0.01;

    let mut u = Array2::<f64>::zeros((ny,nx));
    let mut un = Array2::<f64>::zeros((ny,nx));

    let mut v = Array2::<f64>::zeros((ny,nx));
    let mut vn = Array2::<f64>::zeros((ny,nx));

    let mut p = Array2::<f64>::zeros((ny,nx));
    let mut pn = Array2::<f64>::zeros((ny,nx));
    
    let mut b = Array2::<f64>::zeros((ny,nx));

    let mut udiff = 1.;
    let mut stepcount = 0;

    while udiff > 0.001{
        un = u.clone();
        vn = v.clone();

        b = build_up_b(rho, dt, dx, dy, &u, &v);
        p = pressure_poisson_periodic(&p, dx, dy);

        for i in 1..ny-1{
            for j in 1..nx-1{
                u[[i,j]] = un[[i,j]] 
                            - un[[i,j]]*dt/dx*(un[[i,j]]-un[[i,j-1]])
                            - vn[[i,j]]*dt/dy*(un[[i,j]]-un[[i-1,j]])
                            - dt/rho/2./dx * (p[[i,j+1]]-p[[i,j-1]])
                            + nu*(dt/dx.powi(2)*(un[[i,j+1]]-2.*un[[i,j]]+un[[i,j-1]])
                                  +dt/dy.powi(2)*(un[[i+1,j]]-2.*un[[i,j]]+un[[i-1,j]])
                                 )
                            +dt*F;
                v[[i,j]] = vn[[i,j]] 
                            - un[[i,j]]*dt/dx*(vn[[i,j]]-vn[[i,j-1]])
                            - vn[[i,j]]*dt/dy*(vn[[i,j]]-v[[i-1,j]])
                            - dt/rho/2./dx * (p[[i+1,j]]-p[[i-1,j]])
                            + nu*(dt/dx.powi(2)*(vn[[i,j+1]]-2.*vn[[i,j]]+vn[[i,j-1]])
                                  +dt/dy.powi(2)*(vn[[i+1,j]]-2.*vn[[i,j]]+vn[[i-1,j]])
                                 );
            }
            u[[i,nx-1]] = un[[i,nx-1]] 
                            - un[[i,nx-1]]*dt/dx*(un[[i,nx-1]]-un[[i,nx-2]])
                            - vn[[i,nx-1]]*dt/dy*(un[[i,nx-1]]-un[[i-1,nx-1]])
                            - dt/rho/2./dx * (p[[i,0]]-p[[i,nx-2]])
                            + nu*(dt/dx.powi(2)*(un[[i,0]]-2.*un[[i,nx-1]]+un[[i,nx-2]])
                                  +dt/dy.powi(2)*(un[[i+1,nx-1]]-2.*un[[i,nx-1]]+un[[i-1,nx-1]])
                                 )
                            +dt*F;
            u[[i,0]] = un[[i,0]] 
                            - un[[i,0]]*dt/dx*(un[[i,0]]-un[[i,nx-1]])
                            - vn[[i,0]]*dt/dy*(un[[i,0]]-un[[i-1,0]])
                            - dt/rho/2./dx * (p[[i,1]]-p[[i,nx-1]])
                            + nu*(dt/dx.powi(2)*(un[[i,1]]-2.*un[[i,0]]+un[[i,nx-1]])
                                  +dt/dy.powi(2)*(un[[i+1,0]]-2.*un[[i,0]]+un[[i-1,0]])
                                 )
                            +dt*F;
            v[[i,nx-1]] = vn[[i,nx-1]] 
                            - un[[i,nx-1]]*dt/dx*(vn[[i,nx-1]]-vn[[i,nx-2]])
                            - vn[[i,nx-1]]*dt/dy*(vn[[i,nx-1]]-v[[i-1,nx-1]])
                            - dt/rho/2./dx * (p[[i+1,nx-1]]-p[[i-1,nx-1]])
                            + nu*(dt/dx.powi(2)*(vn[[i,0]]-2.*vn[[i,nx-1]]+vn[[i,nx-2]])
                                  +dt/dy.powi(2)*(vn[[i+1,nx-1]]-2.*vn[[i,nx-1]]+vn[[i-1,nx-1]])
                                 );
            v[[i,0]] = vn[[i,0]] 
                            - un[[i,0]]*dt/dx*(vn[[i,0]]-vn[[i,nx-1]])
                            - vn[[i,0]]*dt/dy*(vn[[i,0]]-v[[i-1,0]])
                            - dt/rho/2./dx * (p[[i+1,0]]-p[[i-1,0]])
                            + nu*(dt/dx.powi(2)*(vn[[i,1]]-2.*vn[[i,0]]+vn[[i,nx-1]])
                                  +dt/dy.powi(2)*(vn[[i+1,0]]-2.*vn[[i,0]]+vn[[i-1,0]])
                                 );
        }
        for j in 0..nx{
            u[[0,j]] = 0.;
            u[[ny-1,j]] = 0.;
            u[[0,j]] = 0.;
            u[[ny-1,j]] = 0.;
        }

        stepcount+=1;
        udiff = (u.sum()-un.sum())/u.sum();
    }
    println!("{:?}",now.elapsed());
    println!("{}",stepcount);

    let mut fg = Figure::new();
    fg.axes3d().surface(&u,nx,ny,Some((0.,0.,2.,2.)),&[]);
    fg.show();
    let mut fg = Figure::new();
    fg.axes3d().surface(&v,nx,ny,Some((0.,0.,2.,2.)),&[]);
    fg.show();
    let mut fg = Figure::new();
    fg.axes3d().surface(&p,nx,ny,Some((0.,0.,2.,2.)),&[]);
    fg.show();
    let mut fileu = fs::File::create("/home/zz/work/rust/cfd_tut/tut15_u.csv").unwrap();
    let mut filev = fs::File::create("/home/zz/work/rust/cfd_tut/tut15_v.csv").unwrap();
    let mut filep = fs::File::create("/home/zz/work/rust/cfd_tut/tut15_p.csv").unwrap();

    for i in 0..ny{
        for j in 0..nx{
            // println!("{},{},{},{},{}",i as f64*dx,j as f64*dy
            //     ,u[[i,j]],v[[i,j]],p[[i,j]]);
            let su:String=format!("{},",u[[i,j]]);
            fileu.write_all(su.as_bytes()).unwrap();
            let sv:String=format!("{},",v[[i,j]]);
            filev.write_all(sv.as_bytes()).unwrap();
            let sp:String=format!("{},",p[[i,j]]);
            filep.write_all(sp.as_bytes()).unwrap();
        }
        fileu.write_all(b"\n").unwrap();
        filev.write_all(b"\n").unwrap();
        filep.write_all(b"\n").unwrap();
    }
}
