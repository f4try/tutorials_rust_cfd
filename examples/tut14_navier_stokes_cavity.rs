use ndarray::*;
use gnuplot::*;
use std::fs;
use std::io::Write;
fn build_up_b(b:&mut Array2<f64>,rho:f64,dt:f64,u:&Array2<f64>,v:&Array2<f64>,dx:f64,dy:f64){
    let row = b.shape()[0];
    let col = b.shape()[1];
    Zip::from(b.slice_mut(s![1..row-1,1..col-1]))
        .and(u.slice(s![2..,1..col-1]))
        .and(u.slice(s![..row-2,1..col-1]))
        // .and(u.slice(s![1..row-2,2..]))
        // .and(u.slice(s![1..row-2,..col-1]))
        // .and(v.slice(s![2..,1..col-2]))
        // .and(v.slice(s![..row-1,1..col-2]))
        .and(v.slice(s![1..row-1,2..]))
        .and(v.slice(s![1..row-1,..col-2]))
        // .apply(|bij,&uip1j,&uim1j,&uijp1,&uijm1,&vip1j,&vim1j,&vijp1,&vijm1|
        .apply(|bij,&uip1j,&uim1j,&vijp1,&vijm1|*bij=
        rho*(1./dt*((uip1j-uim1j)/2./dx+(vijp1-vijm1)/2./dy)
            -(uip1j-uim1j)/2./dx*(uip1j-uim1j)/2./dx
            -2.*(uip1j-uim1j)/2./dx*(vijp1-vijm1)/2./dy
            -(vijp1-vijm1)/2./dy*(vijp1-vijm1)/2./dy
            )
        );
}
fn pressure_poisson(p:&mut Array2<f64>,dx:f64,dy:f64,b:&Array2<f64>,nit:usize,rho:f64){
    let mut pn = p.clone();
    let row = p.shape()[0];
    let col = p.shape()[1];
    for n in 1..nit{
        pn = p.clone();
        for i in 1..row-1{
            for j in 1..col-1{
                p[[i,j]] = ((pn[[i+1,j]]+pn[[i-1,j]])*dy.powi(2)
                    +(pn[[i,j+1]]+pn[[i,j-1]])*dx.powi(2))/2.
                    /(dx.powi(2)+dy.powi(2))
                    -rho*dx.powi(2)*dy.powi(2)/2.
                    /(dx.powi(2)+dy.powi(2))*b[[i,j]];
            }
        }
        for i in 0..row{
            p[[i,col-1]] = p[[i,col-2]];
            p[[i,0]] = p[[i,1]];
        }
        for j in 0..col{
            p[[0,j]]=p[[1,j]];
            p[[row-1,j]] = 0.;
        }
    }
}
fn cavity_flow(nt:usize,u:&mut Array2<f64>,v:&mut Array2<f64>
    ,dt:f64,dx:f64,dy:f64,nit:usize,p:&mut Array2<f64>,rho:f64,nu:f64){
    let mut un = u.clone();
    let mut vn = v.clone();
    let row = u.shape()[0];
    let col = u.shape()[1];
    let mut b = Array2::<f64>::zeros((row,col));
    for n in 0..nt{
        un = u.clone();
        vn = v.clone();

        build_up_b(&mut b, rho, dt, &u, &v, dx, dy);
        pressure_poisson(p, dx, dy, &b, nit, rho);
        for i in 1..row-1{
            for j in 1..col-1{
                u[[i,j]]= un[[i,j]]-un[[i,j]]*dt/dx*(un[[i,j]]-un[[i-1,j]])
                -vn[[i,j]]*dt/dx*(un[[i,j]]-un[[i,j-1]])
                -dt/rho/2./dx*(p[[i+1,j]]-p[[i-1,j]])
                +nu*dt/dx.powi(2)*(un[[i+1,j]]-2.*un[[i,j]]+un[[i-1,j]])
                +nu*dt/dy.powi(2)*(un[[i,j+1]]-2.*un[[i,j]]+un[[i,j-1]]);
                v[[i,j]]= vn[[i,j]]-un[[i,j]]*dt/dx*(vn[[i,j]]-vn[[i-1,j]])
                -vn[[i,j]]*dt/dx*(vn[[i,j]]-vn[[i,j-1]])
                -dt/rho/2./dx*(p[[i,j+1]]-p[[i,j-1]])
                +nu*dt/dx.powi(2)*(vn[[i+1,j]]-2.*vn[[i,j]]+vn[[i-1,j]])
                +nu*dt/dy.powi(2)*(vn[[i,j+1]]-2.*vn[[i,j]]+vn[[i,j-1]]);
            }
        }

        for i in 1..row-1{
            u[[i,0]] = 0.;
            u[[i,col-1]]=0.;
            v[[i,0]] = 0.;
            v[[i,col-1]]=0.;
        }
            
        for j in 1..col-1{
            u[[0,j]] = 0.;
            u[[row-1,j]]=1.;
            v[[0,j]] = 0.;
            v[[row-1,j]]=0.;
        }
    }
}
fn main(){
    let nx = 41;
    let ny = 41;
    let nt = 700;
    let nit = 50;
    let c = 1.;
    let dx = 2./(nx-1) as f64;
    let dy = 2./(ny-1) as f64;
    let rho = 1.;
    let nu = 0.1;
    let dt = 0.001;
    let mut u = Array2::<f64>::zeros((ny,nx));
    let mut v = Array2::<f64>::zeros((ny,nx));
    let mut p = Array2::<f64>::zeros((ny,nx));
    let mut b = Array2::<f64>::zeros((ny,nx));
    cavity_flow(nt, &mut u, &mut v, dt, dx, dy, nit, &mut p, rho, nu);
    let mut fg = Figure::new();
    fg.axes3d().surface(&u,nx,ny,Some((0.,0.,2.,2.)),&[]);
    fg.show();
    let mut fg = Figure::new();
    fg.axes3d().surface(&v,nx,ny,Some((0.,0.,2.,2.)),&[]);
    fg.show();
    let mut fg = Figure::new();
    fg.axes3d().surface(&p,nx,ny,Some((0.,0.,2.,2.)),&[]);
    fg.show();
    let mut fileu = fs::File::create("/home/zz/work/rust/cfd_tut/tut14_u.csv").unwrap();
    let mut filev = fs::File::create("/home/zz/work/rust/cfd_tut/tut14_v.csv").unwrap();
    let mut filep = fs::File::create("/home/zz/work/rust/cfd_tut/tut14_p.csv").unwrap();

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