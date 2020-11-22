use ndarray::*;
use gnuplot::*;
use std::time::{Duration, SystemTime};
fn range(start:f64,end:f64,interval:f64)->Vec<f64>{
    let n = ((end-start)/interval) as usize;
    let mut result = vec![start;n];
    for i in 0..n{
        result[i]+=i as f64*interval;
    }
    result
}
// impl std::ops::Sub for Vec<f64>{
//     type Output = Self;
//     fn sub(self,other:Self)->Self{
//         let result = vec![0.;self.len()];
//         for i in 0..self.len(){
//             result[i]=self[i]-other[i];
//         }
//         result
//     }
// }
fn main(){
    let u = Array1::range(0.,6.,1.);
    println!("{:?}",u);
    let du = &u.slice(s![1..])-&u.slice(s![..u.len()-1]);
    println!("{:?}",du);
    let u = range(0., 6., 1.);
    println!("{:?}",u);
    let mut du = vec![0.;u.len()-1];
    for i in 0..du.len(){
        du[i]=u[i+1]-u[i];
    }
    println!("{:?}",du);

    let nx = 81;
    let ny = 81;
    let nt = 100;
    let c = 1.;
    let dx = 2./(nx-1) as f64;
    let dy = 2./(ny-1) as f64;
    let sigma = 0.2;
    let dt = sigma*dx;

    let x = Array1::linspace(0.,2.,nx);
    let y = Array1::linspace(0.,2.,ny);
    let mut u = Array2::<f64>::ones((ny,nx));
    let mut un = Array2::<f64>::ones((ny,nx));
    let uh = u.slice_mut(s![(0.5/dy) as usize..(1./dy)as usize+1,(0.5/dx) as usize..(1./dx)as usize+1]);
    Zip::from(uh).apply(|ue|*ue=2.);
    let mut fg = Figure::new();
    fg.axes3d().surface(&u,nx,ny,Some((0.,0.,2.,2.)),&[]);
    // println!("{:?}",u);
    let now = SystemTime::now();
    for n in 0..nt+1{
        un = u.clone();
        let row=u.shape()[0];
        let col=u.shape()[1];
        for j in 1..row{
            for i in 1..col{
                u[[j,i]] = un[[j,i]]-(c*dt/dx*(un[[j,i]]-un[[j,i-1]]))-(c*dt/dy*(un[[j,i]]-un[[j-1,i]]));
            }
            Zip::from(u.slice_mut(s![0,..])).apply(|ue|*ue=1.);
            Zip::from(u.slice_mut(s![-1,..])).apply(|ue|*ue=1.);
            Zip::from(u.slice_mut(s![..,0])).apply(|ue|*ue=1.);
            Zip::from(u.slice_mut(s![..,-1])).apply(|ue|*ue=1.);
        }
            // fg.axes3d().surface(&u,nx,ny,Some((0.,0.,2.,2.)),&[]);
    }
    println!("{:?}",now.elapsed());
    // fg.axes3d().surface(&u,nx,ny,Some((0.,0.,2.,2.)),&[]);
    // fg.show();
    let now = SystemTime::now();
    let mut u = Array2::<f64>::ones((ny,nx));
    let mut un = Array2::<f64>::ones((ny,nx));
    let uh = u.slice_mut(s![(0.5/dy) as usize..(1./dy)as usize+1,(0.5/dx) as usize..(1./dx)as usize+1]);
    Zip::from(uh).apply(|ue|*ue=2.);
    let now = SystemTime::now();
    for n in 1..nt+1{
        un = u.clone();
        Zip::from(u.slice_mut(s![1..,1..]))
            .and(un.slice(s![1..,1..]))
            .and(un.slice(s![1..,0..nx-1]))
            .and(un.slice(s![0..nx-1,1..]))
            .apply(|ue, &unij, &uni_1j, &unij_1|*ue=unij-(c*dt/dx*(unij-uni_1j))-c*dt/dy*(unij-unij_1));
        Zip::from(u.slice_mut(s![0,..])).apply(|ue|*ue=1.);
        Zip::from(u.slice_mut(s![-1,..])).apply(|ue|*ue=1.);
        Zip::from(u.slice_mut(s![..,0])).apply(|ue|*ue=1.);
        Zip::from(u.slice_mut(s![..,-1])).apply(|ue|*ue=1.);
        fg.axes3d().surface(&u,nx,ny,Some((0.,0.,2.,2.)),&[]);
    }

    println!("{:?}",now.elapsed());
    fg.axes3d().surface(&u,nx,ny,Some((0.,0.,2.,2.)),&[]);
    fg.show();
}