use ndarray::*;
use gnuplot::*;

fn plot2D(p:&Array2<f64>,nx:usize,ny:usize,x:f64,y:f64){
    let mut fg = Figure::new();
    fg.axes3d().surface(p,nx,ny,Some((0.,0.,x,y)),&[]);
    fg.show();
}

fn abs(array:&Array2<f64>)->Array2<f64>{
    let mut result = array.clone();
    Zip::from(&mut result).apply(|v|*v=v.abs());
    result
}

fn poisson2d(p:&mut Array2<f64>,b:&Array2<f64>,dx:f64,dy:f64,l1norm_target:f64){
    let mut l1norm =1.;
    let row = p.shape()[0];
    let col = p.shape()[1];
    // let mut pn = Array2::zeros((p.shape()[0],p.shape()[1]));
    let mut pn = p.clone();
    while l1norm>l1norm_target{
        pn=p.clone();
    // for _ in 0..nt{
        Zip::from(p.slice_mut(s![1..row-1,1..col-1]))
            .and(pn.slice(s![2..row,1..col-1]))
            .and(pn.slice(s![0..row-2,1..col-1]))
            .and(pn.slice(s![1..row-1,2..col]))
            .and(pn.slice(s![1..row-1,0..col-2]))
            .and(b.slice(s![1..row-1,1..col-1]))
            .apply(|pij,&pip1j,&pim1j,&pijp1,&pijm1,&bij|
                *pij=(dy.powi(2)*(pip1j+pim1j)+dx.powi(2)*(pijp1+pijm1)-dx.powi(2)*dy.powi(2)*bij)
                /2./(dx.powi(2)+dy.powi(2)));
        // let mut l1norm_array = Array2::zeros((p.shape()[0],p.shape()[1]));
        // Zip::from(&mut l1norm_array)
        //     .and(&*p)
        //     .and(&pn)
        //     .apply(|vl:f64, vp:f64, vpn:f64|*vl=vp.abs()-vpn.abs());
        // l1norm = l1norm_array.sum();
        Zip::from(p.slice_mut(s![0,..])).apply(|pb|*pb=0.);
        Zip::from(p.slice_mut(s![row-1,..])).apply(|pb|*pb=0.);
        Zip::from(p.slice_mut(s![..,0])).apply(|pb|*pb=0.);
        Zip::from(p.slice_mut(s![..,col-1])).apply(|pb|*pb=0.);
        l1norm= ((abs(&p)-abs(&pn)).sum()/pn.sum()).abs();
        println!("l1norm:{:?}",l1norm);
    }
}

fn main(){
    let nx = 50;
    let ny = 50;
    // let nt = 100;
    let xmin = 0.;
    let xmax = 2.;
    let ymin = 0.;
    let ymax = 1.;
    let dx = (xmax-xmin)/(nx-1) as f64;
    let dy = (ymax-ymin)/(ny-1) as f64;
    
    let mut p=Array2::<f64>::zeros((nx,ny));
    let pd=p.clone();
    let mut b=Array2::<f64>::zeros((nx,ny));
    Zip::from(b.slice_mut(s![(nx/4) as usize,(ny/4) as usize])).apply(|pb|*pb=100.);
    Zip::from(b.slice_mut(s![(3*nx/4) as usize,(3*ny/4) as usize])).apply(|pb|*pb=-100.);
    plot2D(&b, nx, ny, xmax, ymax);
    let row = p.shape()[0];
    let col = p.shape()[1];
    Zip::from(p.slice_mut(s![0,..])).apply(|pb|*pb=0.);
    Zip::from(p.slice_mut(s![row-1,..])).apply(|pb|*pb=0.);
    Zip::from(p.slice_mut(s![..,0])).apply(|pb|*pb=0.);
    Zip::from(p.slice_mut(s![..,col-1])).apply(|pb|*pb=0.);
    plot2D(&p, nx, ny, xmax, ymax);
    poisson2d(&mut p,&b, dx, dy, 1e-6);
    plot2D(&p, nx, ny, xmax, ymax);
}