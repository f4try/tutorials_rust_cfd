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

fn laplace2d(p:&mut Array2<f64>,dx:f64,dy:f64,l1norm_target:f64){
    let mut l1norm =1.;
    let row = p.shape()[0];
    let col = p.shape()[1];
    // let mut pn = Array2::zeros((p.shape()[0],p.shape()[1]));
    let mut pn = p.clone();
    while l1norm>l1norm_target{
        Zip::from(p.slice_mut(s![1..row-1,1..col-1]))
            .and(pn.slice(s![2..row,1..col-1]))
            .and(pn.slice(s![0..row-2,1..col-1]))
            .and(pn.slice(s![1..row-1,2..col]))
            .and(pn.slice(s![1..row-1,0..col-2]))
            .apply(|pij,&pip1j,&pim1j,&pijp1,&pijm1|
                *pij=(dy.powi(2)*(pip1j+pim1j)+dx.powi(2)*(pijp1+pijm1))/2./(dx.powi(2)+dy.powi(2)));
        // let mut l1norm_array = Array2::zeros((p.shape()[0],p.shape()[1]));
        // Zip::from(&mut l1norm_array)
        //     .and(&*p)
        //     .and(&pn)
        //     .apply(|vl:f64, vp:f64, vpn:f64|*vl=vp.abs()-vpn.abs());
        // l1norm = l1norm_array.sum();
        pn=p.clone();
        Zip::from(p.slice_mut(s![0,..])).apply(|pb|*pb=0.);
        Zip::from(p.slice_mut(s![row-1,..]))
            .and(&Array1::linspace(0.,2.,col))
            .apply(|pb,&b|*pb=b);
        Zip::from(p.slice_mut(s![..,0]))
            .and(pn.slice(s![..,1]))
            .apply(|pb,&b|*pb=b);
        Zip::from(p.slice_mut(s![..,col-1]))
            .and(pn.slice(s![..,col-2]))
            .apply(|pb,&b|*pb=b);
        l1norm= (abs(&p)-abs(&pn)).sum()/pn.sum();
    }
}

fn main(){
    let nx = 31;
    let ny = 31;
    let c = 1;
    let dx = 2./(nx-1) as f64;
    let dy = 2./(ny-1) as f64;
    
    let mut p=Array2::<f64>::zeros((nx,ny));
    let pn=p.clone();
    let row = p.shape()[0];
    let col = p.shape()[1];
    Zip::from(p.slice_mut(s![0,..])).apply(|pb|*pb=0.);
    Zip::from(p.slice_mut(s![row-1,..]))
        .and(&Array1::linspace(0.,2.,col))
        .apply(|pb,&b|*pb=b);
    Zip::from(p.slice_mut(s![..,0]))
        .and(pn.slice(s![..,1]))
        .apply(|pb,&b|*pb=b);
    Zip::from(p.slice_mut(s![..,col-1]))
        .and(pn.slice(s![..,col-2]))
        .apply(|pb,&b|*pb=b);
    plot2D(&p, nx, ny, 2., 2.);
    laplace2d(&mut p, dx, dy, 1e-4);
    plot2D(&p, nx, ny, 2., 2.);
}