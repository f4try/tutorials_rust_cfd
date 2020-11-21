use ndarray::{Array,arr1,s};
// use gnuplot::*;
macro_rules! printVar{
    ($name:expr)=>{
        println!("{:?}",$name);
    };
}

fn main(){
    let myarray = Array::linspace(0.,5.,10);
    printVar!(myarray);
    let myarr = arr1(&[0.,1.,2.,3.]);
    printVar!(myarr);
    let a = 5;
    let b = "five";
    let c = 5.0;
    printVar!(a);
    printVar!(b);
    printVar!(c);
    for i in 0..5{
        println!("Hi");
    }
    for i in 0..3{
        for j in 0..3{
            println!("{},{}",i,j);
        }
        println!("----------");
    }
    let myvals = Array::from(vec![1.,2.,3.,4.,5.]);
    printVar!(myvals);
    printVar!(myvals[0]);
    printVar!(myvals[4]);
    printVar!(myvals.slice(s![0..3]));

    let a = Array::range(1.,6.,1.);
    printVar!(a);
    let mut b = a;
    // a[2]=18.;
    b[2]=18.;
    // printVar!(a);
    printVar!(b);

    
}