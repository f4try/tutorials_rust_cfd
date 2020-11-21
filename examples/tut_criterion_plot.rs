use itertools_num::linspace;
use criterion_plot::prelude::*;

fn main()->Result<Child,Box<(dyn std::error::Error + 'static)>>{
    let ref xs = linspace::<f64>(-10., 10., 51).collect::<Vec<_>>();

    Figure::new()
        .configure(Key, |k| {
            k.set(Boxed::Yes)
            .set(Position::Inside(Vertical::Top, Horizontal::Left))
        })
        .plot(LinesPoints {
                x: xs,
                y: xs.iter().map(|x| x.sin()),
            },
            |lp| {
                lp.set(Color::DarkViolet)
                    .set(Label("sin(x)"))
                    .set(LineType::Dash)
                    .set(PointSize(1.5))
                    .set(PointType::Circle)
            })
        .plot(Steps {
                x: xs,
                y: xs.iter().map(|x| x.atan()),
            },
            |s| {
                s.set(Color::Rgb(0, 158, 115))
                .set(Label("atan(x)"))
                .set(LineWidth(2.))
            })
        .plot(Impulses {
                x: xs,
                y: xs.iter().map(|x| x.atan().cos()),
            },
            |i| {
                i.set(Color::Rgb(86, 180, 233))
                .set(Label("cos(atan(x))"))
            })
        .draw()  // (rest of the chain has been omitted)
}