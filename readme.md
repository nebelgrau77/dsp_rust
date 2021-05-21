# Digital Signal Processing in Rust

This library contains some DSP functions, translated into Rust from the [DSP from the ground up in C](https://www.udemy.com/course/digital-signal-processing-dsp-from-ground-uptm-in-c/) course on Udemy that I'm currently taking.

I'm trying to make them as idiomatic as possible, not a simple translation from C, which could result in lesser performance. The library is intended to be `no_std` ready, hence the use of `micromath` crate.
