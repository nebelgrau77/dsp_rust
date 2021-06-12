use micromath::F32Ext;
const PI: f32 = 3.14159265358979; 
use core::ops::Div;

/// Calculate Discrete Fourier Transform of the signal   - WORKS
pub fn calc_signal_dft(input_signal: &[f32], output_signal_rex: &mut [f32], output_signal_imx: &mut [f32]) {
        
    for (k, out_re) in output_signal_rex.iter_mut().enumerate() {
        for (i, input) in input_signal.iter().enumerate() {

            *out_re += *input*(2.0*PI*(k as f32)*(i as f32)/(input_signal.len() as f32)).cos()

        }
    }

    for (k, out_im) in output_signal_imx.iter_mut().enumerate() {
        for (i, input) in input_signal.iter().enumerate() {

            *out_im -= *input*(2.0*PI*(k as f32)*(i as f32)/(input_signal.len() as f32)).sin()

        }
    }
    
}

/// Calculate Inverse Discrete Fourier Transform of the signal
pub fn calc_signal_idft(input_rex: &mut [f32], input_imx: &mut [f32], output_idft: &mut [f32]) {

    let input_len = input_rex.len();
    let output_len = output_idft.len();

    for k in 0..input_len {
        input_rex[k] = input_rex[k]/(input_len as f32);
        input_imx[k] = -input_imx[k]/(input_len as f32);
        }
    
    input_rex[0] = input_rex[0]/2.0;
    input_imx[0] = -input_imx[0]/2.0;
    
    for k in 0..input_len {
        for i in 0..output_len {
            output_idft[i] += input_rex[k] * ((2.0*PI*(i as f32)*(k as f32))/(output_len as f32)).cos();
            output_idft[i] += input_imx[k] * ((2.0*PI*(i as f32)*(k as f32))/(output_len as f32)).sin();
        }
    }
}


/// TO DO
/// Calculate Inverse Discrete Fourier Transform of the signal
pub fn calc_signal_idft_better(input_rex: &mut [f32], input_imx: &mut [f32], output_idft: &mut [f32]) {

    let input_len = input_rex.len();
    let output_len = output_idft.len();


    for (re, im) in input_rex.iter_mut().zip(input_imx.iter_mut()) {
        *re = *re/input_len as f32;
        *im = *im*(-1.0)/input_len as f32;
    }

    /*
    for k in 0..input_len {
        input_rex[k] = input_rex[k]/(input_len as f32);
        input_imx[k] = -input_imx[k]/(input_len as f32);
        }
     */
    input_rex[0] = input_rex[0]/2.0;
    input_imx[0] = -input_imx[0]/2.0;
    
    for k in 0..input_len {
        for i in 0..output_len {
            output_idft[i] += input_rex[k] * ((2.0*PI*(i as f32)*(k as f32))/(output_len as f32)).cos();
            output_idft[i] += input_imx[k] * ((2.0*PI*(i as f32)*(k as f32))/(output_len as f32)).sin();
        }
    }
}




/// Calculate magnitude of the signal  - WORKS
pub fn calc_signal_magnitude(input_rex: &[f32], input_imx: &[f32], output_mag: &mut [f32]) {
    
    for (rex, imx, mag) in input_rex.iter()
                        .zip(input_imx.iter())
                        .zip(output_mag.iter_mut())
                        .map(|((rex,imx),mag)| (rex,imx,mag)) 
        {                           
            *mag = (rex.powi(2) + imx.powi(2)).powf(0.5);        
        }
        
}


/// A more idiomatic version of the signal mean calculation - WORKS
pub fn calc_signal_mean(signal: &[f32]) -> f32 {

    let mean: f32 = signal.iter()
    .fold(0_f32, |sum, x| sum + x)
    .div(signal.len() as f32);

    mean   
    
}

/// A more idiomatic version of the signal variance calculation - WORKS
pub fn calc_signal_variance(signal: &[f32], signal_mean: f32) -> f32 {
    
    let variance: f32 = signal.iter()
    .fold(0_f32, |sum, x| sum + (x-signal_mean).powi(2))
    .div((signal.len()-1) as f32);
    
    variance
    
}

/// Calculate standard deviation of the signal
pub fn calc_signal_stddev(signal_variance: f32) -> f32 {
    signal_variance.powf(0.5)
}


/// Calculate running sum of the signal
pub fn calc_running_sum(input: &[f32], output: &mut [f32]) {
    
    let mut acc = 0_f32;
    
    for (i,j) in input.iter().zip(output.iter_mut()) {
        *j = *i + acc;
        acc += *i;
    }
    
}


/// Convert complex signal values to polar coordinates (magnitude and phase) - WORKS
pub fn rect_to_polar(input_rex: &mut [f32], 
    input_imx: &[f32],
    output_mag: &mut [f32],
    output_phase: &mut [f32]) 
    {

        let sig_len = input_rex.len();

        for (re,im,mag,ph) in input_rex.iter_mut()
                        .zip(input_imx.iter())
                        .zip(output_mag.iter_mut())
                        .zip(output_phase.iter_mut())
                        .map(|(((re, im), mag), ph)| (re,im,mag,ph)) 
                        
            {
        
                *mag = (re.powi(2) + im.powi(2)).powf(0.5);

                if *re == 0.0 {
                    *re = 10.0.powi(-20);
                    *ph = ((*im)/(*re)).atan();
                }

                if *re < 0.0 && *im < 0.0 {
                    *ph = *ph - PI;
                }

                if *re < 0.0 && *im >= 0.0 {
                    *ph = *ph + PI;
                }

            }

    }


/// Calculate complex DFT of the signal (conversion from time to frequency domain) - WORKS
pub fn complex_dft(time_rex: &[f32],
                        time_imx: &[f32],
                        freq_rex: &mut [f32],
                        freq_imx: &mut [f32])

    {

        let sig_len = time_rex.len();

        for (k, f_re, f_im) in freq_rex.iter_mut()
                            .zip(freq_imx.iter_mut())
                            .enumerate()
                            .map(|(k,(f_re,f_im))| (k,f_re,f_im)) {
            for (i, t_re, t_im) in time_rex.iter()
                                .zip(time_imx.iter())
                                .enumerate()
                                .map(|(i,(t_re,t_im))| (i,t_re,t_im)) {

                                let SR = ((2.0*PI*(k as f32)*(i as f32))/(sig_len as f32)).cos();
                                let SI = -((2.0*PI*(k as f32)*(i as f32))/(sig_len as f32)).sin();

                                *f_re = *f_re + *t_re * SR - *t_im * SI;
                                *f_im = *f_im + *t_im * SI - *t_im * SR;

                            }

                        }
        }


/// convolution of two signals (signal and kernel)
pub fn convolution(input_signal: &[f32], impulse: &[f32], output: &mut [f32]) {
    for i in 0..output.len() {
        //output[i] = 3.0;
        for j in 0..input_signal.len() {
            for k in 0..impulse.len() {
                output[j+k] += input_signal[j] * impulse[k]
            }
        }
    }    
}




#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}


//use std::f32::consts::PI;

/*
/// Calculate Discrete Fourier Transform of the signal
pub fn calc_signal_dft(input: &[f32], output_rex: &mut [f32], output_imx: &mut [f32]) {
        
    for k in 0..output_rex.len() {    
        for i in 0..input.len() {            
            output_rex[k] += input[i]*(2.0*PI*(k as f32)*(i as f32)/(input.len() as f32)).cos()
        }         
    }

    for k in 0..output_imx.len() {    
        for i in 0..input.len() {            
            output_imx[k] -= input[i]*(2.0*PI*(k as f32)*(i as f32)/(input.len() as f32)).sin()
        }         
    }

}
 */

 /*
/// Calculate magnitude of the signal 
pub fn calc_signal_magnitude(input_rex: &mut [f32], input_imx: &mut [f32], output_mag: &mut [f32]) {
    for i in 0..output_mag.len() {
        output_mag[i] = (input_rex[i].powi(2) + input_imx[i].powi(2)).powf(0.5);
    }
}
 */

 /*
/// Calculate mean value of the signal
pub fn calc_signal_mean(signal: &[f32]) -> f32 {
    let mut _mean: f32 = 0.0;
    for i in signal {
        _mean += i;
    }
    _mean / (signal.len() as f32)
    
}
 */

 /*
/// Calculate variance of the signal
pub fn calc_signal_variance(signal: &[f32], signal_mean: f32) -> f32 {
    let mut _variance: f32 = 0.0;
    for i in signal {
        _variance += (i - signal_mean).powi(2); //add squared difference
    }
    _variance / ((signal.len() - 1) as f32)
    
}

 */

/*
/// Calculate running sum of the signal
pub fn calc_running_sum(input: &[f32], output: &mut [f32]) {
    for i in 0..output.len() {    
        if i == 0 {
            output[i] = input[i];
        } else {
            output[i] += output[i-1] + input[i]
        }         
    }
}
 */




/*
/// Convert complex signal values to polar coordinates (magnitude and phase)
pub fn rect_to_polar(input_rex: &mut [f32], 
                    input_imx: &[f32],
                    output_mag: &mut [f32],
                    output_phase: &mut [f32]) 
    {
        
        let sig_len = input_rex.len();

        
        for k in 0..sig_len {
            output_mag[k] = (input_rex[k].powi(2) + input_imx[k].powi(2)).powf(0.5);
        

            if input_rex[k] == 0.0 {
                input_rex[k] = 10.0.powi(-20);
                output_phase[k] = (input_imx[k]/input_rex[k]).atan();
            }

            if input_rex[k]<0.0 && input_imx[k] < 0.0 {
                output_phase[k] = output_phase[k] - PI;
            }

            if input_rex[k]<0.0 && input_imx[k] >= 0.0 {
                output_phase[k] = output_phase[k] + PI;
            }

        }

    }

 */
/*
/// Calculate complex DFT of the signal (conversion from time to frequency domain)
pub fn complex_dft(time_rex: &[f32],
                time_imx: &[f32],
                freq_rex: &mut [f32],
                freq_imx: &mut [f32])
    
    {
        let sig_len = time_rex.len();

        for k in 0..sig_len-1 {
            for i in 0..sig_len-1 {
                let SR = ((2.0*PI*(k as f32)*(i as f32))/(sig_len as f32)).cos();
                let SI = -((2.0*PI*(k as f32)*(i as f32))/(sig_len as f32)).sin();

                freq_rex[k] = freq_rex[k] + time_rex[i] * SR - time_imx[i] * SI;
                freq_imx[k] = freq_imx[k] + time_imx[i] * SI - time_imx[i] * SR;

            }
        }
    }

 */
