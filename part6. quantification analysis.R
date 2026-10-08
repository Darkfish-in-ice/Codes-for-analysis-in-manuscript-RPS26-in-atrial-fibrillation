### qPCR data analysis
control_data <- scan(text = "1.0436796	0.9160781	1.0506511	1.0033555	1.0140665	0.9721693")
angii_data   <- scan(text = "0.424947	0.1216392	0.1922916	0.2621534	0.2948628	0.3741837")
t.test(control_data, angii_data, paired = FALSE)

### Western blot analysis
control_data <- scan(text = "0.9957523	1.0161977	0.9691431	1.0246147	0.9333289	1.0609632")
angii_data   <- scan(text = "0.585819	0.6461114	0.4338944	0.3542234	0.7160133	0.5497199")
t.test(control_data, angii_data, paired = FALSE)