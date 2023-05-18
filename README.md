# ME599_Torque_Estimation
Utilize SINDy to attempt to extrapolate torque from the time derivative approximation of idq. The goal of this project is to identify the governing ODEs relating the electrical current, voltage, and angle to the motor output torque. To do this, we will leverage SINDy in an attempt to find a sparse system of equations that models the time series data and can be used to predict output torques in future time-steps. 

Citations:
[1] Rick Chartrand, "Numerical differentiation of noisy,
nonsmooth data," ISRN Applied Mathematics, Vol. 2011, Article ID 164564, 
2011. 

[2] https://www.kaggle.com/datasets/graxlmaxl/identifying-the-physics-behind-an-electric-motor?select=Time_series_2000rpm.csv
