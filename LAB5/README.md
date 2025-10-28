# EGR 550: Mechatronic Systems (Fall 2025 C)
## Lab 05 – Heater Control and System Identification using Step Response

### Description
This lab focuses on thermal system modeling and control through experimental data collection and step response analysis. Heater 1 is controlled using a Pulse Width Modulated (PWM) signal, and its temperature response is recorded to identify system dynamics. The objective is to determine the heater’s time constant, steady-state gain, and overall behavior using MATLAB or Excel-based modeling.

---

### Objectives
1. Control Heater 1 using a PWM step input in Simulink.  
2. Observe and record temperature response for each input step.  
3. Determine steady-state temperature and time constant.  
4. Collect time, input, and output data in real time.  
5. Develop mathematical models representing the thermal response of the system.

---

### Experimental Procedure

#### Step 1: System Setup
- Connect Heater 1 to the control interface via PWM output.
- Configure Simulink for real-time monitoring using the **Monitor and Tune** function.
- Sampling interval: 0.5 seconds.
- Expected total data points:  
  2 × 500 seconds × 3 steps = 3000 samples.

---

#### Step 2: Step Input Sequence
| Step | PWM Command (%) | Duration (s) | Description |
|------|-----------------|---------------|--------------|
| 1 | 0 → 80 | 500 | Turn on Heater 1 and allow temperature to reach steady state. |
| 2 | 80 → 40 | 500 | Reduce heater power to observe controlled cooling. |
| 3 | 40 → 0 | 500 | Turn off the heater to observe natural cooling behavior. |

*Note:* The easiest way to modify PWM input values during the experiment is to use the **Monitor and Tune** feature in Simulink.

---

#### Step 3: Data Collection
- Record:
  - Time (s)
  - PWM Input (%)
  - Temperature Output (°C)
- Ensure each step reaches steady state before applying the next step input.
- Export all data to MATLAB or Excel for further analysis.

---

### Data Analysis and Modeling

#### Step 3 (0 → 80%)
Model the heating phase using a first-order transfer function of the form:

\[
G(s) = \frac{K}{\tau s + 1}
\]

Where:  
- \( K \) = Steady-state gain  
- \( \tau \) = Time constant

Identify the rise time, settling time, and steady-state temperature.

#### Step 4 (80 → 40%)
Analyze the partial cooling response and compare the derived model parameters to those obtained in Step 3.

#### Step 5 (40 → 0%)
Model the complete cooling response and determine the natural cooling rate. Compare the model fit against the measured data.

---

### Expected Results
- Temperature vs. Time plots for all step changes.
- Comparison of PWM input signal and temperature output.
- System parameters (K, τ) determined for each step.
- Overlay plots showing measured vs. modeled responses.


| Plot No. | Step Description | Input PWM Range (%) | Duration (s) | Plot Type | Purpose |
|:---------:|:----------------|:-------------------:|:-------------:|:-----------|:----------|
| 1 | Step 3: Heater ON (0 → 80%) | 0 to 80 | 500 | Temperature vs. Time | To analyze the heating response and determine steady-state temperature and time constant. |
| 2 | Step 4: Heater Partial Cooling (80 → 40%) | 80 to 40 | 500 | Temperature vs. Time | To model the cooling response and compare it with the heating dynamics. |
| 3 | Step 5: Heater OFF (40 → 0%) | 40 to 0 | 500 | Temperature vs. Time | To model the natural cooling rate and validate system response when heater is turned off. |
| 4 | Combined Input-Output Plot | 0 → 80 → 40 → 0 | 1500 | PWM Input and Temperature Output vs. Time | To show the full experimental cycle of input steps and corresponding thermal response. |


---

### Key Learning Outcomes
- Understand the dynamic behavior of a thermal system under step inputs.
- Implement PWM-based actuation for precise temperature control.
- Apply system identification methods to real-world experimental data.
- Utilize Simulink’s Monitor and Tune capability for real-time control and data acquisition.

---

### References
- MATLAB & Simulink Documentation – System Identification Toolbox  
- Dorf, R.C., & Bishop, R.H., *Modern Control Systems*  
- Course Lecture Notes on Thermal System Dynamics

---

### Recommended Repository Structure
