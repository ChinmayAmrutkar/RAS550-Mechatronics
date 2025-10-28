# EGR 550: Mechatronic Systems (Fall 2025 C)
## Lab 05 – Heater Control and System Identification using Step Response

### Description
This lab focuses on **thermal system modeling and control** through experimental data collection and step response analysis. Heater 1 is controlled using a Pulse Width Modulated (PWM) signal, and its temperature response is recorded to identify system dynamics. The objective is to determine the heater’s **time constant ($\tau$)**, **steady-state gain ($G$ or $K$)**, and overall behavior using MATLAB-based nonlinear least-squares optimization modeling.

---

### Objectives
1. Control Heater 1 using a PWM step input in Simulink. 
2. Observe and record temperature response for each input step. 
3. Determine **steady-state temperature** and **time constant** ($\tau$) for each segment. 
4. Collect time, input, and output data in real time. 
5. Develop and fit a **first-order plus dead time mathematical model** to the thermal response of the system using MATLAB's `lsqcurvefit` function.

---

### Experimental Procedure

#### Step 1: System Setup
* Connect Heater 1 to the control interface via PWM output.
* Configure Simulink for real-time monitoring using the **Monitor and Tune** function.
* Sampling interval: **0.5 seconds**.
* Total data points collected: $2 \times 500 \text{ seconds} \times 3 \text{ steps} / (0.5 \text{ s/sample}) = **3000$ samples** for temperature and duty cycle.

---

#### Step 2: Step Input Sequence
The experiment was divided into three 500-second segments, with the PWM command adjusted using the **Monitor and Tune** feature in Simulink.

| Step | Segment | PWM Command (%) | Duration (s) | Description |
|:-----|:--------|:----------------|:--------------|:-------------|
| 3 | Seg1 | 0 → 80 | 500 | **Heating Phase**: Turn on Heater 1 and allow temperature to reach steady state. |
| 4 | Seg2 | 80 → 40 | 500 | **Partial Cooling**: Reduce heater power to observe controlled cooling. |
| 5 | Seg3 | 40 → 0 | 500 | **Natural Cooling**: Turn off the heater to observe natural cooling behavior. |

---

### Data Analysis and Modeling

#### Step 6-8: MATLAB Model Fitting Results
Each segment's response was modeled using a **first-order plus dead time (FOPDT)** model:
$$y(t) = y_{\text{init}} + G \cdot \text{step}(t - t_0 - \theta_p) \cdot \left(1 - e^{-(t - (t_0 + \theta_p))/\tau}\right)$$
where $t$ is the local time $t_l$, and $t_0$ is 10 seconds.
The four parameters — $y_{\text{init}}$ (**Initial Temperature**), $G$ (**Steady-state Gain**), $\tau$ (**Time Constant**), and $\theta_p$ (**Dead Time**) — were estimated by minimizing the sum of squared differences between the measured data $y(t)$ and the model prediction using MATLAB's nonlinear least-squares optimization function `lsqcurvefit`.

**Estimated Parameters (y\_init, G, tau, theta\_p)**:

| Segment | y\_init (°C) | G (°C/%) | $\tau$ (s) | $\theta_p$ (s) | Step Description |
|:-------:|:--------------|:---------|:----------|:---------------|:-----------------|
| **Seg1** (0 → 80%) | 24.974 | 57.304 | 160.91 | 21.27 | **Heating Phase** (Step 3) |
| **Seg2** (80 → 40%) | 78.001 | -23.799 | 132.2 | 17.486 | **Partial Cooling** (Step 4) |
| **Seg3** (40 → 0%) | **54.412** | **-29.623** | **161.64** | **23.25** | **Natural Cooling** (Step 5) |

**Key Observations from Model Parameters:**

* **Seg1 ($0 \rightarrow 80\%$):** The steady-state gain $G \approx 57.304 \text{ °C}/\text{PWM %}$ is positive, indicating heating, with a time constant $\tau \approx 160.91 \text{ s}$.
* **Seg2 ($80 \rightarrow 40\%$):** The gain $G \approx -23.799 \text{ °C}/\text{PWM %}$ is negative, representing the cooling response, with a faster time constant $\tau \approx 132.2 \text{ s}$.
* **Seg3 ($40 \rightarrow 0\%$):** This segment represents the final natural cooling, having the largest dead time $\theta_p \approx 71.646 \text{ s}$ and a time constant $\tau \approx 162.32 \text{ s}$.

---

### Expected Results and Deliverables

* **Temperature vs. Time plots with PWM input signal** for all step changes (Segments 1, 2, and 3).
<img width="1600" height="850" alt="image" src="https://github.com/user-attachments/assets/307de784-6989-4a98-8a81-5c75240a7025" />

* **Overlay plots** for each segment showing the measured temperature data vs. the FOPDT modeled responses.
<img width="856" height="616" alt="image" src="https://github.com/user-attachments/assets/7f249a7f-bf43-4e2d-8328-855a20e9a6c3" />

* **System parameters** ($y_{\text{init}}$, $G$, $\tau$, $\theta_p$) determined and tabulated for each step.

| Plot No. | Step Description | Input PWM Range (%) | Purpose |
|:---------:|:----------------|:-------------------:|:----------|
| 1 | Step 3: Heater ON (0 → 80%) | 0 to 80 | To analyze the heating response, determine steady-state temperature, gain ($G$), and time constant ($\tau$). |
| 2 | Step 4: Heater Partial Cooling (80 → 40%) | 80 to 40 | To model the controlled cooling response and compare parameters with the heating dynamics. |
| 3 | Step 5: Heater OFF (40 → 0%) | 40 to 0 | To model the natural cooling rate and validate the system response with the heater completely off. |
| 4 | Combined Input-Output Plot | 0 → 80 → 40 → 0 | To show the full experimental cycle of input steps and corresponding thermal response. |

---

### Key Learning Outcomes
* Understand the **dynamic behavior of a thermal system** under step inputs.
* Implement **PWM-based actuation** for precise temperature control.
* Apply **system identification methods** (FOPDT model fitting) to real-world experimental data.
* Utilize Simulink’s **Monitor and Tune** capability for real-time control and data acquisition.
* Use MATLAB's **nonlinear least-squares optimization** to estimate model parameters.

---

### References
* MATLAB & Simulink Documentation – System Identification Toolbox 
* Dorf, R.C., & Bishop, R.H., *Modern Control Systems* * Course Lecture Notes on Thermal System Dynamics
