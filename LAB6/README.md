# EGR 550: Mechatronic Systems (Fall 2025 C)
## Lab 06 – Heater Bang-Bang (On/Off) Control

### Description
This lab focuses on implementing a simple **Bang-Bang (On/Off) controller** for the thermal system (Heater 1). The objective is to maintain a target temperature of **50°C** by switching the heater's power (Duty Cycle) on and off based on the temperature error. A **dead-band** is introduced into the control logic to prevent rapid switching and system chatter, a common issue in Bang-Bang control.

---

### Objectives
1. Implement a **Bang-Bang (On/Off) controller** using a MATLAB Function block in Simulink.
2. Achieve and maintain a **desired temperature target of 50°C**.
3. Incorporate a **dead-band** to stabilize the control signal and heater temperature.
4. Visualize the system response by plotting the **Temperature** and **Duty Cycle** over time.

---

### System Setup and Simulink Model

#### Controller Strategy
The controller logic is based on the temperature error ($u$), where $u = \text{Target Temperature} - \text{Measured Temperature}$. The target temperature is set to **50.0°C**. The MATLAB Function block implements a multi-level dead-band for the control action. The base duty cycle for the heater is 50.0% (representing the 'off' state), and power is added on top of this base.

<img width="408" height="362" alt="Bang_bang_function_block" src="https://github.com/user-attachments/assets/38d2a365-1982-43e8-be36-2c96781fea87" />


**Control Law Implemented in MATLAB Function:** 
| Error ($u = 50 - \text{Temp}$) | Control Action ($y$) | Effective Duty Cycle Command | Description |
|:------------------------------:|:--------------------:|:------------------------------:|:------------|
| $u \geq 20.0$ | $y = 40.0 + 50.0$ | $90.0\%$ | **High Power:** Temperature is far below target. |
| $10.0 \leq u < 20.0$ | $y = 20.0 + 50.0$ | $70.0\%$ | **Medium Power:** Temperature is moderately below target. |
| $0.0 \leq u < 10.0$ | $y = 10.0 + 50.0$ | $60.0\%$ | **Low Power:** Temperature is slightly below target. |
| $u < 0.0$ (Temp $> 50$) | $y = 0.0 + 50.0$ | $50.0\%$ | **Heater Off/Base:** Temperature is above target (cooling). |

*Note: The effective control logic uses the difference between the target (50) and the measured temperature as input ($u$) to the function block, as shown in the Simulink model's Duty Cycle subsystem.*

#### Simulink Model Implementation
The control system was implemented in Simulink using the Arduino Uno hardware board.
<img width="1918" height="987" alt="Model" src="https://github.com/user-attachments/assets/cf951e30-7227-4d80-97ae-8baaf03a2fca" />

* **Temperature Sensing:** The Arduino Analog Input **Pin A2** is used for temperature measurement. The signal is scaled and offset (100.0) to output the temperature in °C.
* **Clock:** A simple clock is used to track the experiment time.
* **Control Loop (Duty Cycle):**
    * The measured temperature is subtracted from the desired setpoint ($50.0^{\circ}\text{C}$) to generate the error signal, which feeds into the **fcn** (MATLAB Function) block.
    * The output of the Bang-Bang function is passed through a **Saturation** block and then scaled to an 8-bit integer (0-255) for the PWM command.
    * The final command is sent to the Arduino PWM output on **Pin 5**.
* **Stop Time:** The simulation was run for **1500** seconds (though the plot shows 500 seconds).

---

### Experimental Results

#### Bang-Bang Control Performance Plot
The experiment successfully demonstrated Bang-Bang control, where the duty cycle switches in response to the error to drive the temperature towards the target .
![Plot](https://github.com/user-attachments/assets/ef0cd4a6-9edc-4346-a1a9-1d67fdd6b382)

**Key Observations from the Plot:**

* **Heating Phase (t < 150s):** The duty cycle starts at the maximum (90%) because the error is large. As the temperature approaches the target, the control logic steps down the power to 70% and then 60%.
* **Target Acquisition (t ≈ 150s):** Once the temperature reaches the target of $50^{\circ}\text{C}$, the error becomes negative, and the Duty Cycle drops to the base level ($50\%$), indicating the heater is effectively turned **off** (or at its minimum power).
* **Steady State (t > 150s):** The temperature successfully settles around the **$50^{\circ}\text{C}$ target**, exhibiting the characteristic oscillation (limit cycle) inherent to Bang-Bang control with a dead-band.
* **Duty Cycle Response:** The control signal shows clear step-down and step-up behavior, demonstrating the multi-level dead-band logic working to modulate the power before the final on/off state is reached.

---

### Key Learning Outcomes
* Successfully implemented a **Bang-Bang (On/Off) controller** for a thermal system.
* Understood how to use the **MATLAB Function block** to implement custom, conditional control logic in Simulink.
* Observed the effect of a **dead-band** in reducing control chatter and stabilizing the temperature around the setpoint.
* Analyzed the **limit cycle oscillation** that is characteristic of this type of discontinuous control.
