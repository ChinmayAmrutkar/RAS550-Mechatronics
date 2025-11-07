# EGR 550: Mechatronic Systems (Fall 2025 C)
## Lab 07 – Heater Proportional-Integral-Derivative (PID) Control

### Description
This lab focuses on implementing a standard **Proportional-Integral-Derivative (PID) controller** for the thermal system (Heater 1). The objective is to achieve and maintain a precise **target temperature of 50°C** by continuously adjusting the heater's power (Duty Cycle). The PID block was used to achieve controlled performance, aiming for faster response and reduced steady-state error compared to previous control methods.

---

### Objectives
1. Implement a **PID controller** using the dedicated Simulink block.
2. Achieve and maintain a **desired temperature target of 50°C** with minimal steady-state error and overshoot.
3. Observe and analyze the continuous, proportional control action on the Duty Cycle signal.
4. Visualize the system response by plotting the **Temperature** and **Duty Cycle** over time.

---

### System Setup and Simulink Model

#### Controller Strategy
The control system utilizes a continuous-time PID algorithm to calculate the necessary change in the Duty Cycle based on the error between the setpoint and the measured temperature. The output of the PID block directly determines the heater's duty cycle command.

#### Simulink Model Implementation
The control system was implemented in Simulink on the Arduino Uno hardware board.
<img width="1918" height="986" alt="Model" src="https://github.com/user-attachments/assets/b67ff23e-21b7-4079-9aa2-025d8f1733d2" />

* **Target Setpoint:** The desired temperature is set at **50** degrees Celsius.
* **Temperature Sensing:** The temperature is sensed via the Arduino Analog Input on **Pin A2**. The signal is scaled using blocks with constants $5.0 / 1023.0$ and offset by $100.0$ to output the temperature in °C.
* **Control Loop (Duty Cycle):**
    * The error signal is calculated by subtracting the measured temperature ($\text{temp}$) from the setpoint ($50^{\circ}\text{C}$).
    * This error feeds into the **PID(z)** block (a discrete-time PID implementation).
    * The signal is then passed through a **Saturation** block to respect the physical limits of the PWM, scaled ($255.0 / 100.0$), and converted to an 8-bit integer ($\text{uint8}$) for the PWM command.
    * The final command is sent to the Arduino PWM output on **Pin 5**.
* **Stop Time:** The simulation was set to run for **1500** seconds.

---

### Experimental Results

#### PID Control Performance Plot 
The plot illustrates the heater's temperature response and the corresponding duty cycle adjustments as calculated by the PID controller.

![Plot](https://github.com/user-attachments/assets/17e727cc-1e61-464a-8725-8fead1f1c2db)

**Key Observations from the Plot:**

* **Heating Phase (t < 250s):** The Duty Cycle starts high ($\approx 90\%$) due to the large initial error, quickly driving the temperature up towards the setpoint.
* **Control Action:** The Duty Cycle exhibits **continuous and varying modulation** (ranging between $\approx 50\%$ and $\approx 70\%$) as the temperature approaches the setpoint, characteristic of PID control.
* **Target Acquisition (t $\approx$ 250s):** The temperature successfully reaches and stabilizes around the **$50^{\circ}\text{C}$ target**.
* **Steady State (t > 250s):** The temperature is maintained near $50^{\circ}\text{C}$ with relatively small fluctuations, indicating that the controller is effectively holding the setpoint.

---

### Key Learning Outcomes
* Successfully implemented a **PID controller** for temperature regulation.
* Understood the characteristics of **continuous PID control** and its advantages over Bang-Bang control.
* Observed how the PID controller's action results in a **stable temperature response** with minimal overshoot and steady-state error.
* Demonstrated the use of the **PID Controller block** and proper signal scaling for hardware output.
---
