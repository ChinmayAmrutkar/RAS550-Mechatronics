# EGR 550: Mechatronic Systems (Fall 2025 C)
## Lab 04 – Closed-Loop Position Control using a PD Controller

### Description
This lab advances our control systems exploration by focusing on the **closed-loop position control** of a DC motor. We will implement a **PD (Proportional–Derivative) controller** in **MATLAB Simulink** to precisely regulate the angular position of the motor's output shaft.

The primary objective is to design and tune a PD controller that commands the motor shaft to smoothly follow a **sine wave reference signal**. This involves translating the desired rotation into a target encoder count and using continuous feedback to minimize position error and test the system's dynamic tracking capabilities. ⚙️

---

### Hardware Diagrams
The hardware for this lab remains the same as in Lab 03.

| Arduino Uno | L298N Motor Driver | GA12-N20 Motor with Encoder |
|:-----------:|:------------------:|:---------------------------:|
|  |  |  |
| [Arduino Uno Datasheet (PDF)](https://docs.arduino.cc/resources/datasheets/A000066-datasheet.pdf) | [L298N Motor Driver Datasheet (PDF)](https://www.handsontec.com/dataspecs/L298N%20Motor%20Driver.pdf) | [GA12-N20 Motor Datasheet (PDF)](https://www.handsontec.com/dataspecs/GA12-N20.pdf) |

### Circuit Diagram
The circuit wiring is identical to the previous lab.



---

### Key Calculations 📐
To control the output shaft's position, we first need to determine the target encoder count corresponding to 4 revolutions.

* **Counts Per Revolution (Motor):** 7 counts/revolution
* **Encoder Type:** Quadrature, which provides 4 signal edges per pulse.
* **Effective Resolution:** $7 \times 4 = 28$ counts per motor revolution.
* **Gear Ratio (G):** 50:1
* **Target Output Shaft Revolutions:** 4

First, we calculate the required motor revolutions:
$$\text{Motor Revolutions} = \text{Output Shaft Revolutions} \times G = 4 \times 50 = \textbf{200 revolutions}$$

Next, we convert the motor revolutions into the total encoder count for the full range of motion:
$$\text{Target Encoder Count} = \text{Motor Revolutions} \times \text{Effective Resolution} = 200 \times 28 = \textbf{5600 counts}$$

This value is used to normalize the position feedback, so a command signal with an amplitude of 1.0 corresponds to the full 4-revolution travel distance.

---

### Tasks Completed
1.  **PD Controller Development**
    * Implemented a Proportional–Derivative (PD) controller in Simulink.
    * Tuned the proportional ($K_p$) and derivative ($K_d$) gains to achieve responsive and stable position tracking with minimal error.

2.  **Encoder Integration**
    * Successfully read quadrature encoder signals from channels A and B.
    * Utilized the Simulink **Encoder block** to decode the real-time angular position of the motor shaft.

3.  **Position Control**
    * Generated a **discrete sine wave** as a reference signal to test continuous tracking.
    * Commanded the motor's output shaft to follow the smooth, oscillating motion defined by the reference signal.
    * Validated the controller's ability to dynamically track the changing position setpoint.
  
---

### Submissions
* **Simulink Model Screenshot:** Block diagram showing the PD controller, encoder interface, normalized feedback, and motor driver logic with a deadband.

<img width="1657" height="562" alt="image" src="https://github.com/user-attachments/assets/eb4db903-f67f-4df0-99c5-a263a2f58d05" />


* **Results:**

* **Video Demonstration (Sine Wave Position Control)**


https://github.com/user-attachments/assets/daa6a3b6-4544-4183-8c5c-330e5cb4c047



* **Position Plot (Sine Wave Position Control)**
<img width="689" height="516" alt="image" src="https://github.com/user-attachments/assets/eaad2439-bfca-4f99-a559-b447e8e95eac" />
