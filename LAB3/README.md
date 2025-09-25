# EGR 550: Mechatronic Systems (Fall 2025 C)  
## Lab 03 – Closed-Loop Velocity Control using a PI Controller  

### Description  
This lab focuses on **closed-loop velocity control** of a DC motor using a **PI (Proportional–Integral) controller** implemented in **MATLAB Simulink**. Unlike the previous lab’s open-loop PWM approach, here we incorporate encoder feedback for velocity estimation and regulate the motor speed using control inputs.  

The objective is to design and test a PI controller to regulate velocity in response to both **step inputs** and a **sine wave reference signal**.  

---

### Hardware Diagrams  

| Arduino Uno | L298N Motor Driver | GA12-N20 Motor with Encoder |
|-------------|--------------------|-----------------------------|
| <img width="250" alt="Arduino Uno" src="https://github.com/user-attachments/assets/2203837f-746f-4b6e-8214-7adb68a8c7d0" /> | <img width="250" alt="L298N Motor Driver" src="https://github.com/user-attachments/assets/048af7fa-cc8a-4546-8d2d-b68c417c9854" /> | <img width="250" alt="GA12-N20 Motor with Encoder" src="https://github.com/user-attachments/assets/df88ccd8-d7fd-4a69-9c8e-56ec57d347f2" /> |  
| [Arduino Uno Datasheet (PDF)](https://docs.arduino.cc/resources/datasheets/A000066-datasheet.pdf) | [L298N Motor Driver Datasheet (PDF)](https://www.handsontec.com/dataspecs/L298N%20Motor%20Driver.pdf) | [GA12-N20 Motor Datasheet (PDF)](https://www.handsontec.com/dataspecs/GA12-N20.pdf) |  

### Circuit Diagram  

<img width="727" height="710" alt="image" src="https://github.com/user-attachments/assets/23582813-5fad-415e-802c-78c769814993" />  

---

### Tasks Completed  
1. **PI Controller Development**  
   - Implemented a Proportional–Integral (PI) controller in Simulink.  
   - Tuned controller gains to achieve stable velocity response.  

2. **Encoder Integration**  
   - Read encoder signals from **channels A and B**.  
   - Used the Simulink **Encoder block** to decode position data.  

3. **Velocity Estimation**  
   - Calculated motor velocity by taking the derivative of position data.  
   - Verified the real-time feedback loop between encoder and controller.  

4. **Velocity Control**  
   - Applied a **step function input** and validated PI tracking performance.  
   - Applied a **sine wave input** to observe dynamic controller behavior.  

---

### Submissions  

- **Simulink Model Screenshot:** Block diagram showing PI controller, encoder interface, and PWM driver.  
<img width="1920" height="1020" alt="Screenshot 2025-09-25 162112" src="https://github.com/user-attachments/assets/b04f4c7e-a3ae-4621-8379-4da885c21403" />


1. Step Input


https://github.com/user-attachments/assets/fb0b1a87-939e-4516-b65a-400038c26b57




| Test Case | Video Output | Velocity Plot |
|-----------|--------------|---------------|
| Step Input | https://github.com/user-attachments/assets/fb0b1a87-939e-4516-b65a-400038c26b57 | <img width="1076" height="864" alt="Screenshot 2025-09-25 162120" src="https://github.com/user-attachments/assets/09570b3f-fd86-4416-b37c-ab9dcac49998" /> |
| Sine Input | *(Add video link here)* | *(Add plot image here)* |

---
