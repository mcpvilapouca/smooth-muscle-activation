# UMAT for smooth muscle with calcium driven activation

🔸 This constitutive model was developed to simulate the behavior of smooth muscle activation due to the increase in the intracelular calcium concentration.

🔸 It was developed to simulate the uterine contrations during labor. This study is available in the literature, with more details regarding the example below.

🔸 The code is written in FORTRAN and the subroutine can be used with the finite element software ABAQUS. More details regarding the implementation and specific equations can be found in the literature [2].


###### [1] Vila Pouca, M.C.P., Ferreira, J.P.S., Oliveira, D.A. et al. Simulation of the uterine contractions and foetus expulsion using a chemo-mechanical constitutive model. Biomech Model Mechanobiol 18, 829–843 (2019). https://doi.org/10.1007/s10237-019-01117-5
###### [2] Sharifimajd B, Thore CJ, Stålhand J (2016) Simulating uterine contraction by using an electro-chemo-mechanical model. Biomech Model Mechanobiol 15:497–510. https ://doi.org/10.1007/s10237-015-0703-z

## Example - Simulation of the uterine contractions causing the fetus expulsion during labor

- This simulation was performed using ABAQUS/Standard (Implicit)
- The uterus was subjected to an activation curve, which represents the increase of calcium concentration with time.
    - We simulated the uterine contractions with 90 seconds followed by resting stages of 60 seconds
    - The fetal head is expulsed during the 5th uterine contraction
- The maternal pushes were not considered and the only load expulsing the fetus is the uterus contraction

![uterus_compressed](https://user-images.githubusercontent.com/95075305/171604627-b596d673-19e2-4838-98bd-ca1bf1b0f358.gif)
