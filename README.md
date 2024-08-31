# Monte_Carlo_Simultation
- This contains **C++** code of **Importance Sampling Monte Carlo Simulation** of **Canonical Ensemble**(N,V,T).
- We use **Lennard Jones Potential** with defined **sigma** and **epsilon** values and a **cutoff radius** for an **ideal monoatomic gas**.
- **Minimum image convention** and **periodic boundary conditions** were implemented to improve accuracy and bring the model closer to real-life systems.
- Moves were accepted/rejected based on a **probabilistic criterion**(*Boltzmann factor*) ensuring that calculation doesn’t get stuck in any **local minima** and continues until it reaches **global minima** sampling out from the phase space. 


## Further Calculations:
Pressure => We had taken an ideal gas therefore: PV = nRT
            n = N/Na
            Since N and T are already defined at the starting of the simulation, Na and R are constants
            Pressure can be easily calculated.

Heat Capacity => Heat required to raise the Temperature of the system by 1K(_ΔT_)
                Qc = (n * Cv * ΔT) (_Since Volume = Constant_)
                Now for an mono atomic gas Cv = 3/2 * R and n = N/Na
                Therefore plugging in these values the Heat Capacity can be easily obtained.
