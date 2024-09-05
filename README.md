# Monte_Carlo_Simultation
- This contains **C++** code of **Importance Sampling Monte Carlo Simulation** of **Canonical Ensemble**(N,V,T).
- We use **Lennard Jones Potential** with defined **sigma** and **epsilon** values and a **cutoff radius** for an **ideal monoatomic gas**.
- **Minimum image convention** and **periodic boundary conditions** were implemented to improve accuracy for bulk properties and simulate a model depicting real-life systems(representing an infinite system using a finite repeating system).
- Moves were accepted/rejected based on a **probabilistic criterion**(*Boltzmann factor*), ensuring that the calculation doesn’t get stuck in any particular **local minima** and continues to search for the **global minima**.
- The probability criterion samples unique configurations from the phase space, ensuring that the higher in energy we go, the less likely we are to accept that change, biasing towards the low-energy structures(as they would contribute more strongly to the average properties).
