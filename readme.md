# Simple 1D Particle in Cell
Particle in cell simulation with poisson solver written in [rust][rust-lang]. The mesh is computed with second-order centered finite difference, force interpolation is computed with Cloud in Cell, time integrator is computed with leapfrog (Stormer-Verlet) method. 

## Readings:
1. https://perswww.kuleuven.be/~u0052182/weather/pic.pdf
2. http://www.sciencedirect.com/science/article/pii/S0010465514001994
3. https://www.gnu.org/software/archimedes/manual/html/node29.html

## TODO:
- [x] - Implement Accleration Function
- [x] - Implement L2Norm Function
- [ ] - Implement Poisson Function
- [ ] - Implement Evolution 
- [ ] - Implement writing to file

[rust-lang]:https://www.rust-lang.org/
