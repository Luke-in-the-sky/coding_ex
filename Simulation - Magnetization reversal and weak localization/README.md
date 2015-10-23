Simulate the magnetization reversal curve for YIG under different Magnetic Field intensities

The core of the simulation is in the function test(float, float, int, double, bool). 
You can then loop the simulation over for different values of the simulation's parameters
by using the function ManyTimesTheMain(). Important parameters of the simulation are allocated
in a dedicated class named GlobalConstants